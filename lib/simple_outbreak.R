require(dplyr)
require(purrr)
require(progress)

simple_pars <- list(r = 0.9,
                    k = 0.1,
                    max_T = 14,
                    N = 1000,
                    c = 10,
                    init_I = 1)

experiment_parameters <- function(r = 0.9, k = 0.1, max_T = 14, N = 1000, init_I = 1, c = 10) {
  
   pars <-  expand.grid(r = r, k = k, max_T = max_T, N = N, init_I = init_I, c = c)
   pars$parset <- 1:nrow(pars)
   return(pars)
  
}

run_experiment <- function(expars, nsim = 1000) {
  # pb <- progress_bar$new(
  #format = "  Running [:bar] :percent eta: :eta",
  #total = nrow(expars), clear = FALSE, width= 60)
  pb <- txtProgressBar(min = 1, max = nrow(expars), style = 3)
  experiments <- list()
  for (i in 1:nrow(expars)) {
    experiments[[i]] <- sample_n_outbreaks(expars[i,], nsim = nsim)
    experiments[[i]]$parset <- i
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  return(bind_rows(experiments))
  
}

run_intervention_experiment <- function(expars, nsim = 1000) {
  # pb <- progress_bar$new(
  #format = "  Running [:bar] :percent eta: :eta",
  #total = nrow(expars), clear = FALSE, width= 60)
  pb <- txtProgressBar(min = 1, max = nrow(expars), style = 3)
  experiments <- list()
  for (i in 1:nrow(expars)) {
    experiments[[i]] <- sample_n_intervention(expars[i,], nsim = nsim)
    experiments[[i]]$parset <- i
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  return(bind_rows(experiments))
  
}

homog_outbreak <- function(pars) {
  with(pars, {
    
    I <- rep(0, max_T)
    I[1] <- init_I
    S<- rep(0, max_T)
    S[1] <- N-I[1]
    lambda <- rep(0, max_T)
    for (i in 1:(max_T-1)) {
      lambda[i] <- rgamma(1, shape = k*I[i], scale = r/k)
      if (lambda[i] > 0) {
        I[i+1] <- rbinom(1, S[i], 1.0-exp(-lambda[i]/N))
        S[i+1] <- S[i] - I[i+1]
      }
    }
    df <- data.frame(I = I,
                     S = S,
                     total_I = cumsum(I),
                     lambda = lambda,
                     R = if_else(I > 0, lambda/I, 0))
    
    df$Reff <- df$R*S/N
    
    return(df)
  })
}

counterfactual_outbreaks <- function(pars, intervention_weights) {
  with(pars, {
    
    I <- matrix(0, nrow = max_T, ncol = length(intervention_weights))
    I[1,] <- init_I
    S<- matrix(0, nrow = max_T, ncol = length(intervention_weights))
    S[1,] <- N-I[1,]
    lambda <- matrix(0, nrow = max_T, ncol = length(intervention_weights))
    for (i in 1:(max_T-1)) {
      total_I <- sum(I[,1])
      if (total_I < c) {
        lambda[i,] <- rgamma(1, shape = k*I[i,1], scale = r/k)
        I[i+1,] <- rbinom(1, S[i,1], 1.0-exp(-lambda[i,1]/N))
        S[i+1,] <- S[i,] - I[i+1,]
      } else {
        lambda[i,] <- rgamma(ncol(I), shape = k*I[i,], scale = intervention_weights*r/k)
        I[i+1,] <- rbinom(ncol(I), S[i,], 1.0-exp(-lambda[i,]/N))
        S[i+1,] <- S[i,] - I[i+1,]       
      }
    }
    # df <- data.frame(I = I,
    #                  S = S,
    #                  total_I = cumsum(I),
    #                  lambda = lambda,
    #                  R = if_else(I > 0, lambda/I, 0))
    # 
    #df$Reff <- df$R*S/N
    colnames(I) <- intervention_weights
    return(I)
  })
}

intervention_protection <- function(pars) {
  
  intervention_weights <- seq(from = 1.0, to = 0.1, by = -0.1)
  I <- counterfactual_outbreaks(pars, intervention_weights)
  
  ## Get the number of cases in the no outbreak scenario
  total_I <- colSums(I)
  
  outbreak <- if_else(total_I[1] >= pars$c,1, 0)
  
  case_delta <- total_I[1] - total_I[2:length(total_I)]
  case_ratio <- total_I[2:length(total_I)]/total_I[1]
  
  out_df <- data.frame(
    cases_prevented = case_delta,
    case_ratio = case_ratio,
    intervention_protection = 1-intervention_weights[2:length(intervention_weights)]
  )
  
  out_df$outbreak <- outbreak
  
  return(out_df)
}


sample_n_intervention <- function(pars, nsim = 1000) {
 
   outbreaks <- list()
  for (i in 1:nsim) {
    ob <- intervention_protection(pars)
    ob$sim <- i
    outbreaks[[i]] <- ob
  }
  
  return(bind_rows(outbreaks))
 
}

sample_n_outbreaks <- function(pars, nsim = 1000) {
  
  outbreaks <- list()
  for (i in 1:nsim) {
    ob <- homog_outbreak(pars)
    ob$sim <- i
    outbreaks[[i]] <- ob
  }
  
  return(bind_rows(outbreaks))
  
}

sample_large_simple_outbreak <- function(pars, c = 30, maxiter = 10000) {
  
  for(i in 1:maxiter) {
    ob <- homog_outbreak(pars)
    if (sum(ob$I) >= c) {
      return(ob)
    }
  }
  
  return(FALSE)
  
}



simple_outbreak_sim <- function(pars) {

  
  with(pars, {
    
    ## Make daily infectiousness increments
    daily_inf <- pgamma(1:max_T, shape = gt_shape, scale = gt_mean/gt_shape)
    daily_inf <- c(daily_inf[1], diff(daily_inf))
    
    ## Create a matrix with lambda values
    lambda <- matrix(0, nrow = N, ncol = max_T)
    I <- matrix(0, nrow = N, ncol = max_T)
    
    ## Begin by drawing infection values for everone
    ind_r <- rgamma(N, shape = k, scale = r/k)
    
    onset_times <- rep(0, N)
    
    ## Create a shuffled vector of IDs we'll use 
    ## to select the next individual to be infected
    S_ID <- sample(1:N)
    
    incidence <- init_I
    inf_ids <- S_ID[1:incidence]
    I[inf_ids,1] <- 1
    onset_times[inf_ids] <- 1
    next_inf_index <- incidence + 1
    S <- N - incidence
    t <- 1
    inf_df <- data.frame()
    for (i_id in inf_ids) {
      lambda[i_id,(t+1):max_T] <- ind_r[i_id]*daily_inf[1:(max_T-t)] 
    }
    
    for (t in 2:max_T) {
      total_foi <- sum(lambda[,t])
      infectious_ids <- which(lambda[,t] > 0)
      inf_probs <- lambda[infectious_ids,t]/sum(lambda[infectious_ids,t])
      ## Draw new infections
      incidence <- rbinom(1, size = S, prob = 1.0 - exp(-total_foi/N))
      S <- S - incidence
      
      ## Steps to be completed only if there are new infections on a given step
      if (incidence > 0) {
        inf_ids <- S_ID[next_inf_index:(next_inf_index+incidence-1)]
        onset_times[inf_ids] <- t+1
        next_inf_index <- next_inf_index + incidence
        I[inf_ids,t] <- 1
        
       if (t < max_T) {

          for (i_id in inf_ids) {
            lambda[i_id,(t+1):max_T] <- ind_r[i_id]*daily_inf[1:(max_T-t)] 
          }
        }
      }
      
      ## Sample number of offspring
      offspring <- rmultinom(1, incidence, inf_probs)
      print("Offspring")
      print(offspring)
    }
    
    
    return(list(I = I,
                lambda = lambda
                ))
  })
  
}

sample_infectors <- function(model_out) {
  
  with(model_out,{
    max_T <- ncol(I)
    
    
    
  })
  
}
