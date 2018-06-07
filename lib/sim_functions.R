require(dplyr)

parameters <- list(beta = 0.9,
                   beta_shape = 0.06,
                   zeta = 0.35,
                   gamma_shape = 3.36,
                   gamma_scale = 1.0, 
                   T = 16,
                   sim_T = 10,
                   I0 = c(1,1,0,0,0,0,0),
                   P = rep(500, 7))

outbreak_sim <- function(pars) {
  with(pars, {
     ## Get number of locations
    C <- length(I0)
    
    ## Generate daily inf distribution
    daily_inf <- rep(0,T)
    daily_inf[1] <- 0.99
    daily_inf[2] <- 0.01 
    
    ## Make matrix of susceptibles at each time step
    S <- matrix(0, nrow = C, ncol = T)
    I <- matrix(0, nrow = C, ncol = T)
    dimnames(I) <- list("group" = 1:length(I0),
                     "time" = 1:T)
    r <- matrix(0, nrow = C, ncol = T)
    dimnames(r) <- list("group" = 1:length(I0),
                     "time" = 1:T)
    
    ## Get initial susceptibles for each location
    S[,1] <- P - I0
    I[,1] <- I0
   
    ## In cell contact matrix
    in_cell_A <- matrix(0, nrow = C, ncol = C)
    diag(in_cell_A) <- 1/P
    
    ## Between cell contact matrix
    between_cell_A <- matrix(0, nrow = C, ncol = C)
    for (i in 1:C) {
      for (j in 1:C) {
        if (i != j) {
          between_cell_A[i,j] = 1/(sum(P)-P[i])
        }
      }
    }
    
    ## Make a matrix to calculate lambda at each timestep
    lambda <- matrix(0, nrow = C, ncol = T)
    
    for (i in 1:sim_T) {
      
      ## Sample FOI from distribution
      cell_r <- rgamma(C, scale = beta/beta_shape, shape = beta_shape*I[,i])
      r[,i] <- cell_r
     
      ## Propagate in-cell per-capita risk of infection across days
      in_cell_inf <- (1-zeta)*as.vector(in_cell_A %*% cell_r) 
      between_cell_inf <- zeta*as.vector(between_cell_A %*% cell_r) 
      total_daily_inf <- (in_cell_inf + between_cell_inf) %o% daily_inf 
      lambda[,i:T] <- lambda[,i:T] + total_daily_inf[,1:(T-i+1)]
      
      ## Draw Infections
      new_inf <- rbinom(C, S[,i], prob = 1.0 - exp(-lambda[,i]))
      
      ## Update susceptibles
      if (i < T) {
        S[,i+1] <- S[,i] - new_inf
      }
      for (j in 1:C) {
        n <- new_inf[j]
        if (n > 0) {
          ## Sample onset days
          camp_incidence <- i + rgeom(n, prob = 0.6) + 1
          if( sum(camp_incidence <= T) > 0) {
            camp_incidence <- table(camp_incidence[camp_incidence <= T])
            onset_days <- as.numeric(names(camp_incidence))
            I[j,onset_days] <- I[j,onset_days] + as.vector(camp_incidence) 
        }
      }

    }
    }
  
    sim_results <- list(I = I,
                        r = r) %>%
      sim_to_df
  return(list(sim = sim_results,
              contact_end = sim_T))
  })
}

sim_to_df <- function(sim) {
  I_df <- sim$I %>% 
    as.table %>% 
    as.data.frame %>%
    rename(I = Freq)
  
  r_df <- sim$r %>%
    as.table %>%
    as.data.frame %>%
    rename(r = Freq)
  
  sim_df <- r_df %>% inner_join(I_df, by = c("group","time"))
  sim_df$time <- as.numeric(sim_df$time)
  return(sim_df)
}

intervention_outbreak <- function(sim, c) {
  
  if (sum(sim$sim$I) < c) {
    return(FALSE)
  }
  ## First, filter out days without incident cases and 
  ## after the end time
  total_df <- sim$sim %>%
    group_by(time) %>%
    summarize(I = sum(I), 
              r = sum(r)) %>%
    mutate(total_I = cumsum(I)) %>%
    filter(time <= sim$contact_end)
  
  ## Get the time we cross the threshhold
  if (sum(total_df$total_I >= c) == 0) {
    return(FALSE)
  } 
  thresh_t <- min(which(total_df$total_I >= c))
  
  total_df$after <- as.numeric(total_df$time >= thresh_t)
  
  
  if (sum(total_df$after) == 0) {
    return(FALSE)
  }
  return(TRUE)
}





before_after_r <- function(sim_df, end_T, c) {
  
  ## First, filter out days without incident cases and 
  ## after the end time
  total_df <- sim_df %>%
    group_by(time) %>%
    summarize(I = sum(I), 
              r = sum(r)) %>%
    mutate(total_I = cumsum(I)) %>%
    filter(time <= end_T)
  
  ## Get the time we cross the threshhold
  thresh_t <- min(which(total_df$total_I >= c))
  
  total_df$after <- as.numeric(total_df$time >= thresh_t)
 
  before_after <- total_df %>%
    group_by(after) %>%
    filter(I > 0) %>%
    summarize(avg_r = sum(r)/sum(I))
  
  return(before_after)
  
}

sample_ratios <- function(pars, c, N = 100) {
  before_r <- rep(0, N)
  after_r <- rep(0, N)
  out_ratios <- rep(0, N)
  for (i in 1:N) {
    print(i)
  while(TRUE) {
    sample_outbreak <- outbreak_sim(pars)
    if (intervention_outbreak(sample_outbreak, c) == TRUE) {
      break
    }
  }
    
    br <- before_after_r(sample_outbreak$sim, 
                         sample_outbreak$contact_end,
                         c)
    before_r[i] <- br$avg_r[1]
    after_r[i] <- br$avg_r[2]
    out_ratios[i] <- br$avg_r[2]/br$avg_r[1]

  }
    
    return(data.frame(before = before_r,
                      after = after_r,
                      ratio = out_ratios))
  
}

sample_outbreaks <- function(pars, N = 500) {
  model_runs <- list()
  for (i in 1:N) {
    model_runs[[i]] <- outbreak_sim(pars)
  }
  
  return(model_runs)
}

sample_large_outbreak <- function(pars, n = 300, maxiter = 500) {
  
  for (i in 1:maxiter) {
    I <- gentime_sim(pars)
    if (sum(I) >= n) {
      return(I)
    }
  }
  
  return(NA)
  
}

outbreak_histogram <- function(pars, n = 300, maxiter = 500) {
  
  total_I <- rep(0,maxiter)
  for (i in 1:maxiter) {
    I <- outbreak_sim(pars)
    total_I[i] <- sum(I)
    }
  
  return(total_I)
  
}

sim_to_infer_input <- function(parameters, I) {
  
  T <- ncol(I)
  C <- nrow(I)
  
  sim_df <- data.frame()
  for (i in 1:C) {
    t <- 1:T
    y <- I[i,]
    camp_df <- data.frame(t = t, y = y)
    camp_df$c <- i
    
    sim_df <- rbind(sim_df, camp_df)
  }
  
  ## Filter out zero days
  sim_df <- sim_df %>% filter(y > 0)
  
    data_in = list(C = nrow(I),
                T = ncol(I),
                N = nrow(sim_df),
                P = parameters$P,
                J = sim_df$c,
                t = sim_df$t,
                Y = sim_df$y,
                epsilon = 0.6)

  return(data_in)
}
               

