require(dplyr)

parameters <- list(beta = 0.89,
                   beta_shape = 0.03,
                   zeta = 0.2,
                   gamma_shape = 3.36,
                   gamma_scale = 1.0, 
                   T = 16,
                   sim_T = 10,
                   I0 = c(1,1,0,0,0,0,0),
                   P = rep(500, 7))

gentime_sim <- function(pars) {
  with(pars, {
     ## Get number of locations
    C <- length(I0)
    
    ## Generate daily inf distribution
    daily_inf <- rep(0,T)
    daily_inf[1] <- 0.99
    daily_inf[2] <- 0.01 
              #pgamma(1:(T+1), shape = gamma_shape, scale = gamma_scale) %>% diff
    
    ## Make matrix of susceptibles at each time step
    S <- matrix(0, nrow = C, ncol = T)
    I <- matrix(0, nrow = C, ncol = T)
    
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
      cell_r <- rgamma(C, rate = beta*beta_shape, shape = beta_shape*I[,i])
      
      ## Propagate in-cell per-capita risk of infection across days
      in_cell_inf <- as.vector(in_cell_A %*% cell_r) 
      between_cell_inf <- as.vector(between_cell_A %*% (zeta*I[,i])) 
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
          camp_incidence <- i + ceiling(rlnorm(n))
          if( sum(camp_incidence <= T) > 0) {
            camp_incidence <- table(camp_incidence[camp_incidence <= T])
            onset_days <- as.numeric(names(camp_incidence))
            I[j,onset_days] <- I[j,onset_days] + as.vector(camp_incidence) 
        }
      }

    }
    }
  
  return(I)
  })
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
    I <- gentime_sim(pars)
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
               

