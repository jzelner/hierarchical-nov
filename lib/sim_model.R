

ind_contribution <- function(beta, zeta, t, c, C, inf_d) {

  ind_inf_density <- inf_d(t)

  all_c <- matrix(0, nrow = C, ncol = length(ind_inf_density))
  for (i in 1:C) {
    all_c[i,] <- beta*ind_inf_density*ifelse(i == c, zeta, 1)
  }

  return(all_c)

}

ind_lambda_factory <- function(Z, gamma) {
  ## Make CDF of infection densities
  inf_dens <- dgeom(0:(Z-1), prob = gamma)

  ind_lambda <- function(t) {
    lam <- rep(0, Z)
    if (t > Z) {
      return(lam)
    }
    lam[t:Z] <- inf_dens[1:(Z-t+1)]
    return(lam)
  }

  return(ind_lambda)

}

onset_time_factory <- function(epsilon) {

  onset_times <- function(t) {
    return(t+1+rgeom(length(t), p = epsilon))
  }

  return(onset_times)

}


initialize_outbreak <- function(gamma, epsilon, zeta, beta_shape, beta_mean, C, N, I0, init_betas, Z, T) {

  require(dplyr)

  inf_dens = ind_lambda_factory(Z, gamma)
  onset_time = onset_time_factory(epsilon)

  ## Matrix of daily betas
  beta <- rgamma(Z*C, shape = beta_shape, scale = beta_mean/beta_shape) %>%
    matrix(ncol=Z, nrow=C)

  ## Set init betas based on inputs
  beta[,1] <- init_betas

  ## Make a matrix of FOIs
  lambda <- matrix(0, C, Z)

  ## Make a matrix of incidences
  incidence <- matrix(0, C, T)

  ## Set the first column of the incidence to be the
  ## number of initial infected
  incidence[,1] <- I0

  return_vals <- list(inf_dens = inf_dens,
                      onset_time = onset_time,
                      beta = beta,
                      zeta = zeta,
                      lambda = lambda,
                      incidence = incidence,
                      N = N,
                      S = N - I0,
                      Z = Z,
                      T = T,
                      C = C)

  return(return_vals)

}

seir_outbreak <- function(pars) {

  incidence <- pars$incidence
  lambda <- pars$lambda
  total_N <- sum(pars$N)
  S <- pars$S
  C <- pars$C
  T <- pars$T
  for (t in 1:pars$Z) {

    ## Calculate FOI from each camp at time T
    for (i in 1:C) {
      camp_incidence <- incidence[i,t]
      if (camp_incidence > 0) {
        lambda <- lambda + (camp_incidence/total_N)*ind_contribution(pars$beta[i,t],
                                                                     pars$zeta,
                                                                     t,
                                                                     i,
                                                                     pars$C,
                                                                     pars$inf_dens)

      }

    
      }
    ## Draw number of new infections in each camp
    new_infections <- rbinom(C, S, 1.0-exp(-lambda[,t]))
    print("UPDATE")
    print(S)
 
    S <- S - new_infections
   print(new_infections)
    print(S)
      for (i in 1:C) {
        if (new_infections[i] > 0) {
          ot <- pars$onset_time(rep(t, new_infections[i]))
          for (j in 1:length(ot)) {
            if (ot[j] <= T) {
              incidence[i,ot[j]] <- incidence[i,ot[j]] + 1
              }
          }
        }
      }




  }

  return(incidence)
}

test_init <- initialize_outbreak(0.9, 0.9, 2.0, 1.0, 2.0, 7, rep(100, 7), c(1,0,0,0,0,0,0), c(2,2,0,0,0,0,0), 10, 14)

test_run <- seir_outbreak(test_init)

