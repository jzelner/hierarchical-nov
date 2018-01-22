sir_outbreak <- function(beta, alpha, gamma, epsilon, C, N, I0 = 1, dt = 0.1, Z, T) {
  ## Set initial lambda values based on inputs
  init_inf_dens <- inf_dens(1)
  for (c1 in 1:C) {
    for (c2 in 1:C) {
      a <- 1
      if (c1 == c2) {
        a <- zeta
      }
      lambda[c2,] <- lambda[c2,] + (a*init_betas[c1]/sum(N))*I0[c1]*init_inf_dens
    }
  }


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
                      lambda = lambda,
                      incidence = incidence,
                      N = N,
                      S = N - I0)

  return(return_vals)

}

test_init <- initialize_outbreak(0.9, 0.5, 2.0, 0.7, 2.0, 7, rep(100, 7), c(1,0,0,0,0,0,0), c(2,2,0,0,0,0,0), 10, 14)
