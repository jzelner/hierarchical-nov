random_incidence_df <- function() {
  ## Sample number of cases
  T <- 16
  end_T <- 10
  C <- 7
  Y <- rpois(C*T, lambda = 10)
  Y[Y==0] <- 1
  t <- rep(1:T, C)
  c <- rep(1:C, each = T)
  beta <- rgamma(C*T, Y*1, 1,)

  df <- data.frame(t = t,
                   Y = Y,
                   camp = c,
                   beta = beta)
  df$beta[df$t > end_T] <- 0

  weights <- matrix(0.25/100, C, C)
  diag(weights) <- 0.75/100
  return(list(
    df = df,
    weights = weights
  ))

}

sample_infectors <- function(incidence_df) {

}
