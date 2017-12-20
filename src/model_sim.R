#! /usr/bin/env Rscript
require(rstan)
require(dplyr)

outbreak_sim <- function(pop, I0, beta, alpha, gamma, T) {

  ## Pre-calculate the distribution of FOI over days of infection
  foi_dist <- pexp(0:T, gamma) %>% diff

  num_sites <- length(pop)
  ## Make an outcome matrix with the number of cases in each
  ## site
  I <- matrix(0, num_sites, T)
  I[,1] <- I0

  S <- matrix(0, num_sites, T)
  S[,1] <- pop - I0

  ## Make a matrix indicating the force of infection
  ## over time at each site
  lambda <- matrix(0, num_sites, T)

  for (t in 1:(T-1)) {
    ## Propagate risk from current step forward
    for (s in 1:num_sites) {
      ## Calculate within-site exposure
      lambda[s,t:T] <- lambda[s,t:T] + (beta[s]*I[s,t]*foi_dist[1:(T-t+1)]/pop[s])

      ## Calculate cross-site exposure
      lambda[s,t:T] <- lambda[s,t:T] + (alpha*(sum(I[,t])-I[s,t])*foi_dist[1:(T-t+1)])
    }

    ## Sample number of new infections in each location
    incidence <- rbinom(num_sites, S[,t], 1.0-exp(-lambda[,t]))

    ## Update susceptibles for following step
    S[,t+1] <- S[,t]-incidence

    ## Draw proportion of cases developing disease at t+1 and t+2
    num_next <- rbinom(num_sites, incidence, 0.5)
    I[,t+1] <- I[,t+1] + num_next
    if (t < (T-1)) {
      I[,t+2] <- I[,t+2] + incidence-num_next
    }
  }

  return(I)
}

## Load traces
z <- readRDS("output/nov_model.Rds") %>% extract

## Load data for initial conditions
d <- readRDS("output/nov_model_input.Rds")

## Reformat into df
df <- data.frame(J = d$J,
                 t = d$t,
                 Y = d$Y) %>% filter(t == 1)

## Simulate
nsamples <- dim(z$beta)[1]

total_cases <- rep(0, nsamples)
for (i in 1:nsamples) {
  zz <- outbreak_sim(d$P, df$Y, z$beta[i,], z$alpha[i], z$gamma[i], max(d$t))
  total_cases[i] <- sum(zz)
  }
