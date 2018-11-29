#! /usr/bin/env Rscript
require(rstan)
require(dplyr)
require(abind)

## Load model outputs
d <- readRDS("output/random_model_runs.Rds")

## Combine predictions across camps
samples <- list()
all_samples <- list()
zeta <- c()
beta_shape <- c()
log_beta_mu <- c()
gamma <- c()
p_within <- c()
for (i in 1:length(d)) {
  z <- extract(d[[i]])
  samples[[i]] <- extract(d[[i]])$camp_r
  all_samples[[i]] <- extract(d[[i]])$daily_avg_r
  zeta <- append(zeta,z$zeta)
  beta_shape <- append(beta_shape, z$beta_shape)
  log_beta_mu <- append(log_beta_mu, z$log_beta_mu)
  gamma <- append(gamma, z$gamma)
  p_within <- rbind(p_within, z$p_within)
}

camp_r_all <- abind(samples, along = 1)

## Combine average predictions
all_r_combined <- abind(all_samples, along = 1)

## combine into a list and write out
out_list <- list(
  camp_r = camp_r_all,
  daily_avg_r = all_r_combined,
  zeta = zeta,
  beta_shape = beta_shape,
  log_beta_mu = log_beta_mu,
  gamma = gamma,
  p_within = p_within
)

saveRDS(out_list, "output/combined_values.Rds")
