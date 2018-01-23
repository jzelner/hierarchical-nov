#! /usr/bin/env Rscript
require(dplyr)
require(rstan)
require(tidyr)
require(ggplot2)
source("lib/sim_model.R")

## Load input data
input_d <- readRDS("output/nov_model_input.Rds")

## Subset out initial vals
df <- data.frame(Y = input_d$Y,
                 J = input_d$J,
                 t = input_d$t) %>% filter(t == 1)

## Load sims
z <- readRDS("output/nov_model.Rds") %>% extract

## Get num_samples
nsamples <- dim(z$beta)[1]
i <- 1


output_sims <- matrix(0, nrow = nsamples, ncol = input_d$T)

for (i in 1:nsamples) {

  test_init <- initialize_outbreak(z$gamma[i], input_d$epsilon, z$zeta[i], z$beta_shape[i], z$beta_shape[i]/exp(z$log_beta[i,1]), input_d$C, input_d$P, df$Y, z$beta[i,,1], 13, input_d$T) 


  test_run <- seir_outbreak(test_init)
  output_sims[i, ] <- colSums(test_run)

}

qvals <- apply(output_sims,2,function(x) quantile(x,probs = c(0.1, 0.25, 0.5, 0.6, 0.7, 0.8, 0.9))) %>% t

qv_df <- data.frame(qvals)
qv_df$t <- 1:nrow(qv_df)
colnames(qv_df) <- c("q10", "q25", "q50", "q60", "q70", "q80", "q90", "t")
qv_df <- qv_df %>% gather(key = "quantile", value =cases, -t)
g <- ggplot(qv_df, aes(x=t, y=cases)) + geom_line(aes(colour = as.factor(quantile)))
