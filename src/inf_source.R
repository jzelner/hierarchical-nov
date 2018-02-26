#! /usr/bin/env Rscript
require(dplyr)
require(rstan)
require(tidyr)

## Load traces
z <- readRDS("output/nov_model.Rds") 
z <- rstan::extract(z)

## Load input data

d <- readRDS("output/nov_model_input.Rds")

## Make weights for time lags
mu <- 1
sigma <- 1.5
max_T <- dim(z$lambda)[3]
num_camps <- dim(z$lambda)[2]
nsamp <- dim(z$lambda)[1]

dp_density <- plnorm(1:max_T, meanlog = mu, sdlog = sigma) - 
  plnorm(0:(max_T-1), meanlog = mu, sdlog = sigma)

dp_mat <- matrix(rep(dp_density, each = nsamp), nrow = nsamp)

all_p <- data.frame()
for (i in 1:num_camps) {
 
    for (t in 2:max_T) {
    lambda_c <- z$lambda[,i,1:(t-1)]
   
    dp_mat <- matrix(rep(rev(dp_density[1:(t-1)]), each = nsamp), nrow = nsamp)
    lambda_w <- dp_mat*lambda_c
    lw <- z$lambda_within[,i,1:(t-1)]
    inf_p <- data.frame(lw*lambda_w/rowSums(lambda_w))
    colnames(inf_p) <- 1:(t-1)
    inf_p <- inf_p %>% gather(key = "inf_day", value = "p") 
    inf_p$camp <- i
    inf_p$onset_day <- t
    all_p <- rbind(all_p, inf_p)
  }
}

## Make a df with the data
data_df <- data.frame(y = d$Y,
                      onset_day = d$t,
                      camp = d$J)

## Merge to get potential infection days
inf_df <- inner_join(all_p, data_df)
inf_df$inf_day <- as.numeric(as.character(inf_df$inf_day))

## Weighted inf day
zz <- inf_df %>% group_by(onset_day, camp) %>% summarize(pp = (sum(inf_day*p/1000)))