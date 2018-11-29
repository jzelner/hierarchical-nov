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
## ADD ITER`` 
    for (t in 2:max_T) {
    lambda_c <- z$lambda[,i,1:(t-1)]
   
    dp_mat <- matrix(rep(rev(dp_density[1:(t-1)]), each = nsamp), nrow = nsamp)
    lambda_w <- dp_mat*lambda_c
    lw <- z$lambda_within[,i,1:(t-1)]
    inf_p <- data.frame(rowSums(lw*lambda_w/rowSums(lambda_w)))
    colnames(inf_p) <- c("p_within")
    inf_p$camp <- i
    inf_p$onset_day <- t
    inf_p$iter <- 1:nrow(inf_p)
    all_p <- rbind(all_p, inf_p)
  }
}

## Make a df with the data
data_df <- data.frame(y = d$Y,
                      onset_day = d$t,
                      camp = d$J)

## Get the total Y on days > 1
total_incident_cases <- sum(data_df$y[data_df$onset_day > 1])

## Merge to get potential infection days
inf_df <- inner_join(all_p, data_df)

## Get estimate of total proportion of cases attributable to 
## within-camp transmission
within_camp_p <- inf_df %>%
  group_by(iter) %>%
  summarize(p_within = sum(y*p_within)/total_incident_cases)

## Break estimates out by camp
cases_by_camp <- data_df %>%
  filter(onset_day > 1) %>% 
  group_by(camp) %>% 
  summarize(total_y = sum(y))

cases_by_camp$p_total <- cases_by_camp$total_y/sum(cases_by_camp$total_y)

within_by_camp <- inf_df %>% 
  group_by(iter, camp) %>%
  summarize(total_within = sum(y*p_within)) %>%
  inner_join(cases_by_camp) %>%
  mutate(p_within = total_within/total_y) %>%
  group_by(camp) %>%
  summarize(median = median(p_within),
            low_ci = quantile(p_within, probs = c(0.025)),
            high_ci = quantile(p_within, probs = c(0.975)))

## Break out by day
cases_by_day <- data_df %>%
  filter(onset_day > 1) %>% 
  group_by(onset_day) %>% 
  summarize(total_y = sum(y))

within_by_day <- inf_df %>% 
  group_by(iter, onset_day) %>%
  summarize(total_within = sum(y*p_within)) %>%
  inner_join(cases_by_day) %>%
  mutate(p_within = total_within/total_y) %>%
  group_by(onset_day) %>%
  summarize(median = median(p_within),
            low_ci = quantile(p_within, probs = c(0.025)),
            high_ci = quantile(p_within, probs = c(0.975)))


data_out <- list(within_camp = within_camp_p,
                 by_camp = within_by_camp,
                 cases_by_camp = cases_by_camp)

saveRDS(data_out, "output/infection_source.Rds")
