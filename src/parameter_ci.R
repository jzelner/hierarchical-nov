#! /usr/bin/env Rscript
require(rstan)
require(readr)
require(dplyr)
require(ggplot2)

z <- readRDS("output/combined_values.Rds") 

analysis_pars <- c("zeta", "log_beta_mu", "gamma", "beta_shape")

all_pars <- data.frame()
for (p in analysis_pars) {
  df <- data.frame(x = z[[p]])
  df$iter <- 1:nrow(df)
  df$par <- p
  all_pars <- rbind(all_pars, df)
}
par_ci <- all_pars %>%
  group_by(par) %>%
  summarize(median = median(x),
            high_ci = quantile(x, probs = c(0.95)),
            low_ci = quantile(x, probs = c(0.05))
            )

write_csv(par_ci, "output/scalar_pars.csv")


## prob_par_ci <- par_ci %>% filter(par != "log_beta_mu")
## g <- ggplot(prob_par_ci, aes(x=par, y = median, ymin = low_ci, ymax = high_ci)) +
##   geom_point() +
##   geom_errorbar()
