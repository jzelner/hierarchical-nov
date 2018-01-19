#! /usr/bin/env Rscript
require(rstan)
require(dplyr)
require(ggplot2)

z <- readRDS("output/homogeneous_nov_model.Rds") %>% extract

day_incr <- apply(z$beta, 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975))) %>% t %>% data.frame

colnames(day_incr) <- c("low_ci", "median", "high_ci")
day_incr$day <- 1:nrow(day_incr)


g <- ggplot(day_incr, aes(x=day)) +
  geom_line(aes(y=median)) +
  geom_line(aes(y=low_ci), linetype = "dashed") +
  geom_line(aes(y=high_ci), linetype = "dashed")+
  geom_hline(yintercept = 1, linetype = "dotted")

ggsave("output/homogeneous_beta.pdf", width = 5, height = 3)
