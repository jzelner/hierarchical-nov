#! /usr/bin/env Rscript
require(rstan)
require(ggplot2)
require(dplyr)

zz <- readRDS("output/nov_model.Rds") %>% extract
df <- data.frame(x=as.vector(zz$beta))
g <- ggplot(df, aes(x, ..density..)) + 
  geom_histogram(binwidth=1) + 
  theme_bw() +
  xlab("Average daily effective reproduction number") +
  ylab("Posterior Density") + 
  scale_x_continuous(breaks = seq(from = 0, to = 30, by = 5), labels = seq(from = 0, to = 30, by = 5)) + 
  xlim(-1, 35)

ggsave("output/figures/marginal_inf_dist.pdf", width = 6, height = 4)
