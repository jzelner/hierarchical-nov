#! /usr/bin/env Rscript
require(rstan)
data_in <- readRDS("output/nov_model_input.Rds")

m <- stan("src/model.stan",
          data = data_in,
          iter = 2000,
          chains = 1)

saveRDS(m, "output/nov_model.Rds")

