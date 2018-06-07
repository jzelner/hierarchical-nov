#! /usr/bin/env Rscript
require(rstan)
data_in <- readRDS("output/nov_model_input.Rds")

m <- stan("src/contact_model.stan",
          data = data_in,
          iter = 2000,
          warmup = 1000,
          chains = 1, 
          control = list(adapt_delta = 0.90))

saveRDS(m, "output/nov_model.Rds")

