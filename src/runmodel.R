#! /usr/bin/env Rscript
require(rstan)
require(readr)
require(dplyr)

## Load outbreak data
od <- read_csv("data/jamboree_data.csv")

## Load pop'n data
pd <- read_csv("data/jamboree_pop.csv")

## Make unique camp identifiers
camp_id <- data.frame(camp = unique(pd$camp)) %>%
                        arrange(camp)
camp_id$id <- 1:nrow(camp_id)
camp_id$camp <- as.character(camp_id$camp)

## Join with camp IDs
od <- od %>% inner_join(camp_id) %>%
  mutate(t = day + 1) %>%
  arrange(id)

pd <- pd %>% inner_join(camp_id) %>%
  arrange(id)

data_in <- list(C = nrow(camp_id),
                T = length(unique(od$day)),
                N = nrow(od),
                P = pd$N,
                J = od$id,
                t = od$t,
                Y = od$cases
                )

m <- stan("src/model.stan",
          data = data_in,
          iter = 2000,
          chains = 1)

saveRDS(m, "output/nov_model.Rds")
