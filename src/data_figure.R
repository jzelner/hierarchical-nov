#! /usr/bin/env Rscript
require(ggplot2)
require(dplyr)
require(readr)

d <- read_csv("data/jamboree_data.csv")

camp_d <- d
camp_d$camp[is.na(camp_d$camp)] <- "UNK"

camp_g <- ggplot(camp_d, aes(x=day, y = cases)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 2, linetype = "dashed") +
  facet_grid( camp ~ .) +
  theme_bw() +
  xlab("Days since start of outbreak") +
  ylab("Cases in camp")

ggsave("output/figures/cases_by_camp.pdf", width = 6, height = 8)

all_d <- camp_d %>%
  group_by(day) %>%
  summarize(cases = sum(cases))

all_g <- ggplot(all_d, aes(x = day, y = cases)) +
  geom_point() +
  geom_line() +
    theme_bw() +
  geom_vline(xintercept = 2, linetype = "dashed") +
  xlab("Days since start of outbreak") +
  ylab("Total cases")



ggsave("output/figures/all_camp_cases.pdf", width = 6, height = 3)
