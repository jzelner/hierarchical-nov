#! /usr/bin/env Rscript
require(ggplot2)
require(dplyr)
require(readr)

camp_nums <- data.frame(camp = toupper(letters)[1:7],
                        campnum = as.character(1:7))


d <- read_csv("data/jamboree_data.csv") %>%
  left_join(camp_nums) %>%
  select(-camp) %>%
  rename(camp = campnum)

d$camp <- as.character(d$camp)



camp_d <- d
camp_d$camp[is.na(camp_d$camp)] <- "UNK"

camp_g <- ggplot(camp_d, aes(x=day, y = cases)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 3, linetype = "dashed") +
  geom_vline(xintercept = 10, linetype = "dotted") + 
  facet_grid( camp ~ .) +
  theme_bw() +
  xlab("Days since start of outbreak") +
  ylab("Incident cases cases in camp") + 
  scale_x_continuous(breaks = 0:15,
                   labels = 0:15)

ggsave("output/figures/cases_by_camp.pdf", width = 5, height = 8)

all_d <- camp_d %>%
  group_by(day) %>%
  summarize(cases = sum(cases))

all_g <- ggplot(all_d, aes(x = day, y = cases)) +
  geom_point() +
  geom_line() +
    theme_bw() +
  geom_vline(xintercept = 3, linetype = "dashed") +
  geom_vline(xintercept = 10, linetype = "dotted") +
  xlab("Days since start of outbreak") +
  ylab("Total cases") +
  scale_x_continuous(breaks = 0:15,
                   labels = 0:15)



ggsave("output/figures/all_camp_cases.pdf", width = 6, height = 3)
