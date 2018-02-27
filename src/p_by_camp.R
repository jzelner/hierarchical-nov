#! /usr/bin/env Rscript

require(ggplot2)
require(dplyr)

z <- readRDS("output/infection_source.Rds")
zz <- inner_join(z$cases_by_camp, z$by_camp)


g <- ggplot(zz, aes(x=p_total, 
                    y = median, 
                    ymin = low_ci, 
                    ymax = high_ci, 
                    label = camp)) + 
  geom_point() + 
  geom_text(nudge_x = 0.005, size = 4) + 
  geom_errorbar() + 
  theme_bw() +
  xlab("Proportion of incident cases") +
  ylab("Proportion of infections required within camp") + 
  scale_x_continuous(breaks = seq(0.0, 0.3, by = 0.05), labels = seq(0.0, 0.3, by = 0.05))

ggsave("output/figures/p_by_camp.pdf", width = 6, height = 4)