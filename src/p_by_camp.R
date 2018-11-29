#! /usr/bin/env Rscript

require(ggplot2)
require(dplyr)

z <- readRDS("output/infection_source.Rds")

## Load combined outputs
m <- readRDS("output/combined_values.Rds")

pw_vals <- data.frame(t(apply(m$p_within, 2, function(x) quantile(x, probs = c(0.05, 0.5, 0.95)))))

colnames(pw_vals) <- c("low_ci", "median", "high_ci")
pw_vals$camp <- 1:nrow(pw_vals)

zz <- inner_join(z$cases_by_camp, pw_vals)

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
