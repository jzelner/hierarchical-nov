#! /usr/bin/env Rscript
require(rstan)
require(dplyr)
require(ggplot2)

## Load data
z <- readRDS("output/nov_model.Rds") %>% extract

daily_r <- z$daily_avg_r %>%
  apply(2, function(x) quantile(x, probs = c(0.1, 0.5, 0.9))) %>%
  t %>%
  data.frame

colnames(daily_r) <- c("low_ci", "median", "high_ci")
daily_r$day <- 1:nrow(daily_r)
daily_r <- daily_r %>% select(day, low_ci, median, high_ci) %>%
  filter(day <= 11)

daily_r$day <- daily_r$day -1

g <- ggplot(daily_r) + geom_line(aes(x=day, y = median)) +
  geom_line(aes(x=day, y = low_ci), linetype = "dashed") +
  geom_line(aes(x=day, y = high_ci), linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_vline(xintercept = 3, linetype = "dashed") + 
  xlab("Outbreak day") +
  ylab("Avg. R across all infectious individuals") +
  theme_bw() +
  scale_x_continuous(breaks = min(daily_r$day):max(daily_r$day),
                     labels = min(daily_r$day):max(daily_r$day))

ggsave("output/figures/daily_avg_r.pdf", width = 6, height = 3)

num_camp <- dim(z$camp_r)[2]

all_camp_r <- data.frame()
for (i in 1:num_camp) { 
daily_r <- z$camp_r[,i,]%>%
  apply(2, function(x) quantile(x, probs = c(0.1, 0.5, 0.9))) %>%
  t %>%
  data.frame 

colnames(daily_r) <- c("low_ci", "median", "high_ci")
daily_r$day <- 1:nrow(daily_r)
daily_r$camp <- i
daily_r <- daily_r %>% select(camp,day, low_ci, median, high_ci) %>% filter(day <= 11)

all_camp_r <- rbind(all_camp_r, daily_r)

}

all_camp_r$day <- all_camp_r$day-1

camp_g <- ggplot(all_camp_r) + 
  geom_point(aes(x=day,y=median)) +
  geom_line(aes(x=day, y=median)) +
  geom_line(aes(x=day,y=high_ci),linetype = "dashed") +
  geom_line(aes(x=day,y=low_ci),linetype = "dashed") +
  geom_vline(xintercept = 3, linetype = "dotted") +
  facet_grid( camp ~ .) +
  theme_bw() +
  xlab("Day") +
  ylab("Daily reproduction number") +
  scale_x_continuous(breaks = min(daily_r$day):max(daily_r$day),
                     labels = min(daily_r$day):max(daily_r$day))

#  coord_cartesian(ylim=c(0,40)) 

ggsave("output/figures/daily_camp_r.pdf", width = 6, height = 8)
