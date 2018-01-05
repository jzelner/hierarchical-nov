#! /usr/bin/env Rscript
require(rstan)
require(dplyr)
require(ggplot2)

## Load data
z <- readRDS("output/nov_model.Rds") %>% extract

daily_r <- z$daily_avg_r %>%
  apply(2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975))) %>%
  t %>%
  data.frame

colnames(daily_r) <- c("low_ci", "median", "high_ci")
daily_r$day <- 1:nrow(daily_r)
daily_r <- daily_r %>% select(day, low_ci, median, high_ci)

g <- ggplot(daily_r) + geom_line(aes(x=day, y = median)) +
  geom_line(aes(x=day, y = low_ci), linetype = "dashed") +
  geom_line(aes(x=day, y = high_ci), linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dotted")

num_camp <- dim(z$camp_r)[2]

all_camp_r <- data.frame
for (i in 1:num_camp) {
daily_r <- z$camp_r[,i,]%>%
  apply(2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975))) %>%
  t %>%
  data.frame

colnames(daily_r) <- c("low_ci", "median", "high_ci")
daily_r$day <- 1:nrow(daily_r)
daily_r$camp <- i
daily_r <- daily_r %>% select(camp,day, low_ci, median, high_ci)

all_camp_r <- rbind(all_camp_r, daily_r)

}
