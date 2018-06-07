require(purrr)
require(dplyr)

mu <- 1
shape <- 1

branch_model <- function(mu, shape) {
  gen_size <- c()
  i <- 0
  scale <- mu/shape
  while(TRUE) {
    
    if (i == 0) {
      total_shape <- shape
    } else {
      total_shape <- gen_size[i]*shape
    }
    
    num_offspring <- rpois(1, rgamma(1, shape = total_shape, scale = scale))
    gen_size <- append(gen_size, num_offspring)
    if(num_offspring == 0) {
      return(gen_size)
    } 
  i <- i + 1
  }
}

sample_runs <- list()
for (i in 1:1000) {
  sample_runs[[i]] <- branch_model(1,0.1)
}

final_sizes <- map_dbl(sample_runs, sum)
num_gen <- map_dbl(sample_runs, length)

## Find the ones with >= 30 cases
large_outbreaks <- sample_runs[which(final_sizes >= 30)]


## Get generation when exceeded thresshold
peak_gen <- function(x, c) {
  cumulative_size <- cumsum(x)
  return(min(which(cumulative_size >= c)[1]))
}

gen_peaks <- map_dbl(large_outbreaks, function(x) peak_gen(x, 30))


## Get generation when exceeded thresshold
pre_post <- function(x, c) {
  cumulative_size <- cumsum(x)
  peak_gen <- min(which(cumulative_size >= c)[1])
  gen_ratios <- x[2:length(x)]/x[1:(length(x)-1)]
  pre_r <- mean(gen_ratios[1:(peak_gen-1)])
  post_r <- mean(gen_ratios[peak_gen:length(gen_ratios)])
  return(post_r/pre_r)
}

pp <- map_dbl(large_outbreaks, function(x) pre_post(x, 50))
ratio <- sum(pp < 1)
