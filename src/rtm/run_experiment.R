#! /usr/bin/env Rscript
pkgconfig::set_config("drake::strings_in_dots" = "literals")
require(drake)
require(dplyr)
require(ggplot2)
require(tidyr)
require(here)
source("lib/simple_outbreak.R")

depletion_drop <- function(r, r_eff) {
  
  r1 <- r %>% 
    select(r,k, p_cross)
  
  r_eff_1 <- r_eff %>%
    select(r, k, p_cross) %>%
    rename(p_cross_eff = p_cross) 
  
  rr <- inner_join(r1, r_eff_1) %>%
    mutate(drop_ar = 1 - (p_cross/p_cross_eff))
  
  
  return(rr)
  
}

plan <- drake_plan(
  int_thresh = 30,
  run_pars = experiment_parameters(
    r = seq(from = 0.5, to = 1.5, by = 0.025),
    k = c(0.01, seq(from = 0.05, to = 1.0, by = 0.05)),
    N = 1000,
    max_T = 14,
    c = 10
  ),
  
  
  exp_runs = run_experiment(run_pars, nsim = 2000),
  intervention_runs = run_intervention_experiment(run_pars, nsim = 1000),
  
  outbreak_prob = exp_runs %>% 
    group_by(parset, sim) %>%
    summarize(size = max(total_I)) %>%
    group_by(parset) %>%
    summarize(p = sum(size > 10)/n(),
              size = max(size)) %>%
    inner_join(run_pars) %>%
    dplyr::filter(r >= 0.7,
                  r <= 1.5,
                  k > 0.01),
  
  outbreaks = exp_runs %>%
    group_by(parset, sim) %>%
    summarize(outbreak = as.numeric(max(total_I) > run_pars$c)) %>%
    dplyr::filter(outbreak == 1) %>%
    select(-outbreak),

  avg_r_outbreak = outbreaks %>%
    inner_join(exp_runs) %>%
    dplyr::filter(I > 0) %>%
    group_by(parset, sim) %>%
    summarize(daily_mean_R = mean(R),
              case_mean_R = sum(lambda)/sum(I)) %>%
    group_by(parset) %>%
    summarize(case_mean_R = mean(case_mean_R)) %>%
    inner_join(run_pars),
  
  avg_r_all = exp_runs %>%
    dplyr::filter(I > 0) %>%
    group_by(parset) %>%
    summarize(daily_mean_R = mean(R),
              case_mean_R = sum(lambda)/sum(I)) %>%
   inner_join(run_pars),
    
  has_intervention = outbreaks %>%
    left_join(exp_runs) %>%
    mutate(intervention = as.numeric(total_I > run_pars$c)) %>%
    group_by(parset, sim) %>%
    summarize(int_outbreak = as.numeric(sum(intervention) > 0)),
  
  parset_rvals = has_intervention %>% 
    left_join(exp_runs) %>%
    dplyr::filter(I > 0) %>%
    mutate(intervention = as.numeric(total_I > run_pars$c)) %>%
    group_by(parset, sim, intervention) %>%
    summarize(avg_R = sum(lambda)/sum(I),
              avg_R_eff = sum(lambda*S/run_pars$N)/sum(I)) %>%
    mutate(intervention = if_else(intervention == 0, "R_pre", "R_post")),
  
  pre_post = parset_rvals %>%
    select(-avg_R_eff) %>% 
    spread(intervention, avg_R) %>%
    mutate(R_drop = R_post-R_pre,
           cross_one = as.numeric((R_post < 1) & (R_pre > 1))),


    
  avg_drop = pre_post %>%
    group_by(parset) %>%
    summarize(avg_pre = median(R_pre),
              avg_post = median(R_post),
              avg_drop = mean(R_drop),
              p_cross = mean(cross_one)) %>%
    inner_join(run_pars) %>% 
    dplyr::filter(r >= 0.7,
                  r <= 1.5,
                  k > 0.01),
  
  pre_post_reff = parset_rvals %>%
    select(-avg_R) %>% 
    spread(intervention, avg_R_eff) %>%
    mutate(R_drop = R_post-R_pre,
           cross_one = as.numeric((R_post < 1) & (R_pre > 1))),
 
 
  avg_drop_reff = pre_post_reff %>%
    group_by(parset) %>%
    summarize(avg_pre = median(R_pre),
              avg_post = median(R_post),
              avg_drop = mean(R_drop),
              p_cross = mean(cross_one)) %>%
    inner_join(run_pars) %>% 
    dplyr::filter(r >= 0.7,
                  r <= 1.5,
                  k > 0.01),
 
  dd = depletion_drop(avg_drop, avg_drop_reff),
  
   
  drop_g = ggplot(avg_drop, aes(x = k,
                                y = r,
                                fill = p_cross)) + 
    geom_tile() + 
    scale_fill_distiller(palette = "Spectral") + 
    geom_hline(yintercept = 1.0) +
    xlab("k") + 
    ylab("R0"),
  
  drop_g_reff = ggplot(avg_drop_reff, aes(x = k,
                                y = r,
                                fill = p_cross)) + 
    geom_tile() + 
    scale_fill_distiller(palette = "Spectral") + 
    geom_hline(yintercept = 1.0) +
    xlab("k") + 
    ylab("R"),
 
  drop_ar_g = ggplot(dd, aes(x = k,
                                y = r,
                                fill = drop_ar)) + 
    geom_tile() + 
    scale_fill_distiller(palette = "Spectral") + 
    geom_hline(yintercept = 1.0) +
    xlab("k") + 
    ylab("R0"),
 
    
  outbreak_g = ggplot(outbreak_prob, aes(x = k,
                                         y = r, 
                                         fill = p)) + 
    geom_tile() +
    scale_fill_distiller(type = "div",
                         palette = "Spectral") + 
    geom_hline(yintercept = 1.0) + 
    xlab("k") + 
    ylab("R0"),

  size_g = ggplot(outbreak_prob, aes(x = r,
                                         y = size, 
                                         colour = k)) + 
    geom_point() +
    scale_colour_distiller(type = "div",
                         palette = "Spectral") + 
    xlab("k") + 
    ylab("R0") + 
    theme_bw(),

  
    pre_g = ggplot(avg_drop, aes(x = k,
                                y = avg_pre,
                                colour = r)) + 
    geom_point() + 
    scale_colour_distiller(palette = "Spectral", direction = 1) + 
    xlab("k") + 
    ylab("Pre-intervention R") + 
    theme_bw(),
 
   post_g = ggplot(avg_drop, aes(x = k,
                                y = avg_post,
                                colour = r)) + 
    geom_point() + 
    scale_colour_distiller(palette = "Spectral", direction = 1) + 
    xlab("k") + 
    ylab("Post-intervention R") + 
    theme_bw(),
 
  rmarkdown::render(
    knitr_in(here("presentations/rtm.Rmd")),
    output_file = file_out(here("output/rtm.pdf"))
  )
  
)

make(plan)
