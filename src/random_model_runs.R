#! /usr/bin/env Rscript
require(drake)
require(dplyr)
pkgconfig::set_config("drake::strings_in_dots" = "literals")
source("lib/helper_functions.R")



plan <- drake_plan(
  outbreak_data = clean_outbreak_data(file_in("data/jamboree_data.csv")),
  pop_data = read_csv(file_in("data/jamboree_pop.csv")),
  camp_ids = make_camp_ids(pop_data),
  known_camp_data = make_baseline_dataset(outbreak_data,pop_data),
  unassigned_cases = filter_unassigned_cases(outbreak_data),
  random_data = sample_camps_for_unassigned(camp_ids,
                                            known_camp_data, 
                                            unassigned_cases, 
                                            unique(known_camp_data$camp)),
  data_in = list(C = nrow(camp_ids),
                T = length(unique(random_data$day)),
                N = nrow(random_data),
                P = pop_data$N,
                J = random_data$id,
                t = random_data$day,
                Y = random_data$cases,
                epsilon = 0.6
                ),
  
  m = stan(file_in("src/gentime.stan"),
          data = data_in,
          iter = 2000,
          chains = 1)
  

)
