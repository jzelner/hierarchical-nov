#! /usr/bin/env Rscript
require(drake)
require(dplyr)
pkgconfig::set_config("drake::strings_in_dots" = "literals")
source("lib/helper_functions.R")

generate_random_input_data <- function(camp_ids, known_camp_data, unassigned_cases, pop_data, end_T) {
  random_data = sample_camps_for_unassigned(camp_ids,
                                            known_camp_data, 
                                            unassigned_cases, 
                                            unique(known_camp_data$camp))
  ids = 1:nrow(random_data)
  data_in = list(C = nrow(camp_ids),
                 T = length(unique(random_data$day)),
                 end_T = 11,
                 before_end = sum(random_data$day <= 11),
                 before_end_id = ids[random_data$day <= 11],
                N = nrow(random_data),
                P = pop_data$N,
                J = random_data$id,
                t = random_data$day,
                Y = random_data$cases,
                epsilon = 0.8
                )

  return(data_in)
}

runmodel <- function(m, d) {
  ms = rstan::sampling(m, data = d,
                       iter = 2000,
                       chains = 1,
                       thin = 5,
                       control = list(
                         adapt_delta = 0.99,
                         max_treedepth = 15))
                       }
plan <- drake_plan(
  outbreak_data = clean_outbreak_data(file_in("data/jamboree_data.csv")),
  pop_data = read_csv(file_in("data/jamboree_pop.csv")),
  camp_ids = make_camp_ids(pop_data),
  known_camp_data = make_baseline_dataset(outbreak_data,pop_data),
  unassigned_cases = filter_unassigned_cases(outbreak_data),
  m = stan_model(file_in("src/contact_model.stan")),
)

random_d_plan <- drake_plan(
  data_in = generate_random_input_data(camp_ids, known_camp_data, unassigned_cases, pop_data, 11)
)

analysis_plan <- drake_plan(
  ms = runmodel(m, dataset__) 
)

random_d_fp <- expand_plan(random_d_plan, values = 1:20)

model_runs <- plan_analyses(analysis_plan, dataset = random_d_fp)
summary_list <- gather_plan(model_runs, target = "all_runs", gather = "list")

big_plan <- rbind(plan, random_d_fp, model_runs, summary_list)
make(big_plan, jobs = 6 )
loadd(all_runs)
saveRDS(all_runs, "output/random_model_runs.Rds")
