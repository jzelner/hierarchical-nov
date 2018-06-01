require(rstan)
require(readr)
require(dplyr)
require(purrr)


make_camp_ids <- function(pd) {
  camp_id <- data.frame(camp = unique(pd$camp)) %>%
    arrange(camp)
  camp_id$id <- 1:nrow(camp_id)
  camp_id$camp <- as.character(camp_id$camp)
  
  return(camp_id)
}

sample_camps_for_unassigned <- function(cdf, kdf, udf, camp_ids) {
  
  total_cases <- sum(udf$cases)
  out_days <- c()
  out_camps <- c() 


  for (i in 1:nrow(udf)) {
    
    out_days <- append(out_days, rep(udf$day[i], udf$cases[i]))
    out_camps <- append(out_camps, sample(camp_ids, udf$cases[i], replace = TRUE))
    
  }
  
  ## Count up cases on each day
  odf <- data.frame(camp = out_camps,
                    day = out_days) %>%
    dplyr::group_by(day, camp) %>%
    dplyr::summarize(cases = n()) %>%
    dplyr::mutate(camp = as.character((camp))) %>%
    dplyr::select(camp, day, cases)
  
  odf <- rbind(data.frame(odf), select(kdf, camp, day, cases)) %>%
    dplyr::group_by(camp, day) %>%
    dplyr::summarize(cases = sum(cases)) %>%
    dplyr::inner_join(cdf) %>%
    dplyr::filter(cases > 0)
  # kdf <- kdf %>% select(camp, day, cases)  
  # 
  #  odf <- rbind(kdf, odf) 
  # #    group_by(camp, day) %>%
  
  #   summarize(cases = sum(cases))
  # 
  return(odf)  
}

filter_unassigned_cases <- function(df) {
  out_df <- df %>%
    dplyr::filter(camp == "UNK")
  return(out_df)
  
}

clean_outbreak_data <- function(odf) {
  
  od <- read_csv(odf)
  od$day <- od$day + 1
  od$camp[is.na(od$camp)] <- "UNK"
  return(od)
  
}

make_baseline_dataset <- function(od, pd) {
## Make unique camp identifiers
camp_id <- data.frame(camp = unique(pd$camp)) %>%
                        arrange(camp)
camp_id$id <- 1:nrow(camp_id)
camp_id$camp <- as.character(camp_id$camp)

## Join with camp IDs
od_out <- od %>% 
  inner_join(camp_id) %>%
  dplyr::arrange(id)

  return(od_out)

}


randomize_assignments <- function(df) {
  
  
}