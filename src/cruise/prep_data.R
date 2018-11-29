require(yaml)
require(readr)
require(drake)
pkgconfig::set_config("drake::strings_in_dots" = "literals")

prep_stan_data <- function(md, od) {
  
  ## Calculate number of survivors
  S <- md$popsize - sum(od$cases)
  
  return(S)
  
}

exposure_matrix <- function(od) {
  
  ## Get outbreak duration
  outbreak_length <- nrow(od)
  
  ## Make an empty matrix
  contact_matrix <- matrix(0, outbreak_length, outbreak_length)
  
  ## Get duration of exposure
  indices <- as.matrix(expand.grid(1:outbreak_length, 1:outbreak_length))
  contact_matrix[indices] <- indices[,2] - indices[,1]
  contact_matrix[contact_matrix < 0] <- 0
  
  return(t(contact_matrix))
}

delay_weights <- function(d) {
  
  outbreak_length <- nrow(od)
}

plan <- drake_plan(md  = read_yaml(file_in("data/wikswo.yaml")),
                   od = read_csv(file_in("data/wikswo.csv")),
                   em = exposure_matrix(od),
                   sd = prep_stan_data(md, od))