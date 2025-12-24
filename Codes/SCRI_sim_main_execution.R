##########################################
### Project: SCRI design for           ###
### vaccine effectiveness              ###
##########################################

###################################
### Script 3: Execute           ###
### the simulation              ###
###################################


if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman")
}


library(pacman)
pacman::p_load(
  foreach,     # foreach loop
  stats,
  extraDistr,
  dplyr,
  SCCS,
  tictoc,      # Measure performance time
  survival,         # Export file  
  here,
  data.table,
  magrittr,     # To use the pipe %>%
  doRNG,        # Reproducible parallel session
  RhpcBLASctl,   # Control threads in parallel session
  future,
  doFuture,
  rio
)

source(here("Codes", "SCRI_sim_main_core_functions.R"))
source(here("Codes", "SCRI_sim_main_workflow_functions.R"))
source(here("Codes", "SCRI_helper_functions.R"))

options(scipen = 999)

##############################
# 1 - Base case scenarios ----
##############################

## Table of the scenarios' information
base_case <- data.frame(sample_size = c(6000, 10000, 20000))
base_case <-base_case %>% mutate(
  scen = "base_case",
  base_infect_shape = 2.5,
  base_infect_mode = 100,
  base_infect_min = 0.0002,
  base_infect_max = 0.002,
  vacc_season = "uniform",
  vacc_mean_d = NA,
  vacc_sd_d = NA, 
  control_start = 3,
  control_end = 7,
  risk_start = 16,
  risk_end = 35,
  scen_id = seq(1,3)
) %>% mutate(
  scen_name = paste0("S", scen_id, "_", scen, "_size", sample_size)
  
)

## Only fit SCRI model without ajustment for seasonality
methods <- c("no_calendar")

## Run the simulation
set.seed(20251218)
plan(multisession, workers = 40)
n_sim <- 1000
full_simulation(scenario_table = base_case, 
                n_sim = n_sim, 
                seeds = get_seeds(n_sim, scenario_table = base_case), 
                methods = methods, 
                output_dir = here("Results","Raw_results_base_case"))

## Summarize the results
base_case_results <- summarise_simulation_results(method_scen = method_scen(method_table = as.data.frame(methods),
                                                                            scenario_table = base_case),
                                                  nsim = n_sim,
                                                  true_VE = 0.6,
                                                  results_dir = here("Results","Raw_results_base_case"),
                                                  summary_dir = file.path(here("Results"), "Summary"),
                                                  summary_file_name = "Summary_base_case_20251219")

##############################
# 2 - Scenarios for bias quantification ----
##############################

## Table of all scenarios
scen_table <- scenarios(vacc_seasonality = c("uniform", "beta"),
                        vacc_mean = c(180, 80),
                        vacc_sd = c(60, 20),
                        baseline_infection_shape = c(2.5,10, 20),
                        baseline_infection_mode = c(100, 300, 200),
                        baseline_infection_min = c(0.0002, 0.0002, 0.0002),
                        baseline_infection_max = c(0.002, 0.003, 0.003),
                        control_start_d = c(3),
                        control_end_d = c(15, 7),
                        risk_start_d = c(16, 8),
                        risk_end_d = c(35,77),
                        cohort_size = c(6000, 10000, 20000))

## Only fit SCRI model without adjustment for seasonality
methods <- c("no_calendar")

## Run the simulation
set.seed(20251218)
plan(multisession, workers = 40)
n_sim <- 1000
full_simulation(scenario_table = scen_table, 
                n_sim = n_sim, 
                seeds = get_seeds(n_sim, scenario_table = scen_table), 
                methods = methods, 
                output_dir = here("Results","Raw_results_all_scens"))

## Summarize the results
results_all_scens <- summarise_simulation_results(method_scen = method_scen(method_table = as.data.frame(methods),
                                                                            scenario_table = scen_table),
                                                  nsim = n_sim,
                                                  true_VE = 0.6,
                                                  results_dir = here("Results","Raw_results_all_scens"),
                                                  summary_dir = file.path(here("Results"), "Summary"),
                                                  summary_file_name = "Summary_all_scens_20251219")

