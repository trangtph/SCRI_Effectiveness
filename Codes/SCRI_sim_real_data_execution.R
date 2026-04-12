##########################################
### Project: SCRI design for           ###
### vaccine effectiveness              ###
##########################################

###################################
### Script 3: Execute           ###
### the simulation for          ###
### real-world-based data       ###
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
source(here("Codes", "SCRI_sim_real_data.R"))


options(scipen = 999)

# ------------------------------------------------------------------------------
# 1. Table of all scenarios ----------------------------------------------------
# ------------------------------------------------------------------------------

scen_real <- expand.grid(cohort_size = c(20000, 50000, 100000, 500000),
                         risk_end = c(21,35))
scen_real <- scen_real %>%
  mutate(
    scen_id   = row_number(),
    scen_name = paste0("S", scen_id, "_riskend", risk_end, "_size", cohort_size)
  )

methods <- c("no_calendar", "calendar_30d", "calendar_7d", "calendar_7df3", "calendar_7df5") #calendar_adjustment

scen_real2 <- scen_real %>% filter(cohort_size == "500000")
methods2 <- "calendar_4d"

# ------------------------------------------------------------------------------
# 2. Perform the simulation ----------------------------------------------------
# ------------------------------------------------------------------------------

set.seed(20251218)
plan(multisession, workers = 40)
n_sim <- 1000

full_simulation_real(scenario_table = scen_real, 
                     n_sim = n_sim, 
                     seeds = get_seeds(n_sim, scenario_table = scen_real), 
                     methods = methods, 
                     output_dir = here("Results","Raw_results_real_20260319"))

set.seed(20251218)
plan(multisession, workers = 40)
n_sim <- 1000

full_simulation_real(scenario_table = scen_real2, 
                     n_sim = n_sim, 
                     seeds = get_seeds(n_sim, scenario_table = scen_real2), 
                     methods = methods2, 
                     output_dir = here("Results","Raw_results_real_3d_20260319"))

set.seed(20251218)
plan(multisession, workers = 40)
n_sim <- 1000

full_simulation_test(scenario_table = scen_real, 
                     n_sim = n_sim, 
                     seeds = get_seeds(n_sim, scenario_table = scen_real), 
                     methods = methods, 
                     output_dir = here("Results","Raw_results_test_20260319"))


# ------------------------------------------------------------------------------
# 3. Summarize the results ----------------------------------------------------
# ------------------------------------------------------------------------------


results_real <- summarise_simulation_results(method_scen = method_scen(method_table = as.data.frame(methods),
                                                                           scenario_table = scen_real),
                                                 nsim = n_sim,
                                                 true_VE = 0.67,
                                                 results_dir = here("Results","Raw_results_real_20260319_2"),
                                                 summary_dir = file.path(here("Results"), "Summary"),
                                                 summary_file_name = "Summary_real_20260319")

results_real_4d <- summarise_simulation_results(method_scen = method_scen(method_table = as.data.frame(methods2),
                                                                          scenario_table = scen_real2),
                                                nsim = n_sim,
                                                true_VE = 0.67,
                                                results_dir = here("Results","Raw_results_real_4d_20260319"),
                                                summary_dir = file.path(here("Results"), "Summary"),
                                                summary_file_name = "Summary_real_4d_20260319")

results_real_test <- summarise_simulation_results(method_scen = method_scen(method_table = as.data.frame(methods),
                                                                            scenario_table = scen_real),
                                                  nsim = n_sim,
                                                  true_VE = 0.67,
                                                  results_dir = here("Results","Raw_results_test_20260319"),
                                                  summary_dir = file.path(here("Results"), "Summary"),
                                                  summary_file_name = "Summary_real_test_20260319")

result_real <- rbind(results_real_all, results_real_3d)