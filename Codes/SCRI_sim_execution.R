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

source(here("Codes", "SCRI_main_functions.R"))
source(here("Codes", "SCRI_sim_workflow_functions.R"))
source(here("Codes", "SCRI_helper_functions.R"))

options(scipen = 999)

# Sample size tables

sample_size_tab <- sample_size_table(size = seq(1000, 10000, by = 1000))

# List of methods

methods <- c("no_calendar", "calendar_adjustment")

plan(multisession, workers = 40)
n_sim <- 1000
full_simulation(scenario_table = sample_size_tab, 
                n_sim = n_sim, 
                seeds = get_seeds(n_sim, scenario_table = sample_size_tab), 
                methods = methods, 
                output_dir = here("Results"))

power_results <- summarise_simulation_results(method_scen = method_scen(method_table = as.data.frame(methods),
                                                                       scenario_table = sample_size_tab),
                                             nsim = n_sim,
                                             results_dir = here("Results"),
                                             summary_dir = file.path(here("Results"), "Summary"),
                                             summary_file_name = "Summary_all_scens")
