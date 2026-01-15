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
source(here("Codes", "SCRI_sim_samplesize_workflow_functions.R"))
source(here("Codes", "SCRI_helper_functions.R"))

options(scipen = 999)

# Sample size tables

sample_size_tab <- sample_size_table()

# List of methods

methods <- c("no_calendar")

plan(multisession, workers = 2)
n_sim <- 1000
full_simulation(scenario_table = sample_size_tab[4:nrow(sample_size_tab),], 
                n_sim = n_sim, 
                seeds = get_seeds(n_sim, scenario_table = sample_size_tab), 
                methods = methods, 
                output_dir = here("Results","Raw_sample_size_20260106"))

power_results <- summarise_simulation_results(method_scen = method_scen(method_table = as.data.frame(methods),
                                                                       scenario_table = sample_size_tab),
                                             nsim = n_sim,
                                             results_dir = here("Results","Raw_sample_size_20260106"),
                                             summary_dir = file.path(here("Results"), "Summary"),
                                             summary_file_name = "Power_results_20260106")
