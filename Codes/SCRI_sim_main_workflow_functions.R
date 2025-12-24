##########################################
### Project: SCRI design for           ###
### vaccine effectiveness              ###
##########################################

###################################
### Script 2: Workflow of       ###
### the simulation              ###
###################################


##############################
# 1 - Table of all scenarios ----
##############################

scenarios <- function(vacc_seasonality = c("uniform", "beta"),
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
                      cohort_size = c(6000, 10000, 20000)) {
  
  # Scenarios of misspecifying risk and control window, no seasonality of vaccination
  scen_misspecify <- expand.grid(
    scen = c("misspecify_control", "misspecify_risk_sta", "misspecify_risk_end"),
    sample_size = cohort_size,
    KEEP.OUT.ATTRS = FALSE
  ) %>%
    mutate(
      base_infect_shape = baseline_infection_shape[1],
      base_infect_mode  = baseline_infection_mode[1],
      base_infect_min = baseline_infection_min[1],
      base_infect_max = baseline_infection_max[1],
      vacc_season        = vacc_seasonality[1],
      vacc_mean_d       = NA_real_,
      vacc_sd_d         = NA_real_,
      control_start     = control_start_d,
      control_end       = c(control_end_d[1], rep(control_end_d[2], 2))[match(scen,
                                                                              c("misspecify_control","misspecify_risk_sta","misspecify_risk_end"))],
      risk_start        = c(risk_start_d[1], risk_start_d[2], risk_start_d[1])[match(scen,
                                                                                     c("misspecify_control","misspecify_risk_sta","misspecify_risk_end"))],
      risk_end          = c(risk_end_d[1], risk_end_d[1], risk_end_d[2])[match(scen,
                                                                               c("misspecify_control","misspecify_risk_sta","misspecify_risk_end"))]
    )
  
  # Scenarios of varying seasonality of baseline infection risk and vaccination date
  shape_mode_pair <- tibble(
    base_infect_shape = baseline_infection_shape,
    base_infect_mode  = baseline_infection_mode,
    base_infect_min = baseline_infection_min,
    base_infect_max = baseline_infection_max,
  )
  
  scen_season <- expand.grid(
    scen = "seasonality",
    base_infect_shape = baseline_infection_shape,
    vacc_mean_d = vacc_mean,
    vacc_sd_d = vacc_sd,
    sample_size = cohort_size,
    KEEP.OUT.ATTRS = FALSE
  ) %>%
    mutate(
      vacc_season    = vacc_seasonality[2],
      control_start = control_start_d,
      control_end   = control_end_d[2],
      risk_start    = risk_start_d[1],
      risk_end      = risk_end_d[1]
    ) %>%
    left_join(shape_mode_pair, by = "base_infect_shape") %>%
    relocate(c(base_infect_mode, base_infect_min, base_infect_max), .after = base_infect_shape)
  
  
  all_scen <- dplyr::bind_rows(scen_misspecify, scen_season)
  
  all_scen <- all_scen %>%
    mutate(
      scen_id   = row_number(),
      scen_name = paste0("S", scen_id, "_", scen, "_size", sample_size)
    )
  
  return(all_scen)
}


##############################
# 2 - Seed for each run ----
##############################

get_seeds <- function(n_sim, scenario_table){
  set.seed(20251115, kind = "L'Ecuyer-CMRG", sample.kind = "Rejection")
  n_seed <- nrow(scenario_table)*n_sim # Independent seed for each run in each scenario
  seed <- sample(1:1e9, 
                 size = n_seed, 
                 replace = FALSE)
}


############################## 
# 3 - Simulation study in layers ----
##############################

## --- Layer 1: Perform one run ------------------------------------------------

perform_one_run <- function(seed, rep, scen, methods, output_dir) {
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  
  # Generate data ---------
  tryCatch({
    data <- cohort_data(n_indiv = scen[['sample_size']], 
                        prop_vacc = 0.8, # proportion of vaccinated
                        vacc_dist = scen[['vacc_season']],
                         mean_vacc_date = scen[['vacc_mean_d']], 
                         sd_vacc_date = scen[['vacc_sd_d']],   
                         max_risk_period = max(scen[['risk_end']]) + 1, # the longest risk period to be considered in SCRI
                         n_days = 365,
                         peak_VE = 0.6,
                         day_start_immun = 8,
                         day_peak_immun = 16,
                         day_start_wane = 36,
                         day_end_immun = 150,
                         baseline_risk_vec = baseline_risk(n_days = n_days, 
                                                           gamma_shape = scen[['base_infect_shape']],
                                                           gamma_mode = scen[['base_infect_mode']], 
                                                           min_risk = scen[['base_infect_min']],
                                                           peak_risk = scen[['base_infect_max']]))

  }, error = function(e) {
    log_error(e, stage = "Data generation", seed = seed, scen_name = scen[['scen_name']], rep = rep)
    return(NULL) 
  })
  
  
  # Ensure directories exist once ----
  invisible(lapply(methods, function(a) {
    create_directory(file.path(output_dir, a))
  }))
  
  # Perform analysis and output results -----
  for (meth in methods) {
    
    file_path <- file.path(output_dir, meth, paste0(scen[['scen_name']], ".csv"))
    
    tryCatch({      
      out <- run_SCRI(dat = data, rep = rep, method = meth,
                      n_days = 365,
                      control_start = scen[['control_start']],
                      control_end   = scen[['control_end']],
                      risk_start    = scen[['risk_start']],
                      risk_end      = scen[['risk_end']],
                      start_calendar = NA,
                      calendar_interval = NA)
                                                           
    }, error = function(e) {
      log_error(e, stage = "Analysis",
                seed = seed, scen_name = scen[['scen_name']], rep = rep, method = meth)
    })
    
    out <- as.data.frame(cbind(out, seed = seed))
    
    append_to_csv(out, file_path)
    
    
  }
}

## --- Layer 2: Repeat 'Perform one run' n_sim times, for each scenario ---------
# (Parallelized)

sim_one_scenario <- function(scen, scenario_table, n_sim, seeds, methods, output_dir){
  
  message(paste0("Running scenario: ", scen[['scen_id']],"/",nrow(scenario_table)," ", scen[['scen_name']], ", at ", Sys.time()))
  
  scen_num <- as.numeric(scen['scen_id'])
  
  foreach(i = 1:n_sim,
          .options.future = list(packages = c("extraDistr","survival","SCCS","data.table","RhpcBLASctl"),
                                 seed = TRUE)) %dofuture% {
                                   
                                   # To avoid issues with native libraries, try keeping native libraries thread safe
                                   data.table::setDTthreads(1L)
                                   RhpcBLASctl::blas_set_num_threads(1L)
                                   RhpcBLASctl::omp_set_num_threads(1L)
                                   
                                   seed = seeds[(n_sim *(scen_num -1) + i)] #S1 gets seed nr 1-1000, S2 get seed nr 1001-2000 and so on
                                   perform_one_run(seed = seed, rep = i, scen = scen, methods = methods, output_dir = output_dir)
                                 }
}

## --- Layer 3: Perform the full simulation ------------------------------------
# For each scenario, generate n_sim datasets and perform the SCRI analysis

full_simulation <- function(scenario_table = sample_size_table(), 
                            n_sim = n_sim, 
                            seeds = get_seeds(n_sim = n_sim, scenario_table = sample_size_table()), 
                            methods = methods, 
                            output_dir = here("Results")) {
  message(paste("Simulation starts at:", Sys.time()))
  create_directory(output_dir)
  for (s in seq_len(nrow(scenario_table))){
    invisible(
      sim_one_scenario(scen = scenario_table[s,], 
                       scenario_table = scenario_table,
                       n_sim = n_sim, 
                       seeds = seeds, 
                       methods = methods, 
                       output_dir = output_dir)
    )
  }
  message(paste("Simulation ends at:", Sys.time()))
}


############################## 
# 4 - Summarise the simulation results ----
##############################

# Tables of all methods and scenarios

method_scen <- function(method_table = as.data.frame(methods), 
                        scenario_table){
  merged_table <- merge(method_table, scenario_table, by = NULL)
  return(merged_table)
}


# Function to read in the .csv simulation result files and summarise the bias

summarise_simulation_results <- function(method_scen = method_scen(),
                                         nsim = n_sim,
                                         true_VE = 0.6,
                                         results_dir = here("Results"),
                                         summary_dir = file.path(here("Results"), "Summary"),
                                         summary_file_name = "Summary_all_scens") {
  create_directory(summary_dir)
  
  all_summaries <- list()
  
  # Read in .csv files
  for (i in seq_len(nrow(method_scen))) {
    method_i <- method_scen$methods[i]
    scen_file <- paste0(method_scen$scen_name[i], ".csv")
    
    file_path <- file.path(results_dir, method_i, scen_file)
    
    if (!file.exists(file_path)) {
      warning(paste("File not found:", file_path))
      next
    }
    
    sim_data <- read.csv(file_path)
    if (is.null(sim_data)) next
    
    # Quantify bias
    summary_row <- tryCatch(
      summary_sim2(true_VE = true_VE, result_table = sim_data, n_sim = nsim),
      error = function(e) {
        warning(paste("Error applying summary_sim2 to", file_path, ":", e$message))
        return(NULL)
      }
    )
    if (is.null(summary_row)) next
    
    # Combine with scenario information
    summary_row$methods <- method_i
    summary_row$scen_name <- method_scen$scen_name[i]
    
    
    
    all_summaries[[length(all_summaries) + 1]] <- summary_row
  }
  
  # Combine all summaries
  summary_table <- do.call(rbind, all_summaries)
  if (!"methods" %in% names(summary_table)) summary_table$methods <- NA
  if (!"scen_name" %in% names(summary_table)) summary_table$scen_name <- NA
  summary_table_merge <- merge(method_scen, summary_table, 
                               by = c("methods", "scen_name"),
                               all.x = TRUE)
  
  # Export file
  output_path <- file.path(summary_dir, paste0(summary_file_name,".xlsx"))
  export(summary_table_merge, output_path)
  
  return(summary_table_merge)
}


