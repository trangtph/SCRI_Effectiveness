
##############################
# 1 - Table of all sample sizes to test 
##############################

sample_size_table <- function(size = seq(10000, 50000, by = 10000)) {
  
  sample_size_ <- expand.grid(size = size)
  
  sample_size_$scen_id <- seq(1:nrow(sample_size_))
  sample_size_$scen_name <- with(sample_size_, 
                              paste0("S",scen_id,"sample_size", size))
  return(sample_size_)
}


##############################
# 2 - Seed for each run 
##############################

get_seeds <- function(n_sim, scenario_table){
  set.seed(20251115, kind = "L'Ecuyer-CMRG", sample.kind = "Rejection")
  n_seed <- nrow(scenario_table)*n_sim # Independent seed for each run in each scenario
  seed <- sample(1:1e9, 
                 size = n_seed, 
                 replace = FALSE)
}


############################## 
# 3 - Simulation study in layers
##############################

# --- Layer 1: Perform one run ------------------------------------------------

perform_one_run <- function(seed, rep, scen, methods, output_dir) {
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  # generate data
  tryCatch({
    # --- Generate data ---
    data <- cohort_data(n_indiv = scen[['size']], 
                         prop_vacc = 0.8, # proportion of vaccinated
                         mean_vacc_date = 180, # Vacc date follows normal distribution
                         sd_vacc_date = 100,   # Vacc date follows normal distribution
                         max_risk_period = 78, # the longest risk period to be considered in SCRI
                         n_days = 365,
                         peak_VE = 0.6,
                         day_start_immun = 8,
                         day_peak_immun = 16,
                         day_start_wane = 36,
                         day_end_immun = 150,
                         baseline_risk_vec = baseline_risk())

  }, error = function(e) {
    log_error(e, stage = "Data generation", seed = seed, scen_name = scen[['scen_name']], rep = rep)
    return(NULL) 
  })
  
  
  # Ensure directories exist once
  invisible(lapply(methods, function(a) {
    create_directory(file.path(output_dir, a))
  }))
  
  # perform 2 analysis methods and output results 
  for (meth in methods) {
    
    file_path <- file.path(output_dir, meth, paste0(scen[['scen_name']], ".csv"))
    
    tryCatch({      
      out <- run_SCRI(dat = data, rep = rep, method = meth,
                      n_days = 365,
                      control_start = 3,
                      control_end   = 7,
                      risk_start    = 17,
                      risk_end      = 35,
                      start_calendar = 1,
                      calendar_interval = 30)
                                                           
    }, error = function(e) {
      log_error(e, stage = "Analysis",
                seed = seed, scen_name = scen[['scen_name']], rep = rep, method = meth)
    })
    
    out <- as.data.frame(cbind(out, seed = seed))
    
    append_to_csv(out, file_path)
    
    
  }
}

# --- Layer 2: Repeat 'Perform one run' n_sim times, for each scenario ---------
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

# --- Layer 3: Perform the full simulation ------------------------------------
# For each scenario, generate n_sim datasets and perform 6 analysis methods

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
# 4 - Summarise the simulation results
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
      summary_sim(true_VE = 0.6, result_table = sim_data, n_sim = nsim),
      error = function(e) {
        warning(paste("Error applying bias_quantification to", file_path, ":", e$message))
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
