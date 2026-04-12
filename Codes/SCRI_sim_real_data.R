##########################################
### Project: SCRI design for           ###
### vaccine effectiveness              ###
##########################################

###################################
### Script: Functions for       ###
### simulating and analyzing    ###
### real-world-based data       ###
###################################

# ------------------------------------------------------------------------------
# 1. Load packages and data ----------------------------------------------------
# ------------------------------------------------------------------------------

library(pacman)
pacman::p_load(
  here,
  rio,
  SCCS,
  survival, 
  data.table,
  dplyr,
  foreach,
  RhpcBLASctl,   # Control threads in parallel session
  future,
  doFuture
)

# Import daily baseline risk of COVID-19 infection
covid_risk_us_22 <- import(here("Data", "covid_risk_us_22.csv"))

covid_risk_base <- covid_risk_us_22$daily_risk_adj

# Import daily probability of COVID-19 vaccination 
covid_vacc_us_22 <- import(here("Data", "covid_vacc_us_22.csv"))

# Import daily vaccination probability

ve_full <- import(here("Data", "ve_daily.csv"))

ve_lookup_dat <- setNames(ve_full$ve, ve_full$day)

# ------------------------------------------------------------------------------
# 2. Function to generate cohort data ------------------------------------------
# ------------------------------------------------------------------------------

cohort_data_real <- function(n_indiv = 10000, 
                        prop_vacc = 0.8, # proportion of vaccinated
                        max_risk_period = 14, # the longest risk period to be considered in SCRI
                        n_days = 365,
                        baseline_risk_vec = covid_risk_base,
                        vacc_prob_dat = covid_vacc_us_22,
                        ve_lookup = ve_lookup_dat
                        ) {
  # Generate exposure information -------------------
  # Number vaccinated
  n_vacc <- round(n_indiv * prop_vacc)
  
  # Sample vaccination dates for vaccinated individuals
  # Trim the upper end of vaccination date so that the risk window is always within the observation period
  vacc_prob_lookup <- vacc_prob_dat[vacc_prob_dat$day <= n_days - max_risk_period,]

    vacc_date <- sample(
      x = vacc_prob_lookup$day,
      size = n_vacc,
      replace = TRUE,
      prob = vacc_prob_lookup$prob
    )
    vacc_date = c(vacc_date, rep(NA, n_indiv - n_vacc)) # add NA vacc_date for unvaccinated individuals

  
  # Generate daily outcome status -------------------
  
  infection <- matrix(0L, nrow = n_indiv, ncol = n_days)
  colnames(infection) <- paste0("day_", 1:ncol(infection))
  
  for (d in 1:n_days) {
    days_since_vacc <- ifelse(is.na(vacc_date), NA, d - vacc_date) #vector with n_indiv values
    VE <- numeric(length(days_since_vacc)) #vector with n_indiv values, default = 0
    
    # identify valid days to look up VE (days outside this range are assigned VE = 0)
    valid <- !is.na(days_since_vacc) &
      days_since_vacc >= 0 &
      days_since_vacc <= as.numeric(tail(names(ve_lookup_dat), 1))
    
    # assign new VE for valid days
    VE[valid] <-
      as.numeric(ve_lookup[as.character(days_since_vacc[valid])])
    prob_infect_d <- pmin(baseline_risk_vec[d] * (1 - VE), 1) #vector with n_indiv values
    infection[, d] <- stats::rbinom(n_indiv, 1, prob_infect_d)
    
  }
  
  # Output data, wide data format
  cohort_dat <- data.table::as.data.table(cbind(id = seq(1:nrow(infection)),
                                                vacc_date,
                                                infection))
  return(cohort_dat)
}

# ------------------------------------------------------------------------------
# 3. Simulation study in layers ------------------------------------------------
# ------------------------------------------------------------------------------

## --- Layer 1: Perform one run ------------------------------------------------

perform_one_run_real <- function(seed, rep, scen, methods, output_dir) {
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  
  # Generate data ---------
  tryCatch({
    data <- cohort_data_real(n_indiv = scen[['cohort_size']], 
                        prop_vacc = 0.8, 
                        max_risk_period = 14,
                        n_days = 365,
                        baseline_risk_vec = covid_risk_base,
                        vacc_prob_dat = covid_vacc_us_22,
                        ve_lookup = ve_lookup_dat
                        )
    
  }, error = function(e) {
    log_error(e, stage = "Data generation", seed = seed, scen_name = scen[['scen_name']], rep = rep)
    return(NULL) 
  })
  if (is.null(data)) return(invisible(NULL))
  
  
  # Ensure directories exist once ----
  invisible(lapply(methods, function(a) {
    create_directory(file.path(output_dir, a))
  }))
  
  # Perform analysis and output results -----
  
  for (meth in methods) {
    
    file_path <- file.path(output_dir, meth, paste0(scen[['scen_name']], ".csv"))
    
    out <- tryCatch({
      
      if (meth == "no_calendar") {
        run_SCRI(
          dat = data, rep = rep, method = "no_calendar",
          n_days = 365,
          control_start = 3,
          control_end   = 7,
          risk_start    = 15,
          risk_end      = scen[['risk_end']],
          start_calendar = NA,
          calendar_interval = NA
        )
      } else if (meth == "calendar_30d") {
        run_SCRI(
          dat = data, rep = rep, method = "calendar_adjustment",
          n_days = 365,
          control_start = 3,
          control_end   = 7,
          risk_start    = 15,
          risk_end      = scen[['risk_end']],
          start_calendar = 31,
          calendar_interval = 30
        )
        
      } else if (meth == "calendar_7d") {
        run_SCRI(
          dat = data, rep = rep, method = "calendar_adjustment",
          n_days = 365,
          control_start = 3,
          control_end   = 7,
          risk_start    = 15,
          risk_end      = scen[['risk_end']],
          start_calendar = 8,
          calendar_interval = 7
        )
        
      } else if (meth == "calendar_7df3") {
        run_SCRI(
          dat = data, rep = rep, method = "calendar_adjustment",
          n_days = 365,
          control_start = 3,
          control_end   = 7,
          risk_start    = 15,
          risk_end      = scen[['risk_end']],
          start_calendar = 3,
          calendar_interval = 7
        )
        
      } else if (meth == "calendar_7df5") {
        run_SCRI(
          dat = data, rep = rep, method = "calendar_adjustment",
          n_days = 365,
          control_start = 3,
          control_end   = 7,
          risk_start    = 15,
          risk_end      = scen[['risk_end']],
          start_calendar = 5,
          calendar_interval = 7
        )
        
      } else if (meth == "calendar_3d") {
        run_SCRI(
          dat = data, rep = rep, method = "calendar_adjustment",
          n_days = 365,
          control_start = 3,
          control_end   = 7,
          risk_start    = 15,
          risk_end      = scen[['risk_end']],
          start_calendar = 4,
          calendar_interval = 3
        )
        
      } else {
        stop("Unknown method: ", meth)
      }
      
    }, error = function(e) {
      log_error(e, stage = "Analysis",
                seed = seed, scen_name = scen[['scen_name']],
                rep = rep, method = meth)
      return(NULL)
    })
    
    if (is.null(out)) next
    
    out <- as.data.frame(cbind(out, seed = seed))
    
    append_to_csv(out, file_path)
  }
  invisible(NULL)
  
}

## --- Layer 2: Repeat 'Perform one run' n_sim times, for each scenario ---------
# (Parallelized)

sim_one_scenario_real <- function(scen, scenario_table, n_sim, seeds, methods, output_dir){
  
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
                                   perform_one_run_real(seed = seed, rep = i, scen = scen, methods = methods, output_dir = output_dir)
                                 }
}

## --- Layer 3: Perform the full simulation ------------------------------------
# For each scenario, generate n_sim datasets and perform the SCRI analysis

full_simulation_real <- function(scenario_table = scen_real, 
                            n_sim = n_sim, 
                            seeds = get_seeds(n_sim = n_sim, scenario_table = scen_real), 
                            methods = methods, 
                            output_dir = here("Results", "Raw_results_real")) {
  message(paste("Simulation starts at:", Sys.time()))
  create_directory(output_dir)
  for (s in seq_len(nrow(scenario_table))){
    invisible(
      sim_one_scenario_real(scen = scenario_table[s,], 
                       scenario_table = scenario_table,
                       n_sim = n_sim, 
                       seeds = seeds, 
                       methods = methods, 
                       output_dir = output_dir)
    )
  }
  message(paste("Simulation ends at:", Sys.time()))
}

# ------------------------------------------------------------------------------
# 4. Test the scenario of no seasonality of infection --------------------------
# ------------------------------------------------------------------------------

covid_risk_const <- rep(0.0001, length = length(covid_risk_base))

## --- Layer 1: Perform one run ------------------------------------------------

perform_one_run_test <- function(seed, rep, scen, methods, output_dir) {
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  
  # Generate data ---------
  tryCatch({
    data <- cohort_data_real(n_indiv = scen[['cohort_size']], 
                             prop_vacc = 0.8, 
                             max_risk_period = 14,
                             n_days = 365,
                             baseline_risk_vec = covid_risk_const,
                             vacc_prob_dat = covid_vacc_us_22,
                             ve_lookup = ve_lookup_dat
    )
    
  }, error = function(e) {
    log_error(e, stage = "Data generation", seed = seed, scen_name = scen[['scen_name']], rep = rep)
    return(NULL) 
  })
  if (is.null(data)) return(invisible(NULL))
  
  
  # Ensure directories exist once ----
  invisible(lapply(methods, function(a) {
    create_directory(file.path(output_dir, a))
  }))
  
  # Perform analysis and output results -----
  
  for (meth in methods) {
    
    file_path <- file.path(output_dir, meth, paste0(scen[['scen_name']], ".csv"))
    
    out <- tryCatch({
      
      if (meth == "no_calendar") {
        run_SCRI(
          dat = data, rep = rep, method = "no_calendar",
          n_days = 365,
          control_start = 3,
          control_end   = 7,
          risk_start    = 15,
          risk_end      = scen[['risk_end']],
          start_calendar = NA,
          calendar_interval = NA
        )
      } else if (meth == "calendar_30d") {
        run_SCRI(
          dat = data, rep = rep, method = "calendar_adjustment",
          n_days = 365,
          control_start = 3,
          control_end   = 7,
          risk_start    = 15,
          risk_end      = scen[['risk_end']],
          start_calendar = 31,
          calendar_interval = 30
        )
        
      } else if (meth == "calendar_7d") {
        run_SCRI(
          dat = data, rep = rep, method = "calendar_adjustment",
          n_days = 365,
          control_start = 3,
          control_end   = 7,
          risk_start    = 15,
          risk_end      = scen[['risk_end']],
          start_calendar = 8,
          calendar_interval = 7
        )
        
      } else if (meth == "calendar_7df3") {
        run_SCRI(
          dat = data, rep = rep, method = "calendar_adjustment",
          n_days = 365,
          control_start = 3,
          control_end   = 7,
          risk_start    = 15,
          risk_end      = scen[['risk_end']],
          start_calendar = 3,
          calendar_interval = 7
        )
        
      } else if (meth == "calendar_7df5") {
        run_SCRI(
          dat = data, rep = rep, method = "calendar_adjustment",
          n_days = 365,
          control_start = 3,
          control_end   = 7,
          risk_start    = 15,
          risk_end      = scen[['risk_end']],
          start_calendar = 5,
          calendar_interval = 7
        )
        
      } else if (meth == "calendar_4d") {
        run_SCRI(
          dat = data, rep = rep, method = "calendar_adjustment",
          n_days = 365,
          control_start = 3,
          control_end   = 7,
          risk_start    = 15,
          risk_end      = scen[['risk_end']],
          start_calendar = 5,
          calendar_interval = 4
        )
        
      } else {
        stop("Unknown method: ", meth)
      }
      
    }, error = function(e) {
      log_error(e, stage = "Analysis",
                seed = seed, scen_name = scen[['scen_name']],
                rep = rep, method = meth)
      return(NULL)
    })
    
    if (is.null(out)) next
    
    out <- as.data.frame(cbind(out, seed = seed))
    
    append_to_csv(out, file_path)
  }
  invisible(NULL)
  
}

## --- Layer 2: Repeat 'Perform one run' n_sim times, for each scenario ---------
# (Parallelized)

sim_one_scenario_test <- function(scen, scenario_table, n_sim, seeds, methods, output_dir){
  
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
                                   perform_one_run_test(seed = seed, rep = i, scen = scen, methods = methods, output_dir = output_dir)
                                 }
}

## --- Layer 3: Perform the full simulation ------------------------------------
# For each scenario, generate n_sim datasets and perform the SCRI analysis

full_simulation_test <- function(scenario_table = scen_real, 
                                 n_sim = n_sim, 
                                 seeds = get_seeds(n_sim = n_sim, scenario_table = scen_real), 
                                 methods = methods, 
                                 output_dir = here("Results", "Raw_results_real")) {
  message(paste("Simulation starts at:", Sys.time()))
  create_directory(output_dir)
  for (s in seq_len(nrow(scenario_table))){
    invisible(
      sim_one_scenario_test(scen = scenario_table[s,], 
                            scenario_table = scenario_table,
                            n_sim = n_sim, 
                            seeds = seeds, 
                            methods = methods, 
                            output_dir = output_dir)
    )
  }
  message(paste("Simulation ends at:", Sys.time()))
}
