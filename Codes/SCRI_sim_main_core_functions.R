##########################################
### Project: SCRI design for           ###
### vaccine effectiveness              ###
##########################################

###################################
### Script 1: Functions for     ###
### simulating and analyzing    ###
### data.                       ###
###################################

# ------------------------------------------------------------------------------
# 1. Function for daily baseline infection risk --------------------------------
# ------------------------------------------------------------------------------
# This function creates daily baseline infection risk for n_days = 365 days. 
# The daily baseline infection risk follows a gamma distribution, peaking at day `gamma_mode`.

baseline_risk <- function(n_days = 365, # number of follow-up days for cohort
                          gamma_shape = 2.5,
                          gamma_mode = 100, # Day of peak infection risk
                          min_risk = 2e-4,
                          peak_risk = 2e-3){
  # scale parameter from mode
  scale <- (gamma_mode - 1) / (gamma_shape - 1)
  quantile <- seq(0, n_days - 1, length.out = n_days)
  pdf_vals <- dgamma(quantile, shape = gamma_shape, scale = scale)
  # scale pdf into [min_risk, peak_risk]
  baseline_risk <- min_risk + (peak_risk - min_risk) * (pdf_vals / max(pdf_vals))
  return(baseline_risk)
}

# ------------------------------------------------------------------------------
# 2. Function for vaccination date----------------------------------------------
# ------------------------------------------------------------------------------
# This function creates vaccination date between days [low, up] for vaccinated individuals. 
# The vaccination date follows a beta distribution with user-specified mean and sd.

sample_beta_scaled <- function(n, mean, sd, low = 1, up = 365) {
  # convert mean/sd on [low, up] to Beta(a,b) parameters
  m <- (mean - low) / (up - low)               # scaled mean in [0,1]
  v <- (sd^2) / ((up - low)^2)               # scaled variance
  
  # Solve for a and b
  tmp <- m * (1 - m) / v - 1
  a <- m * tmp
  b <- (1 - m) * tmp
  
  # sample Beta(a, b) on [0,1]
  x <- rbeta(n, a, b)
  # rescale to [low, up]
  y <- low + x * (up - low)
  return(y)
}


# ------------------------------------------------------------------------------
# 3. Function to calculate daily vaccine effectiveness, day 0-150 post-vaccination 
# ------------------------------------------------------------------------------

# For vaccinated individuals, we assume zero immunity from days 0-7, 
# partial immunity increasing linearly from days 8-16, 
# and full protection at 60% VE during the risk interval period 
# before it begins to wane linearly from days 36- 150.
# Day 0 is day of vaccination.

get_VE <- function(days_since_vacc,
                   peak_VE = 0.6, # highest vaccine effectiveness
                   day_start_immun = 8, # day after vacc that immunity starts to increase
                   day_peak_immun = 16, # day after vacc that immunity reaches the peak
                   day_start_wane = 36, # day after vacc that immunity starts to wane
                   day_end_immun = 150) { # day after vacc that immunity is back to 0
  VE <- numeric(length(days_since_vacc))
  # Before vaccination or unvaccinated
  VE[days_since_vacc < 0 | is.na(days_since_vacc)] <- 0
  
  # From day 0 to [day_start_immun - 1] : no protection
  VE[!is.na(days_since_vacc) &
       days_since_vacc >= 0 & days_since_vacc < day_start_immun] <- 0
  
  # From day_start_immun to day_peak_immun: linear rise to peak_VE
  rise_idx <- !is.na(days_since_vacc) &
    days_since_vacc >= day_start_immun &
    days_since_vacc <= day_peak_immun
  VE[rise_idx] <- peak_VE * (days_since_vacc[rise_idx] - day_start_immun + 1) / (day_peak_immun - day_start_immun + 1)
  
  # From  day_peak_immun to [day_start_wane -1]: plateau at peak_VE
  plateau_idx <- !is.na(days_since_vacc) &
    days_since_vacc > day_peak_immun &
    days_since_vacc < day_start_wane
  VE[plateau_idx] <- peak_VE
  
  # From day_start_wane to day_end_immun: linear waning to 0
  wane_idx <- !is.na(days_since_vacc) &
    days_since_vacc >= day_start_wane &
    days_since_vacc <= day_end_immun
  VE[wane_idx] <- peak_VE * (1 - (days_since_vacc[wane_idx] - day_start_wane+1) / (day_end_immun - day_start_wane+1))
  
  # Beyond day_end_immun: no protection
  VE[!is.na(days_since_vacc) & days_since_vacc > day_end_immun] <- 0
  return(VE)
}

# ------------------------------------------------------------------------------
# 4. Function to generate cohort data ------------------------------------------
# ------------------------------------------------------------------------------

# This function generate data for a cohort of n_indiv individuals observed for 
# n_days (=365 days). Data include vaccination status (assuming 1 dose of vaccine, 
# 80% vaccinated), and daily infection status.

cohort_data <- function(n_indiv = 10000, 
                        prop_vacc = 0.8, # proportion of vaccinated
                        vacc_dist = c("uniform", "beta"), # beta for seasonality of vaccination, uniform for no seasonality
                        mean_vacc_date = 180, # specify this if vacc date follows beta distribution
                        sd_vacc_date = 100,   # specify this if vacc date follows beta distribution
                        max_risk_period = 78, # the longest risk period to be considered in SCRI
                        n_days = 365,
                        peak_VE = 0.6,
                        day_start_immun = 8,
                        day_peak_immun = 16,
                        day_start_wane = 36,
                        day_end_immun = 150,
                        baseline_risk_vec = baseline_risk()) {
  # Generate exposure information -------------------
  # Number vaccinated
  n_vacc <- round(n_indiv * prop_vacc)
  
  # Sample vaccination dates for vaccinated individuals 
  # Case 1: no seasonality of vaccination date - uniform distribution fo vaccination date
  if (vacc_dist == "uniform") {
    vacc_date <- round(stats::runif(n = n_vacc, 
                              min = 1, 
                              max = n_days - max_risk_period))
    vacc_date = c(vacc_date, rep(NA, n_indiv - n_vacc)) # add NA vacc_date for unvaccinated individuals
  }
  # Case 2: vaccination date has seasonality - (truncated) beta distribution 
  else {
    vacc_date <- round(sample_beta_scaled(n = n_vacc, 
                                          mean = mean_vacc_date, 
                                          sd = sd_vacc_date, 
                                          low = 1, 
                                          up = n_days - max_risk_period))
    vacc_date = c(vacc_date, rep(NA, n_indiv - n_vacc)) # add NA vacc_date for unvaccinated individuals
  }
  
  # Generate daily outcome status -------------------
  
  infection <- matrix(0L, nrow = n_indiv, ncol = n_days)
  colnames(infection) <- paste0("day_", 1:ncol(infection))
  
  for (d in 1:n_days) {
    days_since_vacc <- ifelse(is.na(vacc_date), NA, d - vacc_date) #vector with n_indiv values
    VE <- get_VE(days_since_vacc, 
                 peak_VE = peak_VE,
                 day_start_immun = day_start_immun,
                 day_peak_immun = day_peak_immun,
                 day_start_wane = day_start_wane,
                 day_end_immun = day_end_immun) #vector with n_indiv values
    
    prob_infect_d <- baseline_risk_vec[d] * (1 - VE) #vector with n_indiv values
    infection[, d] <- stats::rbinom(n_indiv, 1, prob_infect_d)
    
  }

  # Output data, wide data format
  cohort_dat <- data.table::as.data.table(cbind(id = seq(1:nrow(infection)),
                                                vacc_date,
                                                infection))
  return(cohort_dat)
}

# ------------------------------------------------------------------------------
# 5. Function to reshape data to SCRI-compatible format and fit SCRI model---------
# ------------------------------------------------------------------------------

# When the SCRI model does not include seasonality, the model only includes 
# vaccination status. The data of vaccinated cases
# will be reshaped to long format with information of outcome status in the 
# control window, risk window and washout window. 

# When the SCRI model controls for seasonality, the model include vaccination status
# and calendar time (divided into {365/calendar_interval + 1} intervals).


run_SCRI <- function(dat,
                     method,
                     n_days = 365,
                     rep = 1,
                     control_start = 3,
                     control_end   = 7,
                     risk_start    = 17,
                     risk_end      = 35,
                     start_calendar = 31, #Start of the second interval, the 1st interval start at the start of observation time
                     calendar_interval = 30) {

  infection_cols <- grep("^day_", names(dat), value = TRUE)
  #Keep vaccinated individuals with ≥1 infection
  dat[, infection_count := rowSums(.SD), .SDcols = infection_cols]
  vaccinated_infected <- dat[!is.na(vacc_date) & infection_count > 0]
  
  # Reshape to long format
  long <- data.table::melt(
    vaccinated_infected,
    id.vars = c("id", "vacc_date"),
    measure.vars = infection_cols,
    variable.name = "day",
    value.name   = "infection"
  )
  
  long[, day_infection := as.integer(sub("day_", "", day))]
  
  # Keep only days with event
  long <- long[infection == 1, .(id, vacc_date, day_infection)]
  
  # Define risk and control period
  long[, `:=`(
    control_start_d = vacc_date + control_start,
    control_end_d   = vacc_date + control_end,
    risk_start_d    = vacc_date + risk_start,
    risk_end_d      = vacc_date + risk_end
  )]
  
  # Keep only infection events occurring in control or risk periods
  long <- long[
    (day_infection >= control_start_d & day_infection <= control_end_d) |
      (day_infection >= risk_start_d  & day_infection <= risk_end_d)
  ]
  
  if (nrow(long) == 0) {
    return(data.frame(
      est_V = NA, se_V = NA, IRR_V = NA,
      IRR_V_CI_Lower = NA, IRR_V_CI_Upper = NA,
      VE = NA, VE_CI_Lower = NA, VE_CI_Upper = NA,
      n_event = 0, p_val = NA
    ))
  }
  
  if (method == "no_calendar"){
    # SCRI data format
    long2 <- SCCS::formatdata(
      indiv  = id,
      astart = control_start_d,
      aend   = risk_end_d,
      aevent = day_infection,
      adrug  = vacc_date,
      aedrug = control_end_d,
      expogrp = control_start,
      washout = c(1, risk_start - 1 - control_end),
      data = long
    )
    # Relevel exposure so control period = reference
    long2$vacc_date <- relevel(factor(long2$vacc_date), ref = "1") # Because it is coded as: 0 = risk period, 1 = control period, 2 = wash-out period
    
    # SCRI model
    mod <- summary(
      clogit(event ~ vacc_date + strata(indivL) + offset(log(interval)),
             data = long2))
    
    # Extract estimates and return as data frame
    
    est_V   <- mod$coefficients[1,1]
    se_V    <- mod$coefficients[1,3]
    
    IRR_V   <- exp(est_V)
    IRR_L   <- mod$conf.int[1,3]
    IRR_U   <- mod$conf.int[1,4]
    
    VE      <- 1 - IRR_V
    VE_L    <- 1 - IRR_U
    VE_U    <- 1 - IRR_L
    
    n_event <- mod$nevent
    p_val   <- mod$coefficients[1,5]
    
    return(data.frame(
      rep = rep,
      est_V, se_V,
      IRR_V, IRR_V_CI_Lower = IRR_L, IRR_V_CI_Upper = IRR_U,
      VE, VE_CI_Lower = VE_L, VE_CI_Upper = VE_U,
      n_event, p_val,
      n_calendar_group = NA,
      row.names = NULL
    ))
    
  }
  
  else if (method == "calendar_adjustment"){
    
    calendar_time_group <- seq(start_calendar, n_days-(calendar_interval-1), by = calendar_interval) 
    min_start <- min(long$control_start_d, na.rm = TRUE)
    max_end <- max(long$risk_end_d, na.rm = TRUE)
    
    # Cut the calendar_time_group so that they lie within the observation period
    calendar_time_group_cut <- calendar_time_group[
      calendar_time_group >= min_start+(calendar_interval-1) & calendar_time_group <= max_end -(calendar_interval-1)] # So that the first and last group is adequately long 
    
    long2 <- SCCS::formatdata(indiv = id,
               astart = control_start_d,
               aend = risk_end_d,
               aevent = day_infection,
               adrug = vacc_date,
               aedrug = control_end_d,
               expogrp =control_start,
               washout = c(1,risk_start-1-control_end),
               agegrp = calendar_time_group_cut,
               data=long)
    # Relevel exposure so control period = reference
    long2$vacc_date <- relevel(factor(long2$vacc_date), ref = "1") # Because it is coded as: 0 = risk period, 1 = control period, 2 = wash-out period
    
    # Fit the Conditional Poisson model
    # If there is only one "age" group: do not add "age" term in the model
    n_age <- nlevels(factor(long2$age))
    base_formula <- event ~ vacc_date + strata(indivL) + offset(log(interval))
    
    if (n_age > 1) {
      form <- update(base_formula, . ~ . + age)
    } else {
      form <- base_formula
    }
    
    mod <- summary(clogit(form, data = long2))
    
    # Extract estimates and return as data frame
    
    est_V   <- mod$coefficients[1,1]
    se_V    <- mod$coefficients[1,3]
    
    IRR_V   <- exp(est_V)
    IRR_L   <- mod$conf.int[1,3]
    IRR_U   <- mod$conf.int[1,4]
    
    VE      <- 1 - IRR_V
    VE_L    <- 1 - IRR_U
    VE_U    <- 1 - IRR_L
    
    n_event <- mod$nevent
    p_val   <- mod$coefficients[1,5]
    n_calendar_group <- n_age
    
    return(data.frame(
      rep = rep,
      est_V, se_V,
      IRR_V, IRR_V_CI_Lower = IRR_L, IRR_V_CI_Upper = IRR_U,
      VE, VE_CI_Lower = VE_L, VE_CI_Upper = VE_U,
      n_event, p_val,
      n_calendar_group,
      row.names = NULL
    ))
  }
  
  else {
    stop("Error: Method should be either no_calendar, or calendar_adjustment")
  }

}

# ------------------------------------------------------------------------------
# 6. Function to summarize results (for power calculation) ---------------------
# ------------------------------------------------------------------------------

summary_sim <- function(true_VE, result_table, n_sim=1000)
{

  # Number of missing values of estimated beta1 (e.g due to convergence)
  missing_estimate <- sum(is.na(result_table$est_V)) + n_sim - nrow(result_table)
  convergence_issue <- sum(result_table$IRR_V > 50, na.rm = TRUE )
  
  result_table2 <- filter(result_table, IRR_V < 50)
  # Bias
  VE_hat <- mean(result_table2[,"VE"], na.rm = TRUE)
  bias_VE <- mean(result_table2[,"VE"] -true_VE, na.rm = TRUE)
  mean_n_event <- mean(result_table[,"n_event"], na.rm = TRUE)
  mean_n_calendar_group <- mean(result_table[,"n_calendar_group"], na.rm = TRUE)
  power <- mean(result_table2$p_val < 0.05, na.rm = TRUE)
  
  

  
  performance <- data.frame(missing_estimate,convergence_issue, VE_hat, bias_VE, mean_n_event, mean_n_calendar_group, power)
  
  performance
}

# ------------------------------------------------------------------------------
# 7. Function to summarize results (for main analysis) ----------------------------
# ------------------------------------------------------------------------------

summary_sim2 <- function(true_VE = 0.6, result_table, n_sim=1000)
{
  true_IRR_V = 1 - true_VE
  true_est_V <- log(true_IRR_V)
  # Number of missing values of estimated beta1 (e.g due to convergence)
  missing_estimate <- sum(is.na(result_table$est_V)) + n_sim - nrow(result_table)
  convergence_issue <- sum(result_table$IRR_V > 50, na.rm = TRUE )
  
  result_table2 <- result_table[result_table$IRR_V < 50,]
  nsim2 = nrow(result_table2)
  # Mean number of events
  mean_n_event <- mean(result_table[,"n_event"], na.rm = TRUE)
  
  # Estimates
  est_V_hat <- mean(result_table2[,"est_V"], na.rm = TRUE)
  IRR_V_hat <- mean(result_table2[,"IRR_V"], na.rm = TRUE)
  VE_hat <- mean(result_table2[,"VE"], na.rm = TRUE)
  VE_mean_est=V <- 1 - exp(est_V_hat)
  
  
  # Variance 
  se_est_V_hat <- sqrt(1/(nsim2-1)*sum((result_table2[,"est_V"] - est_V_hat)^2, na.rm = TRUE)) #Empirical standard error
  se_est_V_hat_MCSE <- se_est_V_hat/sqrt(2*(nsim2-1)) # Monte Carlo standard error (MCSE) of empirical SE
  se_est_V_hat_low_CI <- se_est_V_hat - 1.96*se_est_V_hat_MCSE
  se_est_V_hat_up_CI <- se_est_V_hat + 1.96*se_est_V_hat_MCSE
  
  mod_se_est_V_hat <- sqrt(mean((result_table2[,"se_V"])^2, na.rm = TRUE)) # Model-based SE
  
  
  # Bias
  ## Absolute bias log scale
  bias_est_V <- mean(result_table2[,"est_V"]-true_est_V, na.rm = TRUE) 
  bias_est_V_MCSE <- sqrt(1/nsim2)*se_est_V_hat #  MCSE of absolute bias
  bias_est_V_low_CI <- bias_est_V - 1.96*bias_est_V_MCSE
  bias_est_V_up_CI <- bias_est_V + 1.96*bias_est_V_MCSE
  
  est_V_hat_low_CI <-  est_V_hat - 1.96*bias_est_V_MCSE
  est_V_hat_up_CI <-  est_V_hat + 1.96*bias_est_V_MCSE
  
  VE_mean_est_low_CI <- 1 - exp(est_V_hat_up_CI)
  VE_mean_est_up_CI <- 1 - exp(est_V_hat_low_CI)
  
  ## Relative bias log scale
  relative_bias_est_V <- mean((result_table2[,"est_V"]-true_est_V)/true_est_V, na.rm = TRUE)
  relative_bias_est_V_MCSE <- sqrt(1/(nsim2*(nsim2-1))*sum(((result_table2[,"est_V"]-true_est_V)/true_est_V-relative_bias_est_V)^2, na.rm = TRUE))
  relative_bias_est_V_low_CI <- relative_bias_est_V - 1.96*relative_bias_est_V_MCSE
  relative_bias_est_V_up_CI <- relative_bias_est_V + 1.96*relative_bias_est_V_MCSE
  
  ## Absolute bias IRR scale
  bias_IRR_V <- mean(result_table2[,"IRR_V"]-true_IRR_V, na.rm = TRUE) 
  bias_IRR_V_MCSE <- sqrt(1/((nsim2-1)*nsim2)*sum((result_table2[,"IRR_V"] - IRR_V_hat)^2, na.rm = TRUE))
  bias_IRR_V_low_CI <- bias_IRR_V - 1.96*bias_IRR_V_MCSE
  bias_IRR_V_up_CI <- bias_IRR_V + 1.96*bias_IRR_V_MCSE
  IRR_V_hat_low_CI <- IRR_V_hat - 1.96*bias_IRR_V_MCSE
  IRR_V_hat_up_CI <- IRR_V_hat + 1.96*bias_IRR_V_MCSE
  
  ## Relative bias IRR scale
  relative_bias_IRR_V <- mean((result_table2[,"IRR_V"]-true_IRR_V)/true_IRR_V, na.rm = TRUE) 
  relative_bias_IRR_V_MCSE <- sqrt(1/(nsim2*(nsim2-1))*sum(((result_table2[,"IRR_V"]-true_IRR_V)/true_IRR_V-relative_bias_IRR_V)^2, na.rm = TRUE))
  
  ## Absolute bias VE
  bias_VE <- mean(result_table2[,"VE"]-true_VE, na.rm = TRUE)
  bias_VE_MCSE <- sqrt(1/((nsim2-1)*nsim2)*sum((result_table2[,"VE"] - VE_hat)^2, na.rm = TRUE))
  bias_VE_low_CI <- bias_VE - 1.96*bias_VE_MCSE
  bias_VE_up_CI <- bias_VE + 1.96*bias_VE_MCSE
  VE_hat_low_CI <- VE_hat - 1.96*bias_VE_MCSE
  VE_hat_up_CI <- VE_hat + 1.96*bias_VE_MCSE

  # Coverage
  result_table2$coverage_IRR_V <- with(result_table2,
                                      IRR_V_CI_Lower <= true_IRR_V & IRR_V_CI_Upper >= true_IRR_V)
  coverage_irr_V <- mean(result_table2$coverage_IRR_V)
  coverage_irr_V_MCSE <- sqrt(coverage_irr_V*(1-coverage_irr_V)/nsim2)
  coverage_irr_V_low_CI <- coverage_irr_V - 1.96*coverage_irr_V_MCSE
  coverage_irr_V_up_CI <- coverage_irr_V + 1.96*coverage_irr_V_MCSE
  
  
  # Table of performance metrics
  performance <- data.frame(missing_estimate, convergence_issue, mean_n_event,
                            est_V_hat, est_V_hat_low_CI, est_V_hat_up_CI,
                            IRR_V_hat, IRR_V_hat_low_CI, IRR_V_hat_up_CI,
                            VE_hat, VE_hat_low_CI, VE_hat_up_CI,
                            VE_mean_est, VE_mean_est_low_CI, VE_mean_est_up_CI,
                            bias_est_V, bias_est_V_MCSE, bias_est_V_low_CI, bias_est_V_up_CI,
                            mod_se_est_V_hat, se_est_V_hat, se_est_V_hat_MCSE, se_est_V_hat_low_CI, se_est_V_hat_up_CI,
                            relative_bias_est_V, relative_bias_est_V_MCSE, relative_bias_est_V_low_CI, relative_bias_est_V_up_CI, 
                            bias_IRR_V, bias_IRR_V_MCSE, bias_IRR_V_low_CI, bias_IRR_V_up_CI,
                            relative_bias_IRR_V, relative_bias_IRR_V_MCSE,
                            bias_VE, bias_VE_MCSE, bias_VE_low_CI, bias_VE_up_CI,
                            coverage_irr_V, coverage_irr_V_MCSE, coverage_irr_V_low_CI, coverage_irr_V_up_CI
                            )  
  return(performance)
}

