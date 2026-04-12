##########################################
### Project: SCRI design for           ###
### vaccine effectiveness              ###
##########################################

###################################
### Functions for               ###
### simulating and analyzing    ###
### data.                       ###
###################################

# ------------------------------------------------------------------------------
# Function for daily baseline infection risk -----------------------------------
# ------------------------------------------------------------------------------
# This function create daily baseline infection risk for n_days = 365 days. 
# The daily baseline infection risk follows a gamma distribution, peaking at day 140.

baseline_risk <- function(n_days = 365, # number of follow-up days for cohort
                          gamma_shape = 2.5,
                          gamma_mode = 140, # Day of peak infection risk
                          min_risk = 2e-4,
                          peak_risk = 2e-3){
  scale <- (gamma_mode - 1) / (gamma_shape - 1)
  quantile <- seq(0, n_days - 1, length.out = n_days)
  pdf_vals <- dgamma(quantile, shape = gamma_shape, scale = scale)
  baseline_risk <- min_risk + (peak_risk - min_risk) * (pdf_vals / max(pdf_vals))
  return(baseline_risk)
}

# ------------------------------------------------------------------------------
# Function to calculate daily vaccine effectiveness, day 0-150 post-vaccination 
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
# Function to generate cohort data ---------------------------------------------
# ------------------------------------------------------------------------------

# This function generate data for a cohort of n_indiv individuals observed for 
# n_days (=365 days). Data include vaccination status (assuming 1 dose of vaccine, 
# 80% vaccinated), and daily infection status.

cohort_data <- function(n_indiv = 10000, 
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
                        baseline_risk_vec = baseline_risk()) {
  # Generate exposure information -------------------
  # Number vaccinated
  n_vacc <- round(n_indiv * prop_vacc)
  
  # Sample vaccination dates for vaccinated individuals from truncated normal distribution
  vacc_date <- round(rtnorm(n = n_vacc,
                            a = 1,
                            b = n_days - max_risk_period, # so that SCRI risk period is always fully observed
                            mean = mean_vacc_date,
                            sd = sd_vacc_date))
  vacc_date = c(vacc_date, rep(NA, n_indiv - n_vacc)) # add NA vacc_date for unvaccinated
  
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
# Function to reshape data to SCRI-compatible format and fit SCRI model---------
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
                     start_calendar = 1,
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
      row.names = NULL
    ))
    
  }
  
  else if (method == "calendar_adjustment"){
    
    long2 <- SCCS::formatdata(indiv = id,
               astart = control_start_d,
               aend = risk_end_d,
               aevent = day_infection,
               adrug = vacc_date,
               aedrug = control_end_d,
               expogrp =control_start,
               washout = c(1,risk_start-1-control_end),
               agegrp = seq(start_calendar, n_days, by = calendar_interval),
               data=long)
    # Relevel exposure so control period = reference
    long2$vacc_date <- relevel(factor(long2$vacc_date), ref = "1") # Because it is coded as: 0 = risk period, 1 = control period, 2 = wash-out period 
    mod <- summary(clogit(event ~ vacc_date + age + strata(indivL) + offset(log(interval)), data = long2))
    
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
      row.names = NULL
    ))
  }
  
  else {
    stop("Error: Method should be either no_calendar, or calendar_adjustment")
  }

}

# Function to summarise results

summary_sim <- function(true_VE, result_table, n_sim=1000)
{

  # Number of missing values of estimated beta1 (e.g due to convergence)
  missing_estimate <- sum(is.na(result_table$est_V)) + n_sim - nrow(result_table)
  
  # Bias
  VE_hat <- mean(result_table[,"VE"], na.rm = TRUE)
  bias_VE <- mean(result_table[,"VE"] -true_VE, na.rm = TRUE)
  mean_n_event <- mean(result_table[,"n_event"], na.rm = TRUE)
  power <- mean(result_table$p_val < 0.05, na.rm = TRUE)
  

  
  performance <- data.frame(missing_estimate, VE_hat, bias_VE, mean_n_event, power)
  
  performance
}

