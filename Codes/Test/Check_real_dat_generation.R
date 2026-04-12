set.seed(493615992, kind = "L'Ecuyer-CMRG")
scen <- scen_table[20,]
data_test <- cohort_data(n_indiv = scen[['sample_size']], 
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

infection_cols <- grep("^day_", names(data_test), value = TRUE)
#Keep vaccinated individuals with ≥1 infection
data_test[, infection_count := rowSums(.SD), .SDcols = infection_cols]
vaccinated_infected <- data_test[!is.na(vacc_date) & infection_count > 0]

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
  control_start_d = vacc_date + 3,
  control_end_d   = vacc_date + 7,
  risk_start_d    = vacc_date + 17,
  risk_end_d      = vacc_date + 35
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
fixed_agegrp <- seq(31, 365 - 29, by = 30)
min_start <- min(long$control_start_d, na.rm = TRUE)
max_end   <- max(long$risk_end_d, na.rm = TRUE)
agegrp <- fixed_agegrp[
  fixed_agegrp >= min_start +29 & fixed_agegrp <= max_end -29
]
agegrp <- NULL

long2 <- SCCS::formatdata(indiv = id,
                          astart = control_start_d,
                          aend = risk_end_d,
                          aevent = day_infection,
                          adrug = vacc_date,
                          aedrug = control_end_d,
                          expogrp =3,
                          washout = c(1,17-1-7),
                          agegrp = agegrp,
                          data=long)

long2$vacc_date <- relevel(factor(long2$vacc_date), ref = "1") # Because it is coded as: 0 = risk period, 1 = control period, 2 = wash-out period 

n_age <- nlevels(factor(long2$age))
base_formula <- event ~ vacc_date + strata(indivL) + offset(log(interval))

if (n_age > 1) {
  form <- update(base_formula, . ~ . + age)
} else {
  form <- base_formula
}

mod <- summary(clogit(form, data = long2))

mod <- summary(clogit(event ~ vacc_date + age + strata(indivL) + offset(log(interval)), data = long2))


out <- run_SCRI(dat = data_test, rep = 1, method = "calendar_adjustment",
                n_days = 365,
                control_start = scen[['control_start']],
                control_end   = scen[['control_end']],
                risk_start    = scen[['risk_start']],
                risk_end      = scen[['risk_end']],
                start_calendar = 31,
                calendar_interval = 30)