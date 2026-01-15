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
  doFuture
)

source(here("Codes", "SCRI_main_functions.R"))
source(here("Codes", "SCRI_sim_workflow_functions.R"))
source(here("Codes", "SCRI_helper_functions.R"))

set.seed(123)
test_data <- cohort_data()

#Some diagnostic plots

## Vaccination date
hist(test_data$vacc_date,
     breaks = 40,
     col = "skyblue",
     main = "Distribution of vaccination days (truncated normal)",
     xlab = "Day of year",
     border = "white")

## mean daily infection probabilities
vaccinated_ids <- which(!is.na(test_data$vacc_date))
unvaccinated_ids <- which(is.na(test_data$vacc_date))

# Compute mean daily infection probabilities
mean_inf_vacc <- colMeans(test_data[vaccinated_ids, grep("^day_", names(test_data)), with = FALSE])
mean_inf_unvacc <- colMeans(test_data[unvaccinated_ids, grep("^day_", names(test_data)), with = FALSE])

# Plot
plot(1:365, mean_inf_unvacc, type = "l", lwd = 2, col = "red",
     xlab = "Day of year", ylab = "Mean daily infection probability",
     main = "Simulated infection risk over time",
     ylim = range(c(mean_inf_vacc, mean_inf_unvacc)))
lines(1:365, mean_inf_vacc, col = "blue", lwd = 2)

legend("topleft",
       legend = c("Unvaccinated", "Vaccinated"),
       col = c("red", "blue"),
       lty = 1, lwd = 2, bty = "n")

# Test the function to reformat and fit SCRI model
control_start =3 
control_end = 7 
risk_start = 17 
risk_end = 35 
infection_cols <- grep("^day_", names(test_data), value = TRUE) 
test_data[, infection_count := rowSums(.SD), .SDcols = infection_cols] 
vaccinated_infected <- test_data[!is.na(vacc_date) & infection_count > 0] 
vaccinated_infected_long <- data.table::melt( vaccinated_infected, id.vars = c("id", "vacc_date"), measure.vars = infection_cols, variable.name = "day", value.name = "infection" ) 
vaccinated_infected_long[, day_infection := as.integer(sub("day_", "", day))] 
vaccinated_infected_long2 <- vaccinated_infected_long[infection ==1, .(id, vacc_date, day_infection)] 
vaccinated_infected_long2[, `:=`(control_start_d = vacc_date + control_start, control_end_d = vacc_date + control_end, risk_start_d = vacc_date + risk_start, risk_end_d = vacc_date + risk_end)] 
vaccinated_infected_long2 <- vaccinated_infected_long2[ (day_infection >= control_start_d & day_infection <= control_end_d) | (day_infection >= risk_start_d & day_infection <= risk_end_d) ] 
vaccinated_infected_long3 <- formatdata(indiv = id, astart = control_start_d, # start of observation time = start of control period 
                                        aend = risk_end_d, # end of observation time = end of risk period 
                                        aevent = day_infection, adrug = vacc_date, aedrug = control_end_d, # end of the control period 
                                        expogrp =control_start, # start of the control period 
                                        washout = c(1,risk_start-1-control_end), #washout between the control and risk 
                                        data=vaccinated_infected_long2) 

vaccinated_infected_long3$vacc_date <- relevel(vaccinated_infected_long3$vacc_date, ref = "1") # 0 = risk period, 1 = control period, 2 = wash-out period 

# Model adjust for calendar time

long2 <- SCCS::formatdata(indiv = id,
                          astart = control_start_d,
                          aend = risk_end_d,
                          aevent = day_infection,
                          adrug = vacc_date,
                          aedrug = control_end_d,
                          expogrp =control_start,
                          washout = c(1,risk_start-1-control_end),
                          agegrp = seq(1, 365, by = 14),
                          data=vaccinated_infected_long2)
# Relevel exposure so control period = reference
long2$vacc_date <- relevel(factor(long2$vacc_date), ref = "1") # Because it is coded as: 0 = risk period, 1 = control period, 2 = wash-out period 
mod <- summary(clogit(event ~ vacc_date + age + strata(indivL) + offset(log(interval)), data = long2))
mod