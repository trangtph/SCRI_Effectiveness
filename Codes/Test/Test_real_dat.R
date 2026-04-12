set.seed(1001)
test_cohort <- cohort_data_real(n_indiv = 20000)
test_cohort_sccs <- run_SCRI(dat = test_cohort,
                             method = "no_calendar",
                             n_days = 365,
                             rep = 1,
                             control_start = 3,
                             control_end   = 7,
                             risk_start    = 15,
                             risk_end      = 21,
                             start_calendar = 31, #Start of the second interval, the 1st interval start at the start of observation time
                             calendar_interval = 30)

hist(test_cohort$vacc_date)


dt <- as.data.table(test_cohort)

# reshape infection columns
long_inf <- melt(
  dt,
  id.vars = c("id", "vacc_date"),
  measure.vars = patterns("^day_"),
  variable.name = "day",
  value.name = "infection"
)

# convert "day_23" -> 23
long_inf[, infect_date := as.integer(sub("day_", "", day))]

# keep infection events only
sccs_dat <- long_inf[infection == 1,
                     .(vacc_date = first(vacc_date),
                       infect_date),
                     by = .(id, infect_date)]

setorder(sccs_dat, id, infect_date)

sccs_dat2 <- sccs_dat %>%
  select(1:3) %>%
  mutate(inf_since_vac = infect_date - vacc_date) %>% filter(inf_since_vac >=3 & inf_since_vac <=7)

# number of infections per individual
inf_count <- sccs_dat[, .N, by = id]

# distribution of number of infections
infection_distribution <- inf_count[, .N, by = N][order(N)]

infection_distribution

hist(sccs_dat$infect_date)

# Test the simulation workflow -----
perform_one_run_real(seed=1001, rep = 1, scen = scen_real[1,], 
                     methods = c("no_calendar", "calendar_30d", "calendar_7d"), 
                     output_dir = here("Results","Raw_results_real"))


## Run the simulation
set.seed(20251218)
plan(multisession, workers = 4)
n_sim <- 2

sim_one_scenario_real(scen = scen_real[1,], scenario_table = scen_real, n_sim= n_sim, 
                      seeds = get_seeds(n_sim = n_sim, scenario_table = scen_real), 
                      methods = c("no_calendar", "calendar_30d", "calendar_7d"), 
                      output_dir = here("Results","Raw_results_real"))
methods <- c("no_calendar", "calendar_30d", "calendar_7d")
full_simulation_real(scenario_table = scen_real, 
                n_sim = n_sim, 
                seeds = get_seeds(n_sim, scenario_table = scen_real), 
                methods = methods, 
                output_dir = here("Results","Raw_results_real"))

