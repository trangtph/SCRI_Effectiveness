if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman")
}


library(pacman)
pacman::p_load(
  here,
  lubridate,
  dplyr,
  tidyr,
  readr, 
  ggplot2,
  rio
)

options(scipen = 999)

# COVID-19 infection (old version with CDC data) -----------------------------------------------------------

covid_test_us <- read.csv(here("Data", "Percent_Positivity_of_COVID-19_US.csv"), sep=";")

covid_test_us <- covid_test_us %>% 
  mutate(report_week_end = as_date(ymd_hms(mmwrweek_end)),
         date_posted = as_date(mdy_hms(posted)), 
         number_tested = parse_number(number_tested))

covid_test_us_22 <- covid_test_us %>%
  filter(level == "National", 
         report_week_end > as.Date("2021-12-31") & report_week_end < as.Date("2023-01-01")) %>%
  mutate(nr_pos = round(percent_pos/100*number_tested)) %>%
  group_by(report_week_end) %>% slice_max(date_posted)

# Multiply the daily risk by 10 to account for under-reporting
covid_test_us_22 <- covid_test_us_22 %>%
  mutate(day = as.numeric(report_week_end) - as.numeric(as.Date("2022-01-01")) + 1,
         daily_risk_naive = nr_pos/332000000/7,
         daily_risk_adj = nr_pos*10/332000000/7) %>%
  select(report_week_end, number_tested, nr_pos, percent_pos, day, daily_risk_naive, daily_risk_adj) %>%
  ungroup()

covid_test_us_22 %>%
  ggplot(aes(x = report_week_end, y = nr_pos)) +
  #  geom_point() + 
  geom_line() +
  labs(x = "Week",
       y = "Number of SARS-CoV-2 positive tests") +
  theme_bw() +
  scale_x_date(date_breaks = "4 weeks", date_labels = "%d %b") + 
  scale_y_continuous(n.breaks = 10)

ggsave("covid_test_us_22.png", path = here("Plots", "Distribution"))

covid_test_us_22 %>%
  ggplot(aes(x = report_week_end, y = daily_risk_adj)) +
  #  geom_point() + 
  geom_line() +
  labs(x = "Calendar time (2022)",
       y = "Daily risk of SARS-CoV-2 infection (adjusted)") +
  theme_bw() +
  scale_x_date(date_breaks = "4 weeks", date_labels = "%d %b") + 
  scale_y_continuous(n.breaks = 10)
ggsave("covid_dist_us_22.png", 
       width = 18, height = 9, units = "cm",
       path = here("Plots", "Distribution"))

# Assign daily risk = risk of the last day of that week
covid_risk_us_22 <- covid_test_us_22 %>%
  select(day, daily_risk_adj) %>%
  arrange(day) %>%
  # add all days
  complete(day = 1:365) %>%
  # carry values backward
  fill(daily_risk_adj, .direction = "up")

export(covid_risk_us_22, here("Data", "covid_risk_us_22.csv"))

# interpolate linear change of daily risk (not run yet)
d <- import(here("Data", "covid_risk_us_22.csv"))
# anchor = first day of each run of equal values
is_anchor <- c(TRUE, diff(d$daily_risk_adj) != 0)
anc <- d[is_anchor, ]
# linear interpolation; flat extrapolation past the last anchor (rule = 2)
d$daily_risk_adj <- approx(x = anc$day, y = anc$daily_risk_adj,
                           xout = d$day, rule = 2)$y
export(d, here("Data", "covid_risk_us_22.csv"))


# COVID-19 infection (new version with ourworldindata data) ----
covid_incidence_us_22 <- import(here("Data", "daily-covid-19-incidence-US-2022.csv"))

covid_incidence_us_22 <- covid_incidence_us_22 %>%
  filter(Day >= as.Date("2022-01-01")) %>%
  rename(case_per_1M =  'New cases (per 1M)') 

covid_incidence_us_22 <- covid_incidence_us_22 %>% 
  mutate(
    total_case = case_per_1M/(10^6)*332*10^6, 
    day = seq(from = 1, to = nrow(covid_incidence_us_22)), 
    daily_risk = case_per_1M/(10^6))
export(covid_incidence_us_22, here("Data", "covid_incidence_us_22.csv"))

covid_incidence_us_22 %>%
  ggplot(aes(x = Day, y = daily_risk)) +
  #  geom_point() + 
  geom_line() +
  labs(x = "Calendar time (2022)",
       y = "Daily risk of COVID-19 infection") +
  theme_bw() +
  scale_x_date(date_breaks = "4 weeks", date_labels = "%d %b") + 
  scale_y_continuous(n.breaks = 10)

ggsave("covid_risk_us_22.png", 
       width = 18, height = 9, units = "cm",
       path = here("Plots", "Distribution"))


# Vaccine effectiveness --------------------------------------------------------

ve_us <- data.frame(
  day = c(0, seq(8, 162, by =7)),
  ve = c(0, 0.43, 0.67, 0.59, 0.48, 0.47, 0.47, 0.46, 0.45, 0.45, 
         0.44, 0.44, 0.43, 0.43, 0.42, 0.41, 0.41, 0.40, 0.40, 0.39, 0.38, 0.38, 0.37, 0.37)
)

ve_us %>%
  ggplot(aes(x = day, y = ve)) +
  geom_point() + 
  geom_line() +
  labs(x = "Days after vaccination",
       y = "VE on COVID-19 Hospitalization or death") +
  theme_bw() +
  scale_x_continuous(n.breaks = 13) + 
  scale_y_continuous(n.breaks = 10)

ggsave("ve_us.png", path = here("Plots", "Distribution"))

# Interpolate daily VE
# Create daily sequence
ve_daily <- data.frame(day = 0:162) %>%
  left_join(ve_us, by = "day")

# Linear interpolation between observed points
ve_daily <- ve_daily %>%
  mutate(
    ve = approx(
      x = ve_us$day,
      y = ve_us$ve,
      xout = day,
      rule = 2
    )$y
  )

# Enforce rule: VE = 0 from day 0 to 7 & = 67 from day 15 to 21
ve_daily <- ve_daily %>%
  mutate(
    ve = case_when(day <= 7 ~ 0,
                   day >= 15 & day <= 21 ~ ve_daily$ve[ve_daily$day == 15],
                   TRUE ~ ve)
  )

# 1. Fit linear waning model using last k points
k <- 15
ve_tail <- ve_us %>% tail(k)

fit <- lm(ve ~ day, data = ve_tail)

slope <- coef(fit)[["day"]]
intercept <- coef(fit)[["(Intercept)"]]

# 2. Day when VE reaches zero
day_zero <- -intercept / slope
day_zero <- ceiling(day_zero)

# 3. Extrapolate VE after day 162
ve_future <- data.frame(
  day = 163:day_zero
) %>%
  mutate(
    ve = intercept + slope * day,
    ve = pmax(ve, 0)
  )

# 4. Combine observed and extrapolated data
ve_full <- bind_rows(ve_daily, ve_future)

export(ve_full, here("Data", "ve_daily.csv"))

ve_full %>%
  ggplot(aes(x = day, y = ve)) +
  geom_point(data = ve_us, aes(x = day, y = ve)) + 
  geom_line() +
  labs(x = "Days after vaccination",
       y = "VE on COVID-19 Hospitalization or death") +
  theme_bw() +
  scale_x_continuous(n.breaks = 30, limits = c(0, 600)) + 
  scale_y_continuous(n.breaks = 10)

ggsave("ve_us2.png", path = here("Plots", "Distribution"),
       width = 20, height = 10, units = "cm")


# Vaccination administration --------------------------------------------------
covid_vacc_us_22 <- read.csv(here("Data", "daily-covid-19-vaccine-doses-administered-per-million-people.filtered", 
                                  "daily-covid-19-vaccine-doses-administered-per-million-people.csv"), 
                             header=TRUE)
covid_vacc_us_22 <- covid_vacc_us_22 %>%
  select(c(3,4)) %>%
  rename(date = Day, daily_dose_per1M = COVID.19.doses..daily..7.day.average..per.million.people.) %>%
  mutate(date = as.Date(date)) %>%
  filter(date >= as.Date("2022-01-01") & date <= as.Date("2022-12-31")) %>%
  mutate(day = as.numeric(date) - as.numeric(as.Date("2022-01-01")) + 1,
         prob = daily_dose_per1M/1000000)

export(covid_vacc_us_22, here("Data", "covid_vacc_us_22.csv"))

covid_vacc_us_22 %>%
  ggplot(aes(x = date, y = prob)) +
  #  geom_point() + 
  geom_line() +
  labs(x = "Calendar time (2022)",
       y = "Daily probability of COVID-19 vaccine administration") +
  theme_bw() +
  scale_x_date(date_breaks = "4 weeks", date_labels = "%d %b") + 
  scale_y_continuous(n.breaks = 10)

ggsave("covid_vacc_us_22.png", 
       width = 18, height = 9, units = "cm",
       path = here("Plots", "Distribution"))



  

  