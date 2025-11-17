##########################################
### Project: SCRI design for           ###
### vaccine effectiveness              ###
##########################################

###################################
### Functions for               ###
### producing plots             ###
###################################

if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman")
}


library(pacman)
pacman::p_load(
  here,
  ggplot2,
  rio,
  tidyverse)

options(scipen = 999)

# Daily infection risk: Gamma distribution -------------------------------------
n_days <- 365
min_risk <- 2e-4
peak_risk <- 2e-3

# Choose gamma shape (k) and target peak day (mode_day)
shape <- 2.5        # >1 to have a mode. smaller -> longer right tail
mode_day <- 140     # desired peak day (must be >= 1)

# For gamma, mode = (shape - 1) * scale  (if shape > 1)
scale <- (mode_day - 1) / (shape - 1)

# Build x from 0..(n_days-1)
x_raw <- seq(0, n_days - 1, length.out = n_days)
pdf_vals <- dgamma(x_raw, shape = shape, scale = scale)

# Scale to desired min and peak
scaled_risk <- min_risk + (peak_risk - min_risk) * (pdf_vals / max(pdf_vals))

df <- data.frame(
  day = 1:n_days,
  quantile_raw = x_raw,
  pdf = pdf_vals,
  daily_risk = scaled_risk,
  cum_incidence = 1 - cumprod(1 - scaled_risk)
)

# Plot cumulative incidence and distribution of daily infection risk

png(here("Plots", "CI_infection_gamma.png"), width = 800, height = 800, units = "px", res = 150)
plot(df$day, df$cum_incidence, type = "l", lwd = 2, col = "red",
     xlab = "Day", ylab = "Cumulative incidence", main = "Cumulative incidence of infection \n Gamma distribution")
dev.off()

png(here("Plots", "Infection_risk_gamma.png"), width = 1024, height = 1094, units = "px", res = 150)
plot(1:n_days, scaled_risk, type = "l", lwd = 2,
     xlab = "Day", ylab = "Daily risk",
     main = sprintf("Gamma distribution: shape=%.2f scale=%.2f \n (mode≈%d days)", shape, scale, mode_day))
dev.off()

# Vaccination date: Normal distribution ----------------------------------------

n_days <- 365-78
mean_day <- 180
sd_day <- 100
day <- 1:n_days

# Compute normal density truncated at [1, 287]
dens <- dnorm(day, mean = mean_day, sd = sd_day)

# Scale so it fits in [0, 1] range or desired probability scale
dens <- dens / max(dens)  # normalized so peak = 1

# Plot
png(here("Plots", "Vaccination_normal"), width = 1024, height = 1094, units = "px", res = 150)
plot(
  day, dens, type = "l", lwd = 2,
  xlab = "Day", ylab = "Scaled density",
  main = sprintf("Vaccination date, Normal distribution \n(mean = %d, sd = %d)", mean_day, sd_day)
)
dev.off()

# Results of the simulation for sample size determination ---------------------

power_results <- import(here("Results", "Summary", "Power_results.xlsx"))
power_results <- power_results %>% mutate(n_event_size = paste0(ceiling(mean_n_event), "\n(", size,")"),
                                          VE_hat = round(VE_hat, digits = 3))
#Power plot
png(here("Plots", "Power.png"), width = 15, height = 10, units = "cm", res = 150)

power_results %>% ggplot( mapping = aes(x = size, y = VE_hat, colour = methods)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  geom_hline(yintercept=0.8, linetype="dashed") + 
  labs(x = "Mean number of events (corresponding cohort size)", y = "Power") +
  scale_x_continuous(
    breaks = power_results$size,          
    labels = power_results$n_event_size 
  ) + 
  scale_color_discrete(name = "Model", labels = c("Adjust for calendar", "No seasonality adjustment")) + 
  theme_bw() + 
  theme(legend.position = "bottom")
dev.off()

# Bias plot
png(here("Plots", "Bias.png"), width = 15, height = 10, units = "cm", res = 150)

power_results %>% ggplot( mapping = aes(x = size, y = VE_hat, colour = methods)) +
  geom_point(size = 1.5) +
  geom_line(size = 1) +
  geom_hline(yintercept=0.6, linetype="dashed") + 
  labs(x = "Mean number of events (corresponding cohort size)", y = "Estimated Vaccine Effectiveness \n True value = 0.6") +
  scale_x_continuous(
    breaks = power_results$size,          # the actual values used as breaks
    labels = power_results$n_event_size         # labels pulled from another variable
  ) + 
  scale_color_discrete(name = "Model", labels = c("Adjust for calendar", "No seasonality adjustment")) + 
  theme_bw() + 
  theme(legend.position = "bottom")
dev.off()