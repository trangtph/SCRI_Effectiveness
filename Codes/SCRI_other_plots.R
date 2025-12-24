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

# ------------------------------------------------------------------------------
# 1. Daily infection risk: Gamma distribution ----------------------------------
# ------------------------------------------------------------------------------
if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
library(pacman)
pacman::p_load(here,
               rio
)

# Generate daily infection risk curve from a gamma distribution
gamma_risk_curve <- function(shape,
                             mode_day,
                             min_risk,
                             peak_risk,
                             n_days = 365,
                             dist_n) {
  
  if (shape <= 1) stop("Shape must be > 1 for a valid mode.")
  
  # scale parameter from mode
  scale <- (mode_day - 1) / (shape - 1)
  
  # domain: 0 to n_days-1 for gamma
  x_raw <- seq(0, n_days - 1, length.out = n_days)
  
  # gamma pdf on [0, ∞)
  pdf_vals <- dgamma(x_raw, shape = shape, scale = scale)
  
  # scale pdf into [min_risk, peak_risk]
  daily_risk <- min_risk + (peak_risk - min_risk) * (pdf_vals / max(pdf_vals))
  
  # cumulative incidence
  cum_incidence <- 1 - cumprod(1 - daily_risk)
  
  data.frame(
    day = 1:n_days,
    x_raw = x_raw,
    pdf = pdf_vals,
    daily_risk = daily_risk,
    cum_incidence = cum_incidence,
    shape = shape,
    mode_day = mode_day,
    scale = scale,
    dist_n = dist_n,
    min_risk = min_risk,
    peak_risk = peak_risk
  )
}

# Function to plot daily infection risk
plot_gamma_risk <- function(..., n_days = 365,
                            xlab = "Day",
                            ylab = "Daily risk",
                            main = "Daily infection risk (Gamma)") {
  
  scenarios <- list(...)
  colors <- c("black", "red", "blue", "darkgreen", "purple", "orange")
  
  # Plot first scenario
  df <- scenarios[[1]]
  plot(df$day, df$daily_risk, type = "l", lwd = 2, col = colors[1],
       ylim = c(0,0.005),
       xlab = xlab, ylab = ylab,
       main = main)
  
  # Add others
  if (length(scenarios) > 1) {
    for (i in 2:length(scenarios)) {
      df_i <- scenarios[[i]]
      lines(df_i$day, df_i$daily_risk, lwd = 2, col = colors[i])
    }
  }
  
  # Legend
  legend_labels <- sapply(scenarios, function(df)
    sprintf("Dist %d,shape=%.2f, mode=%d, peak risk=%.3f", df$dist_n[1], df$shape[1], df$mode_day[1], df$peak_risk[1])
  )
  
  legend("topright", legend = legend_labels, cex = 0.7,
         col = colors[1:length(scenarios)], lwd = 2)
}

# Function to plot cumulative incidence of infection
plot_gamma_ci <- function(...,
                          xlab = "Day",
                          ylab = "Cumulative incidence",
                          main = "Cumulative incidence (Gamma)") {
  
  scenarios <- list(...)
  colors <- c("black", "red", "blue", "darkgreen", "purple", "orange")
  
  df <- scenarios[[1]]
  plot(df$day, df$cum_incidence, type = "l", lwd = 2, col = colors[1],
       ylim = c(0, 0.8),
       xlab = xlab, ylab = ylab,
       main = main)
  
  if (length(scenarios) > 1) {
    for (i in 2:length(scenarios)) {
      df_i <- scenarios[[i]]
      lines(df_i$day, df_i$cum_incidence, lwd = 2, col = colors[i])
    }
  }
  
  legend_labels <- sapply(scenarios, function(df)
    sprintf("Dist %d,shape=%.2f, mode=%d, peak risk=%.3f", df$dist_n[1], df$shape[1], df$mode_day[1], df$peak_risk[1])
  )
  
  legend("topleft", legend = legend_labels, cex = 0.7,
         col = colors[1:length(scenarios)], lwd = 2)
}

# Specify scenarios
## Scenario 1
g1 <- gamma_risk_curve(
  shape = 2.5,
  mode_day = 100,
  min_risk = 2e-4,
  peak_risk = 2e-3,
  n_days = 365,
  dist_n = 1
)

## Scenario 2
g2 <- gamma_risk_curve(
  shape = 20,
  mode_day = 200,
  min_risk = 2e-4,
  peak_risk = 3e-3,
  n_days = 365,
  dist_n = 2
)
## Scenario 3
g3 <- gamma_risk_curve(
  shape = 10,
  mode_day = 300,
  min_risk = 2e-4,
  peak_risk = 3e-3,
  dist_n = 3
)

# Make the plots
png(here("Plots", "Distribution", "Infection_CI_3dist1.png"), width = 1000, height = 800, units = "px", res = 150)
plot_gamma_ci(g1, g2, g3)
dev.off()
png(here("Plots", "Distribution", "Infection_risk_3dist1.png"), width = 1000, height = 800, units = "px", res = 150)
plot_gamma_risk(g1, g2, g3)
dev.off()
png(here("Plots", "Distribution", "Infection_risk_basecase.png"), width = 1000, height = 800, units = "px", res = 150)
plot_gamma_risk(g1)
dev.off()
# ------------------------------------------------------------------------------
# 2. Vaccination date: Beta distribution ---------------------------------------
# ------------------------------------------------------------------------------

# --- Function to convert mean/sd on [L,U] to Beta(a,b) parameters ---
beta_params <- function(mean, sd, L = 1, U = 365) {
  m <- (mean - L) / (U - L)               # scaled mean in [0,1]
  v <- (sd^2) / ((U - L)^2)               # scaled variance
  
  # Solve for a and b
  tmp <- m * (1 - m) / v - 1
  a <- m * tmp
  b <- (1 - m) * tmp
  
  return(list(a = a, b = b))
}

# --- Function to generate random sample from truncated beta distribution
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


# --- Define scenarios ---
scenarios <- list(
  list(mean = 180, sd = 60),
  list(mean = 180, sd = 20),
  list(mean = 80, sd = 60),
  list(mean = 80, sd = 20)
)


# --- Plot the four scenarios ---
# Get parameters for all scenarios
params <- lapply(scenarios, function(s) beta_params(s$mean, s$sd))

# --- Create x-axis on [1,365] ---
x <- seq(1, 365, length.out = 1000)

# --- Compute densities ---
dens <- lapply(params, function(p) {
  dbeta((x - 1) / (365 - 1), p$a, p$b) / (365 - 1)
})


png(here("Plots", "Vaccination_4dist.png"), width = 1000, height = 800, units = "px", res = 150)
plot(x, dens[[1]], type = "l", lwd = 2, ylim=c(0, 0.02),
     col = "blue",
     ylab = "Density", xlab = "Day of vaccination",
     main = "Beta Distributions of Vaccination Date")

lines(x, dens[[2]], lwd = 2, col = "red")
lines(x, dens[[3]], lwd = 2, col = "darkgreen")
lines(x, dens[[4]], lwd = 2, col = "purple")

legend("topright",
       legend = c(
         "Dist 1: Mean=180, SD=60",
         "Dist 2: Mean=180, SD=20",
         "Dist 3: Mean=80, SD=60",
         "Dist 4: Mean=80, SD=20"
       ),
       col = c("blue", "red", "darkgreen", "purple"),
       lwd = 2,
       cex = 0.8)
dev.off()

# ------------------------------------------------------------------------------
# 3. Results of the simulation for sample size determination -------------------
# ------------------------------------------------------------------------------

power_results <- import(here("Results", "Summary", "Power_results.xlsx"))
power_results <- power_results %>% mutate(n_event_size = paste0(ceiling(mean_n_event), "\n(", size,")"),
                                          VE_hat = round(VE_hat, digits = 3))
#Power plot
png(here("Plots", "Power.png"), width = 15, height = 10, units = "cm", res = 150)

power_results %>% ggplot( mapping = aes(x = size, y = power, colour = methods)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  geom_hline(yintercept=0.8, linetype="dashed") + 
  labs(x = "Mean number of events (corresponding cohort size)", y = "Power") +
  scale_x_continuous(
    breaks = power_results$size,          
    labels = power_results$n_event_size 
  ) + 
  scale_color_discrete(name = "Model", labels = c("Adjust for calendar month", "No seasonality adjustment")) + 
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
  scale_color_discrete(name = "Model", labels = c("Adjust for calendar month", "No seasonality adjustment")) + 
  theme_bw() + 
  theme(legend.position = "bottom")
dev.off()