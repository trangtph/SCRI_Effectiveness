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
                            xlab = "Days since start of follow-up",
                            ylab = "Daily risk",
                            main = "Daily infection risk (Gamma distribution)") {
  
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
                          xlab = "Days since start of follow-up",
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


png(here("Plots", "Distribution", "Vaccination_4dist.png"), width = 1000, height = 800, units = "px", res = 150)
plot(x, dens[[1]], type = "l", lwd = 2, ylim=c(0, 0.02),
     col = "blue",
     ylab = "Density", xlab = "Days since the start of follow-up",
     main = "Vaccination Date (Beta Distributions)")

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

png(here("Plots", "Distribution", "Vaccination_4dist2.png"),
    width = 1000, height = 800, units = "px", res = 150)

plot(x, dens[[1]], type = "l", lwd = 2, ylim = c(0, 0.02),
     col = "black",
     ylab = "Density", xlab = "Days since the start of follow-up",
     main = "Vaccination Date (Uniform/Beta Distributions)")

lines(x, dens[[2]], lwd = 2, col = "red")
lines(x, dens[[3]], lwd = 2, col = "blue")
lines(x, dens[[4]], lwd = 2, col = "darkgreen")

## ---- ADD uniform distribution ---------------------------------------------
uniform_height <- 1 / 287

lines(c(1, 287),
      c(uniform_height, uniform_height),
      lwd = 2,
      col = "grey",
      lty = 2)

legend("topright",
       legend = c(
         "Dist 0: Uniform [1–287]",
         "Dist 1: Mean=180, SD=60",
         "Dist 2: Mean=180, SD=20",
         "Dist 3: Mean=80, SD=60",
         "Dist 4: Mean=80, SD=20"
       ),
       col = c("grey", "black", "red", "blue", "darkgreen"),
       lwd = 2,
       lty = c(2, 1, 1, 1, 1),
       cex = 0.8)

dev.off()

# ------------------------------------------------------------------------------
# 3. Results of the simulation for sample size determination -------------------
# ------------------------------------------------------------------------------

power_results <- import(here("Results", "Summary", "Power_results_20260106.xlsx"))
power_results <- power_results %>% mutate(n_event_size = paste0(ceiling(mean_n_event), "\n(", size,")"),
                                          VE_hat = round(VE_hat, digits = 3))
#Power plot
png(here("Plots", "Sample size calculation", "Power_20260106.png"), width = 20, height = 10, units = "cm", res = 150)

power_results %>% ggplot( mapping = aes(x = size, y = power)) +
  geom_point(size = 3, color = "#2f357c") +
  geom_line(linewidth = 1, color = "#2f357c") +
  geom_hline(yintercept=0.8, linetype="dashed") + 
  labs(x = "Mean number of events (corresponding cohort size)", y = "Power") +
  scale_x_continuous(
    breaks = power_results$size,          
    labels = power_results$n_event_size 
  ) +
  theme_bw() 
dev.off()

# Bias plot
png(here("Plots", "Sample size calculation", "Bias_20260106.png"), width = 20, height = 10, units = "cm", res = 150)

power_results %>% ggplot( mapping = aes(x = size, y = VE_hat)) +
  geom_point(size = 1.5, color = "#2f357c") +
  geom_line(linewidth = 1, color = "#2f357c") +
  geom_hline(yintercept=0.6, linetype="dashed") + 
  labs(x = "Mean number of events (corresponding cohort size)", y = "Estimated Vaccine Effectiveness \n True value = 0.6") +
  scale_x_continuous(
    breaks = power_results$size,        
    labels = power_results$n_event_size         
  ) + 
  scale_y_continuous(breaks = seq(0.48, 0.6, by = 0.01), limits = c(0.48, 0.6)) +
  theme_bw() 
dev.off()

# ------------------------------------------------------------------------------
# 4. Combination grid: 3 infection-risk curves x 5 vaccination-date curves -----
# ------------------------------------------------------------------------------


Tmax <- 365
days <- 1:Tmax

## Infection risk: scaled Gamma daily-risk curves -----
# Gamma mode = (shape - 1) * scale  ->  scale = mode / (shape - 1)
# Curve is the Gamma density rescaled so its peak equals the target peak risk.
inf_params <- list(
  A = list(shape = 2.5,  mode = 100, peak = 0.002, col = "black"),
  B = list(shape = 20.0, mode = 200, peak = 0.003, col = "#bc371b"),
  C = list(shape = 10.0, mode = 300, peak = 0.003, col = "#6c5d9e")
)

infection_risk <- function(t, shape, mode, peak) {
  scale    <- mode / (shape - 1)
  dens     <- dgamma(t,    shape = shape, scale = scale)
  peakdens <- dgamma(mode, shape = shape, scale = scale)
  peak * dens / peakdens
}

## ----- Vaccination date: Uniform + Beta (method of moments on [0, Tmax]) -----
# For Y = Tmax * X, X ~ Beta(a,b):
#   mu = mean/Tmax ,  v = (sd/Tmax)^2
#   a + b = mu(1-mu)/v - 1 ,  a = mu(a+b) ,  b = (1-mu)(a+b)
beta_from_moments <- function(mean, sd, T) {
  mu <- mean / T
  v  <- (sd / T)^2
  ab <- mu * (1 - mu) / v - 1
  c(a = mu * ab, b = (1 - mu) * ab)
}

vax_params <- list(
  V0 = list(type = "unif", lo = 1, hi = 365,           col = "#4d4d4d", lty = 1, lab = "V0: Uniform [1-365]"),
  V1 = list(type = "beta", mean = 180, sd = 60,        col = "#ea9e0a", lty = 1, lab = "V1: Mean 180, SD 60"),
  V2 = list(type = "beta", mean = 180, sd = 20,        col = "#d47261", lty = 1, lab = "V2: Mean 180, SD 20"),
  V3 = list(type = "beta", mean = 80,  sd = 60,        col = "#2f5328", lty = 1, lab = "V3: Mean 80, SD 60"),
  V4 = list(type = "beta", mean = 80,  sd = 20,        col = "#17154f", lty = 1, lab = "V4: Mean 80, SD 20")
)

vax_density <- function(t, p, T) {
  if (p$type == "unif") {
    ifelse(t >= p$lo & t <= p$hi, 1 / (p$hi - p$lo), 0)
  } else {
    ab <- beta_from_moments(p$mean, p$sd, T)
    dbeta(t / T, ab["a"], ab["b"]) / T
  }
}

## ----- Fixed axis limits so every panel is directly comparable -----
risk_max <- 0.005    # left axis (daily infection risk)
dens_max <- 0.021    # right axis (vaccination density)

inf_names <- names(inf_params)
vax_names <- names(vax_params)

png(here("Plots", "Distribution", "inf_vac_combination.png"), width = 2100, height = 1350, res = 300)

par(mfrow = c(3, 5), # fill by row: 3 rows x 5 columns
    mar = c(1, 1, 1, 1), # inner margin between cells
    oma = c(3, 3.5, 1, 3.5)) # outer margin

panel_no <- 0
for (i in seq_along(inf_names)) {
  for (j in seq_along(vax_names)) {
    panel_no <- panel_no + 1
    ip <- inf_params[[inf_names[i]]]
    vp <- vax_params[[vax_names[j]]]
    
    risk  <- infection_risk(days, ip$shape, ip$mode, ip$peak)
    vdens <- vax_density(days, vp, Tmax)
    
    # --- Infection risk (left axis) ---
    plot(days, risk, type = "n",
         ylim = c(0, risk_max), xlim = c(0, Tmax),
         axes = FALSE, xlab = "", ylab = "")
    
    # shaded vaccination density (rescaled to left axis for visual fill)
    vfill <- vdens / dens_max * risk_max
    polygon(c(days[1], days, days[length(days)]),
            c(0, vfill, 0),
            col = adjustcolor(vp$col, alpha.f = 0.18), border = NA)
    
    # infection risk line on top
    lines(days, risk, col = ip$col, lwd = 2.4)
    # vaccination density outline
    lines(days, vfill, col = vp$col, lwd = 1.8, lty = vp$lty)
    
    box(col = "grey70")
    title(main = paste0("Scenario ", 4+(i-1)*5+(j-1),": Inf ", inf_names[i], " \u00d7 ", vax_names[j]),
          cex.main = 0.7, font.main = 2)
    
    # axes only on edges to keep the grid clean
    if (j == 1) axis(2, at = seq(0, risk_max, 0.001), las = 1, cex.axis = 0.7, col = "grey50")
    if (j == length(vax_names)) {
      axis(4, at = seq(0, risk_max, length.out = 5),
           labels = signif(seq(0, dens_max, length.out = 5), 2),
           las = 1, cex.axis = 0.7, col = "grey50")
    }
    if (i == length(inf_names)) axis(1, at = seq(0, 300, 100), cex.axis = 0.7, col = "grey50")
  }
}

# outer labels
mtext("Days since start of follow-up", side = 1, outer = TRUE, line = 1.6, cex = 0.7)
mtext("Daily infection risk (left axis)", side = 2, outer = TRUE, line = 2.6, cex = 0.7)
mtext("Vaccination density (right axis)", side = 4, outer = TRUE, line = 2.6, cex = 0.7)
#mtext("15 scenarios: infection-risk curve (line) overlaid with vaccination-date curve (shaded)",
#      side = 3, outer = TRUE, line = 1.8, cex = 1, font = 2)

dev.off()

## ----- Quick numeric check that the Beta fits reproduce target moments -----
chk <- do.call(rbind, lapply(vax_names[-1], function(nm) {
  p  <- vax_params[[nm]]
  ab <- beta_from_moments(p$mean, p$sd, Tmax)
  m  <- Tmax * ab["a"] / (ab["a"] + ab["b"])
  v  <- Tmax^2 * ab["a"] * ab["b"] /
    ((ab["a"] + ab["b"])^2 * (ab["a"] + ab["b"] + 1))
  data.frame(dist = nm, target_mean = p$mean, fit_mean = round(m, 1),
             target_sd = p$sd, fit_sd = round(sqrt(v), 1),
             a = round(ab["a"], 2), b = round(ab["b"], 2))
}))
print(chk, row.names = FALSE)