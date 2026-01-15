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
  
  # cumulative incidence (optional but useful)
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