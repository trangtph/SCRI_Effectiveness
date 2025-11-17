# --- Parameters ---
n_days <- 365-78
mean_day <- 180
sd_day <- 100

# Define day sequence
day <- 1:n_days

# Compute normal density truncated at [1, 365]
dens <- dnorm(day, mean = mean_day, sd = sd_day)

# Scale so it fits nicely in [0, 1] range or desired probability scale
dens <- dens / max(dens)  # normalized so peak = 1

# Plot
plot(
  day, dens, type = "l", lwd = 2,
  xlab = "Day", ylab = "Scaled density",
  main = sprintf("Vaccination date, Normal distribution \n(mean = %d, sd = %d)", mean_day, sd_day)
)
grid()

# Optional: check tails
cat("Density near day 1:", dens[1], "\n")
cat("Density near day 365:", dens[365], "\n")
