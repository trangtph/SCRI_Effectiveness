# --- Parameters ---
n_days <- 365
min_risk <- 5e-4      # lowest daily risk
peak_risk <- 5e-3     # highest daily risk

mode_day <- 180       # day where the peak occurs (1–365)
concentration <- 10   # controls sharpness (higher -> sharper peak)

# --- Compute Beta parameters (alpha, beta) from mode and concentration ---
mode <- (mode_day - 1) / (n_days - 1)  # convert day to [0,1]
s <- concentration                     # shorthand
if (s <= 2) stop("concentration must be > 2 to have an interior mode")

alpha <- mode * (s - 2) + 1
beta  <- (1 - mode) * (s - 2) + 1

# --- Generate daily risks ---
x <- seq(0, 1, length.out = n_days)
pdf_vals <- dbeta(x, alpha, beta)

# Scale so min = min_risk and max = peak_risk
scaled_risk <- min_risk + (peak_risk - min_risk) * (pdf_vals / max(pdf_vals))

# --- Create data frame ---
df <- data.frame(
  day = 1:n_days,
  x = x,
  alpha = alpha,
  beta = beta,
  pdf = pdf_vals,
  daily_risk = scaled_risk
)

# --- Plot ---
plot(
  df$day, df$daily_risk, type = "l", lwd = 2,
  xlab = "Day", ylab = "Daily risk (probability)",
  main = sprintf("365-day COVID-19 infection risk shaped by Beta(%.2f, %.2f)\nPeak ≈ day %d, concentration = %d",
                 alpha, beta, mode_day, concentration)
)
grid()

# --- Display summary ---
cat("Alpha =", round(alpha, 2), "\n")
cat("Beta  =", round(beta, 2), "\n")
cat("Observed min risk =", min(df$daily_risk), "\n")
cat("Observed max risk =", max(df$daily_risk), "\n")

# ---- Skewed beta -----------------------------------------------------------
# ---------------------------------------------------------------------------

# Skewed Beta example: long right tail by setting alpha < beta
n_days <- 365
min_risk <- 5e-4
peak_risk <- 5e-3

# Peak day and "concentration" (sharpness)
mode_day <- 140        # peak earlier in year
concentration <- 10    # larger => narrower overall mass

# compute alpha/beta from mode and concentration (same formula as before)
mode <- (mode_day - 1) / (n_days - 1)
s <- concentration
if (s <= 2) stop("concentration must be > 2")
alpha <- mode * (s - 2) + 1
beta  <- (1 - mode) * (s - 2) + 1

# adjust to make alpha smaller (increase skew) if desired
# e.g. multiply beta by a factor to lengthen right tail:
beta <- beta * 1.6   # increase this factor to make tail longer

x <- seq(0, 1, length.out = n_days)
pdf_vals <- dbeta(x, alpha, beta)
scaled_risk <- min_risk + (peak_risk - min_risk) * (pdf_vals / max(pdf_vals))

plot(1:n_days, scaled_risk, type = "l", lwd = 2,
     xlab = "Day", ylab = "Daily risk", main = sprintf("Skewed Beta: alpha=%.2f beta=%.2f", alpha, beta))
grid()
cat("alpha=", round(alpha,2), " beta=", round(beta,2), "\n")
cat("min, max:", min(scaled_risk), max(scaled_risk), "\n")


# Gamma ----------------------------------------------------------------

# Gamma example: use a shifted gamma PDF to get a long slow waning tail
n_days <- 365
min_risk <- 2e-4
peak_risk <- 2e-3

# Choose gamma shape (k) and target peak day (mode_day)
shape <- 2.5        # >1 to have a mode (1.5-5 typical). smaller -> longer right tail
mode_day <- 140     # desired peak day (must be >= 1)

# For gamma, mode = (shape - 1) * scale  (if shape > 1)
scale <- (mode_day - 1) / (shape - 1)   # solve for scale so mode ≈ mode_day-1 (x starting at 0)

# Build x from 0..(n_days-1) so that "mode_day" maps to x = mode_day-1
x_raw <- seq(0, n_days - 1, length.out = n_days)
pdf_vals <- dgamma(x_raw, shape = shape, scale = scale)

# Optionally shift slightly if you want the mode exactly at an integer day:
# (not necessary if using continuous days)

# Scale to requested min and peak
scaled_risk <- min_risk + (peak_risk - min_risk) * (pdf_vals / max(pdf_vals))

df <- data.frame(
  day = 1:n_days,
  quantile_raw = x_raw,
  pdf = pdf_vals,
  daily_risk = scaled_risk,
  cum_incidence = 1 - cumprod(1 - scaled_risk)
)


plot(df$day, df$cum_incidence, type = "l", lwd = 2, col = "red",
     xlab = "Day", ylab = "Cumulative incidence", main = "Cumulative incidence of infection")


plot(1:n_days, scaled_risk, type = "l", lwd = 2,
     xlab = "Day", ylab = "Daily risk",
     main = sprintf("Gamma distribution: shape=%.2f scale=%.2f \n (mode≈%d days)", shape, scale, mode_day))
grid()
cat("gamma shape=", shape, " scale=", round(scale,2), "\n")
cat("min, max:", min(scaled_risk), max(scaled_risk), "\n")

