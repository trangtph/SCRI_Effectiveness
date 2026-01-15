library(pacman)
pacman::p_load(
  here,
  ggplot2,
  rio,
  tidyverse)

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
  # rescale to [low, U]
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


set.seed(123)

s1 <- sample_beta_scaled(10000, mean = 180, sd = 60)
s2 <- sample_beta_scaled(10000, mean = 180, sd = 20)
s3 <- sample_beta_scaled(10000, mean = 80, sd = 60)
s4 <- sample_beta_scaled(10000, mean = 80, sd = 20)
par(mfrow=c(2,2))

hist(s1, main="Mean=180, SD=60", xlab="Day", breaks = 20)
hist(s2, main="Mean=180, SD=20",  xlab="Day", breaks = 20)
hist(s3, main="Mean=80, SD=60", xlab="Day", breaks = 20)
hist(s4, main="Mean=80, SD=20",  xlab="Day", breaks = 20)
