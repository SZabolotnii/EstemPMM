# Test script for SARMA/SARIMA PMM2 implementation
# This script tests the unified approach for seasonal models

library(EstemPMM)

cat("=================================================================\n")
cat("Testing SARMA/SARIMA PMM2 Implementation\n")
cat("=================================================================\n\n")

# Set seed for reproducibility
set.seed(123)

# Test 1: SARMA(1,0)×(1,0)_12 - Pure AR+SAR model
cat("Test 1: SARMA(1,0)×(1,0)_12 model (AR + SAR only)\n")
cat("-----------------------------------------------------------\n")

# Simulate data
n <- 200
phi <- 0.5   # AR(1) coefficient
Phi <- 0.6   # SAR(1) coefficient
s <- 12      # Seasonal period

# Generate with asymmetric gamma errors
innov <- rgamma(n, shape = 2, scale = 1) - 2
y1 <- numeric(n)
for (t in (s+2):n) {
  y1[t] <- phi * y1[t-1] + Phi * y1[t-s] + innov[t]
}
y1 <- y1[(s+2):n]  # Remove burn-in

# Fit SARMA model
fit1 <- sarma_pmm2(y1, order = c(1, 1, 0, 0), season = list(period = 12),
                   method = "pmm2", verbose = TRUE)

cat("\nModel summary:\n")
summary(fit1)

cat("\n\nTrue coefficients: AR1 =", phi, ", SAR1 =", Phi, "\n")
cat("PMM2 estimates: AR1 =", fit1@coefficients[1], ", SAR1 =", fit1@coefficients[2], "\n\n")

# Test 2: SARMA(0,0)×(1,1)_12 - Pure seasonal model
cat("\n=================================================================\n")
cat("Test 2: SARMA(0,0)×(1,1)_12 model (SAR + SMA only)\n")
cat("-----------------------------------------------------------\n")

# Simulate pure seasonal SARMA
Phi2 <- 0.7
Theta2 <- 0.4

innov2 <- rgamma(n, shape = 2, scale = 1) - 2
y2 <- numeric(n)
eps <- numeric(n)
for (t in 1:n) {
  sar_term <- if (t > s) Phi2 * y2[t-s] else 0
  sma_term <- if (t > s) Theta2 * eps[t-s] else 0
  eps[t] <- innov2[t]
  y2[t] <- sar_term + sma_term + eps[t]
}
y2 <- y2[(s+10):n]  # Remove burn-in

# Fit model
fit2 <- sarma_pmm2(y2, order = c(0, 1, 0, 1), season = list(period = 12),
                   method = "pmm2", verbose = FALSE)

cat("\nModel summary:\n")
summary(fit2)

cat("\n\nTrue coefficients: SAR1 =", Phi2, ", SMA1 =", Theta2, "\n")
cat("PMM2 estimates: SAR1 =", fit2@coefficients[1], ", SMA1 =", fit2@coefficients[2], "\n\n")

# Test 3: Full SARMA(1,1)×(1,1)_12 model
cat("\n=================================================================\n")
cat("Test 3: Full SARMA(1,1)×(1,1)_12 model\n")
cat("-----------------------------------------------------------\n")

# Use arima.sim for complex model
y3 <- arima.sim(n = 300,
                list(ar = 0.5, ma = 0.3,
                     seasonal = list(sar = 0.6, sma = 0.4, period = 12)),
                innov = rgamma(300, shape = 2, scale = 1) - 2)

# Fit full SARMA model
fit3 <- sarma_pmm2(y3, order = c(1, 1, 1, 1), season = list(period = 12),
                   method = "pmm2", verbose = FALSE)

cat("\nModel summary:\n")
summary(fit3)

cat("\n\nExpected coefficients: AR1 = 0.5, SAR1 = 0.6, MA1 = 0.3, SMA1 = 0.4\n")
cat("PMM2 estimates:\n")
cat("  AR1  =", fit3@coefficients[1], "\n")
cat("  SAR1 =", fit3@coefficients[2], "\n")
cat("  MA1  =", fit3@coefficients[3], "\n")
cat("  SMA1 =", fit3@coefficients[4], "\n\n")

# Test 4: SARIMA(1,1,1)×(1,1,1)_12 with differencing
cat("\n=================================================================\n")
cat("Test 4: SARIMA(1,1,1)×(1,1,1)_12 model with differencing\n")
cat("-----------------------------------------------------------\n")

# Generate non-stationary series
y4 <- arima.sim(n = 300,
                list(order = c(1, 1, 1),
                     ar = 0.4, ma = 0.3,
                     seasonal = list(order = c(1, 1, 1), period = 12,
                                    sar = 0.5, sma = 0.4)),
                innov = rgamma(300, shape = 2, scale = 1) - 2)

# Fit SARIMA model
fit4 <- sarima_pmm2(y4, order = c(1, 1, 1, 1),
                    seasonal = list(order = c(1, 1), period = 12),
                    method = "pmm2", verbose = FALSE)

cat("\nModel summary:\n")
summary(fit4)

cat("\n\nExpected coefficients: AR1 = 0.4, SAR1 = 0.5, MA1 = 0.3, SMA1 = 0.4\n")
cat("PMM2 estimates:\n")
cat("  AR1  =", fit4@coefficients[1], "\n")
cat("  SAR1 =", fit4@coefficients[2], "\n")
cat("  MA1  =", fit4@coefficients[3], "\n")
cat("  SMA1 =", fit4@coefficients[4], "\n\n")

# Test 5: Compare with CSS method
cat("\n=================================================================\n")
cat("Test 5: Comparison PMM2 vs CSS for SARMA(1,0)×(1,1)_12\n")
cat("-----------------------------------------------------------\n")

# Use asymmetric data
y5 <- arima.sim(n = 250,
                list(ar = 0.6,
                     seasonal = list(sar = 0.5, sma = 0.6, period = 12)),
                innov = rgamma(250, shape = 2, scale = 1) - 2)

# Fit with CSS
fit5_css <- sarma_pmm2(y5, order = c(1, 1, 0, 1), season = list(period = 12),
                       method = "css", verbose = FALSE)

# Fit with PMM2
fit5_pmm2 <- sarma_pmm2(y5, order = c(1, 1, 0, 1), season = list(period = 12),
                        method = "pmm2", verbose = FALSE)

cat("\nCSS estimates:\n")
cat("  AR1  =", fit5_css@coefficients[1], "\n")
cat("  SAR1 =", fit5_css@coefficients[2], "\n")
cat("  SMA1 =", fit5_css@coefficients[3], "\n")
cat("  Variance factor g =",
    pmm2_variance_factor(fit5_css@m2, fit5_css@m3, fit5_css@m4)$g, "\n")

cat("\nPMM2 estimates:\n")
cat("  AR1  =", fit5_pmm2@coefficients[1], "\n")
cat("  SAR1 =", fit5_pmm2@coefficients[2], "\n")
cat("  SMA1 =", fit5_pmm2@coefficients[3], "\n")
cat("  Variance factor g =",
    pmm2_variance_factor(fit5_pmm2@m2, fit5_pmm2@m3, fit5_pmm2@m4)$g, "\n")

cat("\n=================================================================\n")
cat("All tests completed successfully!\n")
cat("=================================================================\n")
