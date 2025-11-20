source("R/optimized_direct_pmm2.R")
source("R/sarimax_wrapper.R")
library(numDeriv)

# Test on SARIMA model which failed previously
model_name <- "SARIMA(0,0,1)(0,0,1)[12]"
order <- c(0, 0, 1)
seasonal <- list(order = c(0, 0, 1), period = 12)
params <- c(ma1 = 0.5, sma1 = 0.4, intercept = 0)

n <- 200
set.seed(123)
innovations <- (rchisq(n + 100, df = 3) - 3) / sqrt(6)

theta <- params["ma1"]
Theta <- params["sma1"]
x <- numeric(n + 100)
e <- innovations
for (t in 14:(n + 100)) {
    x[t] <- e[t] + theta * e[t - 1] + Theta * e[t - 12] + theta * Theta * e[t - 13]
}
sim_data <- ts(x[101:(n + 100)])

# MLE
fit_mle <- arima(sim_data, order = order, seasonal = seasonal, include.mean = TRUE)
theta_mle <- coef(fit_mle)

cat("True Params: ma1=0.5, sma1=0.4\n")
cat("MLE Params: ", paste(names(theta_mle), round(theta_mle, 4), sep = "=", collapse = ", "), "\n")

# Optimized Direct PMM2
cat("\nRunning Optimized Direct PMM2...\n")
opt_res <- optimized_direct_pmm2(theta_mle, sim_data, order = order, seasonal = seasonal)
cat("Optimized PMM2 Params: ", paste(round(opt_res$coefficients, 4), collapse = ", "), "\n")

# Calculate MSE for this single run
true_vec <- c(ma1 = 0.5, sma1 = 0.4, intercept = 0)
# Align names
theta_opt <- opt_res$coefficients
mse_mle <- sum((theta_mle - true_vec[names(theta_mle)])^2)
mse_opt <- sum((theta_opt - true_vec[names(theta_opt)])^2)

cat("\nMSE MLE: ", mse_mle, "\n")
cat("MSE Optimized PMM2: ", mse_opt, "\n")
