# Compare PMM2 Strategies for SMA and Mixed Models
# 1. MLE (Baseline)
# 2. Linear PMM2 (EstemPMM-style via sarima_pmm2)

set.seed(123)
library(EstemPMM)

# Data Generation: Additive MA+SMA
generate_additive_series <- function(n, theta, Theta, s, innov) {
    x <- numeric(n)
    q <- length(theta)
    Q <- length(Theta)

    max_lag <- max(q, Q * s)
    x[1:max_lag] <- innov[1:max_lag]

    for (t in (max_lag + 1):n) {
        ma_term <- if (q > 0) sum(theta * innov[(t - 1):(t - q)]) else 0

        sma_term <- 0
        if (Q > 0) {
            for (J in 1:Q) {
                lag <- J * s
                sma_term <- sma_term + Theta[J] * innov[t - lag]
            }
        }
        x[t] <- innov[t] + ma_term + sma_term
    }
    return(x)
}

# Data Generation: Multiplicative SARIMA
generate_multiplicative_series <- function(n, theta, Theta, s, innov) {
    x <- numeric(n)
    max_lag <- s + 1
    x[1:max_lag] <- innov[1:max_lag]

    for (t in (max_lag + 1):n) {
        term1 <- theta * innov[t - 1]
        term2 <- Theta * innov[t - s]
        term3 <- theta * Theta * innov[t - s - 1]
        x[t] <- innov[t] + term1 + term2 + term3
    }
    return(x)
}

# Simulation Parameters
n <- 100
s <- 12
theta <- 0.4
Theta <- 0.6
n_sim <- 50

# ---------------------------------------------------------
# Experiment 1: Mixed MA(1)+SMA(1) - Additive Data
# ---------------------------------------------------------
cat("\n=== Experiment 1: Mixed MA(1)+SMA(1) - Additive Data ===\n")
mse_mle <- numeric(n_sim)
mse_lin <- numeric(n_sim)

for (i in 1:n_sim) {
    innov <- rgamma(n, shape = 2, scale = 1) - 2
    x <- generate_additive_series(n, theta, Theta, s, innov)

    # MLE
    fit_mle <- tryCatch(
        {
            stats::arima(x, order = c(0, 0, 1), seasonal = list(order = c(0, 0, 1), period = s), include.mean = FALSE, method = "CSS-ML")
        },
        error = function(e) NULL
    )

    if (!is.null(fit_mle)) {
        mse_mle[i] <- (coef(fit_mle)["ma1"] - theta)^2 + (coef(fit_mle)["sma1"] - Theta)^2
    } else {
        mse_mle[i] <- NA
    }

    # Linear PMM2 (via sarima_pmm2)
    fit_lin <- tryCatch(
        {
            sarima_pmm2(x,
                order = c(0, 0, 1, 1), seasonal = list(order = c(0, 0), period = s),
                ma_method = "pmm2", include.mean = FALSE, verbose = FALSE, multiplicative = TRUE
            )
        },
        error = function(e) NULL
    )

    if (!is.null(fit_lin) && fit_lin@convergence) {
        coefs <- fit_lin@coefficients
        mse_lin[i] <- (coefs[1] - theta)^2 + (coefs[2] - Theta)^2
    } else {
        mse_lin[i] <- NA
    }
}

cat(sprintf("MLE MSE: %.6f\n", mean(mse_mle, na.rm = TRUE)))
cat(sprintf("Lin PMM2 MSE: %.6f\n", mean(mse_lin, na.rm = TRUE)))


# ---------------------------------------------------------
# Experiment 2: Mixed MA(1)+SMA(1) - Multiplicative Data
# ---------------------------------------------------------
cat("\n=== Experiment 2: Mixed MA(1)+SMA(1) - Multiplicative Data ===\n")
mse_mle <- numeric(n_sim)
mse_lin <- numeric(n_sim)

for (i in 1:n_sim) {
    innov <- rgamma(n, shape = 2, scale = 1) - 2
    x <- generate_multiplicative_series(n, theta, Theta, s, innov)

    # MLE
    fit_mle <- tryCatch(
        {
            stats::arima(x, order = c(0, 0, 1), seasonal = list(order = c(0, 0, 1), period = s), include.mean = FALSE, method = "CSS-ML")
        },
        error = function(e) NULL
    )

    if (!is.null(fit_mle)) {
        mse_mle[i] <- (coef(fit_mle)["ma1"] - theta)^2 + (coef(fit_mle)["sma1"] - Theta)^2
    } else {
        mse_mle[i] <- NA
    }

    # Linear PMM2
    fit_lin <- tryCatch(
        {
            sarima_pmm2(x,
                order = c(0, 0, 1, 1), seasonal = list(order = c(0, 0), period = s),
                ma_method = "pmm2", include.mean = FALSE, verbose = FALSE, multiplicative = TRUE
            )
        },
        error = function(e) NULL
    )

    if (!is.null(fit_lin) && fit_lin@convergence) {
        coefs <- fit_lin@coefficients
        mse_lin[i] <- (coefs[1] - theta)^2 + (coefs[2] - Theta)^2
    } else {
        mse_lin[i] <- NA
    }
}

cat(sprintf("MLE MSE: %.6f\n", mean(mse_mle, na.rm = TRUE)))
cat(sprintf("Lin PMM2 MSE: %.6f\n", mean(mse_lin, na.rm = TRUE)))

# ---------------------------------------------------------
# Experiment 3: Full SARIMA(1,0,1)x(1,0,1)_12
# ---------------------------------------------------------
cat("\n=== Experiment 3: Full SARIMA(1,0,1)x(1,0,1)_12 ===\n")

# Parameters
phi <- 0.5
Phi <- 0.3
theta <- 0.4
Theta <- 0.4
s <- 12
n <- 200 # Increased to 200
R <- 50

# Expand polynomials for simulation
# AR: (1 - phi*B)(1 - Phi*B^s) = 1 - phi*B - Phi*B^s + phi*Phi*B^{s+1}
ar_poly <- numeric(s + 2)
ar_poly[2] <- phi
ar_poly[s + 1] <- Phi
ar_poly[s + 2] <- -phi * Phi
ar_coefs <- ar_poly[-1]

# MA: (1 + theta*B)(1 + Theta*B^s) = 1 + theta*B + Theta*B^s + theta*Theta*B^{s+1}
ma_poly <- numeric(s + 2)
ma_poly[2] <- theta
ma_poly[s + 1] <- Theta
ma_poly[s + 2] <- theta * Theta
ma_coefs <- ma_poly[-1]

mse_mle_3 <- 0
mse_pmm2_3 <- 0

for (i in 1:R) {
    set.seed(123 + i)
    innov <- rgamma(n, shape = 2, scale = 1) - 2

    # Simulate SARIMA
    x <- arima.sim(n = n, list(ar = ar_coefs, ma = ma_coefs), innov = innov)
    x <- x + 10 # Add mean

    # 1. MLE
    fit_mle <- tryCatch(
        {
            stats::arima(x, order = c(1, 0, 1), seasonal = list(order = c(1, 0, 1), period = s))
        },
        error = function(e) NULL
    )

    if (!is.null(fit_mle)) {
        # Check stationarity/invertibility
        if (abs(fit_mle$coef["ar1"]) < 1 && abs(fit_mle$coef["sar1"]) < 1) {
            mse_mle_3 <- mse_mle_3 + mean(fit_mle$residuals^2)
        } else {
            mse_mle_3 <- mse_mle_3 + var(x)
        }
    } else {
        mse_mle_3 <- mse_mle_3 + var(x)
    }

    # 2. PMM2 (Full SARIMA)
    fit_pmm2 <- tryCatch(
        {
            EstemPMM::sarima_pmm2(x,
                order = c(1, 1, 1, 1), seasonal = list(order = c(0, 0), period = s),
                multiplicative = TRUE, verbose = FALSE
            )
        },
        error = function(e) {
            print(e)
            NULL
        }
    )

    if (i == 1) {
        cat("True Coefs: ar1=", phi, " sar1=", Phi, " ma1=", theta, " sma1=", Theta, "\n")
        if (!is.null(fit_pmm2)) {
            cat("PMM2 Coefs:\n")
            print(fit_pmm2@coefficients)
        }
    }

    if (!is.null(fit_pmm2)) {
        mse_pmm2_3 <- mse_pmm2_3 + mean(fit_pmm2@residuals^2)
    } else {
        mse_pmm2_3 <- mse_pmm2_3 + var(x)
    }
}

cat(sprintf("MLE MSE: %.6f\n", mse_mle_3 / R))
cat(sprintf("Lin PMM2 MSE: %.6f\n", mse_pmm2_3 / R))
