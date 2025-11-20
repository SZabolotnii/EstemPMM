source("R/unified_pmm2.R")
source("R/sarimax_wrapper.R")

# --- Monte Carlo Simulation Settings ---
R <- 50 # Number of replications (small for speed, increase for final results)
n <- 200 # Sample size

# --- Models to Test ---
models <- list(
    list(
        name = "AR(1)", order = c(1, 0, 0), seasonal = list(order = c(0, 0, 0), period = NA),
        params = c(ar1 = 0.7, intercept = 0)
    ),
    list(
        name = "MA(1)", order = c(0, 0, 1), seasonal = list(order = c(0, 0, 0), period = NA),
        params = c(ma1 = -0.6, intercept = 0)
    ),
    list(
        name = "SARIMA(0,0,1)(0,0,1)[12]", order = c(0, 0, 1), seasonal = list(order = c(0, 0, 1), period = 12),
        params = c(ma1 = 0.5, sma1 = 0.4, intercept = 0)
    )
)

results <- data.frame()

cat("Starting Monte Carlo Simulation (R =", R, ", n =", n, ")...\n\n")

for (model in models) {
    cat("Processing Model:", model$name, "\n")

    mse_mle <- 0
    mse_pmm_onestep <- 0
    mse_pmm_iter <- 0

    for (i in 1:R) {
        set.seed(123 + i)

        # 1. Generate Data (Asymmetric Errors)
        # Using Chi-squared(3) centered and scaled
        innovations <- (rchisq(n + 100, df = 3) - 3) / sqrt(6)

        sim_data <- arima.sim(
            n = n,
            model = list(
                ar = if (model$order[1] > 0) model$params[grep("ar", names(model$params))] else NULL,
                ma = if (model$order[3] > 0) model$params[grep("ma", names(model$params))] else NULL
            ),
            innov = innovations
        )

        # Add seasonal component manually if needed (arima.sim is limited for seasonal)
        # For simplicity in this test, we rely on arima.sim for non-seasonal or simple seasonal if supported.
        # For the SARIMA case, let's construct it manually to be sure.

        if (model$name == "SARIMA(0,0,1)(0,0,1)[12]") {
            # Manual generation for SARIMA
            # x_t = (1 + theta B)(1 + Theta B^12) e_t
            theta <- model$params["ma1"]
            Theta <- model$params["sma1"]
            x <- numeric(n + 100)
            e <- innovations
            for (t in 14:(n + 100)) {
                x[t] <- e[t] + theta * e[t - 1] + Theta * e[t - 12] + theta * Theta * e[t - 13]
            }
            sim_data <- ts(x[101:(n + 100)])
        }

        # 2. MLE Estimation
        fit_mle <- tryCatch(
            {
                arima(sim_data, order = model$order, seasonal = model$seasonal, include.mean = TRUE)
            },
            error = function(e) NULL
        )

        if (is.null(fit_mle)) next

        theta_mle <- coef(fit_mle)

        # 3. PMM2 Estimation
        # Define wrapper closures for this specific dataset
        fn_res <- function(th) get_sarimax_residuals(th, sim_data, order = model$order, seasonal = model$seasonal)
        fn_jac <- function(th) get_sarimax_jacobian(th, sim_data, order = model$order, seasonal = model$seasonal)

        # One-step
        pmm_os <- pmm2_nonlinear_onestep(theta_mle, fn_res, fn_jac)
        theta_os <- pmm_os$coefficients

        # Iterative
        pmm_it <- pmm2_nonlinear_iterative(theta_mle, fn_res, fn_jac, verbose = FALSE)
        theta_it <- pmm_it$coefficients

        # Calculate Squared Errors (sum of squared errors for all params)
        # Align parameters with true values
        # Note: arima() returns intercept last usually.

        # Construct true vector matching MLE output names
        true_vec <- numeric(length(theta_mle))
        names(true_vec) <- names(theta_mle)

        for (nm in names(true_vec)) {
            if (nm %in% names(model$params)) {
                true_vec[nm] <- model$params[nm]
            } else if (nm == "intercept") {
                true_vec[nm] <- 0 # We simulated with 0 mean
            }
        }

        mse_mle <- mse_mle + sum((theta_mle - true_vec)^2)
        mse_pmm_onestep <- mse_pmm_onestep + sum((theta_os - true_vec)^2)
        mse_pmm_iter <- mse_pmm_iter + sum((theta_it - true_vec)^2)
    }

    # Average MSE
    mse_mle <- mse_mle / R
    mse_pmm_onestep <- mse_pmm_onestep / R
    mse_pmm_iter <- mse_pmm_iter / R

    # Store Results
    results <- rbind(results, data.frame(
        Model = model$name,
        MSE_MLE = mse_mle,
        MSE_PMM_OneStep = mse_pmm_onestep,
        MSE_PMM_Iter = mse_pmm_iter,
        Improvement_OneStep = (mse_mle - mse_pmm_onestep) / mse_mle * 100,
        Improvement_Iter = (mse_mle - mse_pmm_iter) / mse_mle * 100
    ))
}

print(results)
