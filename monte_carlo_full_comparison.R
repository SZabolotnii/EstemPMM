source("R/unified_pmm2.R")
source("R/sarimax_wrapper.R")
source("new_sarima/experimental/05_nonlinear_pmm2_for_ma.R")
source("new_sarima/experimental/06_estpmm_style_ma.R")

library(numDeriv)

# --- Monte Carlo Simulation Settings ---
R <- 50
n <- 200

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

cat("Starting Full Monte Carlo Comparison (R =", R, ", n =", n, ")...\n\n")

for (model in models) {
    cat("Processing Model:", model$name, "\n")

    mse_mle <- 0
    mse_pmm_unified_os <- 0
    mse_pmm_unified_it <- 0
    mse_pmm_direct <- 0
    mse_pmm_linear <- 0

    count_linear <- 0 # To track if we actually ran Linearized PMM2

    for (i in 1:R) {
        set.seed(123 + i)

        # 1. Generate Data
        innovations <- (rchisq(n + 100, df = 3) - 3) / sqrt(6)

        if (model$name == "SARIMA(0,0,1)(0,0,1)[12]") {
            theta <- model$params["ma1"]
            Theta <- model$params["sma1"]
            x <- numeric(n + 100)
            e <- innovations
            for (t in 14:(n + 100)) {
                x[t] <- e[t] + theta * e[t - 1] + Theta * e[t - 12] + theta * Theta * e[t - 13]
            }
            sim_data <- ts(x[101:(n + 100)])
        } else {
            sim_data <- arima.sim(
                n = n,
                model = list(
                    ar = if (model$order[1] > 0) model$params[grep("ar", names(model$params))] else NULL,
                    ma = if (model$order[3] > 0) model$params[grep("ma", names(model$params))] else NULL
                ),
                innov = innovations
            )
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

        # Construct true vector
        true_vec <- numeric(length(theta_mle))
        names(true_vec) <- names(theta_mle)
        for (nm in names(true_vec)) {
            if (nm %in% names(model$params)) {
                true_vec[nm] <- model$params[nm]
            } else if (nm == "intercept") true_vec[nm] <- 0
        }

        mse_mle <- mse_mle + sum((theta_mle - true_vec)^2)

        # 3. Unified PMM2
        fn_res <- function(th) get_sarimax_residuals(th, sim_data, order = model$order, seasonal = model$seasonal)
        fn_jac <- function(th) get_sarimax_jacobian(th, sim_data, order = model$order, seasonal = model$seasonal)

        # One-step
        pmm_os <- pmm2_nonlinear_onestep(theta_mle, fn_res, fn_jac)
        mse_pmm_unified_os <- mse_pmm_unified_os + sum((pmm_os$coefficients - true_vec)^2)

        # Iterative
        pmm_it <- pmm2_nonlinear_iterative(theta_mle, fn_res, fn_jac, verbose = FALSE)
        mse_pmm_unified_it <- mse_pmm_unified_it + sum((pmm_it$coefficients - true_vec)^2)

        # 4. Direct Nonlinear PMM2 (Article)
        # Only supports MA/SMA structure easily via the provided script
        # nonlinear_pmm2_ma(x, q, Q, s, ...)
        if (model$name %in% c("MA(1)", "SARIMA(0,0,1)(0,0,1)[12]")) {
            q_val <- model$order[3]
            Q_val <- model$seasonal$order[3]
            s_val <- if (is.na(model$seasonal$period)) 1 else model$seasonal$period

            direct_fit <- nonlinear_pmm2_ma(sim_data,
                q = q_val, Q = Q_val, s = s_val,
                include.mean = TRUE, init_method = "mle", verbose = FALSE
            )

            # Map results to vector
            theta_direct <- theta_mle # Start with MLE structure
            # Fill in
            idx <- 1
            if (q_val > 0) {
                theta_direct[grep("ma", names(theta_direct))] <- direct_fit$ma_coef
            }
            if (Q_val > 0) {
                theta_direct[grep("sma", names(theta_direct))] <- direct_fit$sma_coef
            }
            if ("intercept" %in% names(theta_direct)) {
                theta_direct["intercept"] <- direct_fit$mean
            }

            mse_pmm_direct <- mse_pmm_direct + sum((theta_direct - true_vec)^2)
        } else {
            # For AR(1), Direct PMM2 script doesn't support it (it's MA specific)
            # Use MLE MSE as placeholder to not skew results (or mark as NA)
            mse_pmm_direct <- mse_pmm_direct + sum((theta_mle - true_vec)^2)
        }

        # 5. Linearized PMM2 (EstemPMM-style)
        # Only supports pure MA or pure SMA in the file we sourced
        if (model$name == "MA(1)") {
            lin_fit <- estpmm_style_ma(sim_data, q = 1, include.mean = TRUE, verbose = FALSE)
            theta_lin <- theta_mle
            theta_lin["ma1"] <- lin_fit$ma_coef
            theta_lin["intercept"] <- lin_fit$mean

            mse_pmm_linear <- mse_pmm_linear + sum((theta_lin - true_vec)^2)
            count_linear <- count_linear + 1
        } else {
            # Not supported or mixed
            mse_pmm_linear <- mse_pmm_linear + sum((theta_mle - true_vec)^2)
        }
    }

    # Averages
    mse_mle <- mse_mle / R
    mse_pmm_unified_os <- mse_pmm_unified_os / R
    mse_pmm_unified_it <- mse_pmm_unified_it / R
    mse_pmm_direct <- mse_pmm_direct / R
    mse_pmm_linear <- mse_pmm_linear / R

    # For models where methods weren't applicable, set to NA to avoid confusion
    if (model$name == "AR(1)") {
        mse_pmm_direct <- NA
        mse_pmm_linear <- NA
    }
    if (model$name == "SARIMA(0,0,1)(0,0,1)[12]") {
        mse_pmm_linear <- NA # Mixed not supported by simple script
    }

    results <- rbind(results, data.frame(
        Model = model$name,
        MLE = mse_mle,
        Unified_OneStep = mse_pmm_unified_os,
        Unified_Iter = mse_pmm_unified_it,
        Direct_Nonlinear = mse_pmm_direct,
        Linearized_Specific = mse_pmm_linear
    ))
}

print(results)
