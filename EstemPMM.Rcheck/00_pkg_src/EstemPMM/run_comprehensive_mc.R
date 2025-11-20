#!/usr/bin/env Rscript

# ==============================================================================
# Comprehensive Monte Carlo Simulation for EstemPMM
# Comparing PMM2 vs CSS/ML for Seasonal Models
# ==============================================================================

suppressPackageStartupMessages({
    library(EstemPMM)
    library(parallel)
    library(methods)
})

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

CONFIG <- list(
    n_replications = 1000,
    sample_sizes = c(100, 200, 500, 1000),
    seasonal_period = 12,
    n_cores = max(1, parallel::detectCores() - 1),
    output_dir = "mc_results_comprehensive",
    seed = 20250115
)

# Create output directory
if (!dir.exists(CONFIG$output_dir)) {
    dir.create(CONFIG$output_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# Data Generation Functions
# ------------------------------------------------------------------------------

generate_innovations <- function(n, type) {
    switch(type,
        "normal" = rnorm(n),
        "gamma" = rgamma(n, shape = 2, scale = 1) - 2, # Centered Gamma(2,1)
        "lognormal" = rlnorm(n, meanlog = 0, sdlog = 0.5) - exp(0.125), # Centered LogNormal
        stop("Unknown innovation type")
    )
}

generate_series <- function(model_type, n, params, innov_type, s = 12) {
    n_burn <- 200
    n_total <- n + n_burn
    innov <- generate_innovations(n_total, innov_type)

    y <- switch(model_type,
        "SAR" = {
            # SAR(1,1)_12: y_t = phi*y_{t-1} + Phi*y_{t-s} + eps_t
            # Note: This is a simplified additive SAR, not multiplicative
            y_sim <- numeric(n_total)
            phi <- params$ar
            Phi <- params$sar
            for (t in (s + 1):n_total) {
                y_sim[t] <- phi * y_sim[t - 1] + Phi * y_sim[t - s] + innov[t]
            }
            y_sim
        },
        "SMA" = {
            # SMA(1)_12: y_t = eps_t + Theta*eps_{t-s}
            y_sim <- numeric(n_total)
            Theta <- params$sma
            for (t in (s + 1):n_total) {
                y_sim[t] <- innov[t] + Theta * innov[t - s]
            }
            y_sim
        },
        "SARMA" = {
            # SARMA(1,0)x(1,1)_12
            model_list <- list(
                ar = params$ar,
                seasonal = list(period = s, sar = params$sar, sma = params$sma)
            )
            arima.sim(n = n_total, model = model_list, innov = innov)
        },
        "SARIMA" = {
            # SARIMA(1,1,0)x(0,1,1)_12
            model_list <- list(
                order = c(1, 1, 0),
                ar = params$ar,
                seasonal = list(order = c(0, 1, 1), period = s, sma = params$sma)
            )
            arima.sim(n = n_total, model = model_list, innov = innov)
        }
    )

    return(as.numeric(tail(y, n)))
}

# ------------------------------------------------------------------------------
# Scenarios
# ------------------------------------------------------------------------------

SCENARIOS <- list(
    list(
        name = "SAR_1_1",
        type = "SAR",
        params = list(ar = 0.5, sar = 0.6),
        fit_args = list(order = c(1, 1), season = list(period = 12))
    ),
    list(
        name = "SMA_1",
        type = "SMA",
        params = list(sma = 0.6),
        fit_args = list(order = 1, season = list(period = 12))
    ),
    list(
        name = "SARMA_1_0_1_1",
        type = "SARMA",
        params = list(ar = 0.5, sar = 0.6, sma = 0.4),
        fit_args = list(order = c(1, 1, 0, 1), season = list(period = 12))
    ),
    list(
        name = "SARIMA_1_1_0_0_1_1",
        type = "SARIMA",
        params = list(ar = 0.4, sma = 0.6),
        fit_args = list(order = c(1, 1, 0, 0), seasonal = list(order = c(0, 1, 1), period = 12))
    )
)

INNOVATION_TYPES <- c("normal", "gamma", "lognormal")

# ------------------------------------------------------------------------------
# Worker Function
# ------------------------------------------------------------------------------

run_replication <- function(rep_id, scenario, n, innov_type, config) {
    # Generate data
    y <- generate_series(scenario$type, n, scenario$params, innov_type, config$seasonal_period)

    # Select fitting function based on scenario type
    fit_func <- switch(scenario$type,
        "SAR" = EstemPMM::sar_pmm2,
        "SMA" = EstemPMM::sma_pmm2,
        "SARMA" = EstemPMM::sarma_pmm2,
        "SARIMA" = EstemPMM::sarima_pmm2
    )

    # Fit PMM2
    res_pmm2 <- tryCatch(
        {
            fit <- do.call(fit_func, c(list(x = y, method = "pmm2", verbose = FALSE), scenario$fit_args))
            list(coef = fit@coefficients, converged = fit@convergence, error = NULL)
        },
        error = function(e) list(coef = NULL, converged = FALSE, error = e$message)
    )

    # Fit CSS (Control)
    res_css <- tryCatch(
        {
            fit <- do.call(fit_func, c(list(x = y, method = "css", verbose = FALSE), scenario$fit_args))
            list(coef = fit@coefficients, converged = fit@convergence, error = NULL)
        },
        error = function(e) list(coef = NULL, converged = FALSE, error = e$message)
    )

    return(list(rep_id = rep_id, pmm2 = res_pmm2, css = res_css))
}

# ------------------------------------------------------------------------------
# Main Loop
# ------------------------------------------------------------------------------

cat("Starting Comprehensive Monte Carlo Simulation\n")
cat(sprintf("Cores: %d, Replications: %d\n", CONFIG$n_cores, CONFIG$n_replications))

for (scenario in SCENARIOS) {
    cat(sprintf("\nScenario: %s\n", scenario$name))

    for (innov in INNOVATION_TYPES) {
        cat(sprintf("  Innovation: %s\n", innov))

        for (n in CONFIG$sample_sizes) {
            cat(sprintf("    Sample size: %d... ", n))

            # Run parallel simulation
            results <- mclapply(1:CONFIG$n_replications, function(i) {
                run_replication(i, scenario, n, innov, CONFIG)
            }, mc.cores = CONFIG$n_cores)

            # Save intermediate results
            filename <- sprintf("%s/%s_%s_n%d.rds", CONFIG$output_dir, scenario$name, innov, n)
            saveRDS(results, filename)
            cat("Done. Saved to", filename, "\n")
        }
    }
}

cat("\nSimulation Complete!\n")
