#!/usr/bin/env Rscript
# ==============================================================================
# Comprehensive Monte Carlo Testing for EstemPMM Seasonal Models
# Testing: SAR, SMA, SARMA, SARIMA with Gaussian and Asymmetric Innovations
# Replications: 100, 200, 500
# Author: Serhii Zabolotnii
# Date: 2025-11-20
# ==============================================================================

suppressPackageStartupMessages({
    library(EstemPMM)
    library(parallel)
    library(ggplot2)
    library(tidyr)
    library(dplyr)
})

# ==============================================================================
# Configuration
# ==============================================================================

CONFIG <- list(
    n_replications = c(100, 200, 500),  # Different MC replication counts
    sample_sizes = c(100, 200, 500),     # Sample sizes to test
    seasonal_period = 12,
    n_cores = max(1, detectCores() - 1),
    output_dir = "mc_results_final_comprehensive",
    seed = 20251120
)

cat("==============================================================================\n")
cat("Comprehensive Monte Carlo Testing - EstemPMM Seasonal Models\n")
cat("==============================================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Cores available:", CONFIG$n_cores, "\n")
cat("Sample sizes:", paste(CONFIG$sample_sizes, collapse = ", "), "\n")
cat("Replication counts:", paste(CONFIG$n_replications, collapse = ", "), "\n")
cat("Seasonal period:", CONFIG$seasonal_period, "\n")
cat("==============================================================================\n\n")

# Create output directory
if (!dir.exists(CONFIG$output_dir)) {
    dir.create(CONFIG$output_dir, recursive = TRUE)
    cat("Created output directory:", CONFIG$output_dir, "\n\n")
}

set.seed(CONFIG$seed)

# ==============================================================================
# Innovation Generation Functions
# ==============================================================================

#' Generate innovations with different distributions
#' @param n Sample size
#' @param type Innovation type: "gaussian", "gamma", "lognormal", "exponential"
#' @return Numeric vector of innovations (centered to have mean 0)
generate_innovations <- function(n, type = "gaussian") {
    innovations <- switch(type,
        "gaussian" = rnorm(n, mean = 0, sd = 1),
        "gamma" = {
            # Gamma(shape=2, scale=1) has mean=2, center to 0
            rgamma(n, shape = 2, scale = 1) - 2
        },
        "lognormal" = {
            # LogNormal(0, 0.5) centered
            rlnorm(n, meanlog = 0, sdlog = 0.5) - exp(0.125)
        },
        "exponential" = {
            # Exponential(rate=1) has mean=1, center to 0
            rexp(n, rate = 1) - 1
        },
        stop("Unknown innovation type: ", type)
    )
    
    return(innovations)
}

# ==============================================================================
# Time Series Generation Functions
# ==============================================================================

#' Generate SAR(p,P)_s time series
#' @param n Sample size
#' @param ar Non-seasonal AR coefficients
#' @param sar Seasonal AR coefficients
#' @param s Seasonal period
#' @param innov_type Innovation distribution type
#' @return Time series vector
generate_sar_series <- function(n, ar = NULL, sar = NULL, s = 12, 
                                innov_type = "gaussian") {
    p <- length(ar)
    P <- length(sar)
    n_burn <- 200
    n_total <- n + n_burn
    
    innov <- generate_innovations(n_total, innov_type)
    y <- numeric(n_total)
    
    max_lag <- max(p, P * s)
    
    for (t in (max_lag + 1):n_total) {
        ar_term <- if (p > 0) sum(ar * y[(t - 1):(t - p)]) else 0
        sar_term <- if (P > 0) sum(sar * y[t - s * (1:P)]) else 0
        y[t] <- ar_term + sar_term + innov[t]
    }
    
    return(y[(n_burn + 1):n_total])
}

#' Generate SMA(Q)_s time series
#' @param n Sample size
#' @param sma Seasonal MA coefficients
#' @param s Seasonal period
#' @param innov_type Innovation distribution type
#' @return Time series vector
generate_sma_series <- function(n, sma = NULL, s = 12, innov_type = "gaussian") {
    Q <- length(sma)
    n_burn <- 200
    n_total <- n + n_burn
    
    innov <- generate_innovations(n_total, innov_type)
    y <- numeric(n_total)
    eps <- numeric(n_total)
    
    for (t in 1:n_total) {
        sma_term <- if (t > s * Q && Q > 0) {
            sum(sma * eps[t - s * (1:Q)])
        } else {
            0
        }
        eps[t] <- innov[t]
        y[t] <- sma_term + eps[t]
    }
    
    return(y[(n_burn + 1):n_total])
}

#' Generate SARMA(p,P,q,Q)_s time series
#' @param n Sample size
#' @param ar Non-seasonal AR coefficients
#' @param sar Seasonal AR coefficients
#' @param ma Non-seasonal MA coefficients (not used in this simple version)
#' @param sma Seasonal MA coefficients
#' @param s Seasonal period
#' @param innov_type Innovation distribution type
#' @return Time series vector
generate_sarma_series <- function(n, ar = NULL, sar = NULL, ma = NULL, 
                                  sma = NULL, s = 12, innov_type = "gaussian") {
    n_burn <- 200
    n_total <- n + n_burn
    
    innov <- generate_innovations(n_total, innov_type)
    
    model_list <- list()
    if (!is.null(ar)) model_list$ar <- ar
    if (!is.null(ma)) model_list$ma <- ma
    
    if (!is.null(sar) || !is.null(sma)) {
        seasonal_list <- list(period = s)
        if (!is.null(sar)) seasonal_list$sar <- sar
        if (!is.null(sma)) seasonal_list$sma <- sma
        model_list$seasonal <- seasonal_list
    }
    
    y <- arima.sim(n = n_total, model = model_list, innov = innov)
    return(as.numeric(y[(n_burn + 1):n_total]))
}

#' Generate SARIMA(p,d,q)×(P,D,Q)_s time series
#' @param n Sample size
#' @param ar Non-seasonal AR coefficients
#' @param d Non-seasonal differencing order
#' @param sar Seasonal AR coefficients
#' @param D Seasonal differencing order
#' @param sma Seasonal MA coefficients
#' @param s Seasonal period
#' @param innov_type Innovation distribution type
#' @return Time series vector
generate_sarima_series <- function(n, ar = NULL, d = 1, sar = NULL, D = 1, 
                                   sma = NULL, s = 12, innov_type = "gaussian") {
    n_burn <- 200
    n_total <- n + n_burn
    
    innov <- generate_innovations(n_total, innov_type)
    
    model_list <- list(
        order = c(length(ar), d, 0),
        ar = ar
    )
    
    seasonal_list <- list(
        order = c(length(sar), D, length(sma)),
        period = s
    )
    if (!is.null(sar)) seasonal_list$sar <- sar
    if (!is.null(sma)) seasonal_list$sma <- sma
    
    model_list$seasonal <- seasonal_list
    
    y <- arima.sim(n = n_total, model = model_list, innov = innov)
    return(as.numeric(y[(n_burn + 1):n_total]))
}

# ==============================================================================
# Model Fitting Functions
# ==============================================================================

#' Safe model fitting wrapper
#' @param fit_func Fitting function (e.g., sar_pmm2, sma_pmm2)
#' @param ... Arguments to pass to fit_func
#' @return List with coefficients, convergence status, moments, and error info
safe_fit <- function(fit_func, ...) {
    result <- tryCatch({
        fit <- fit_func(...)
        list(
            coefficients = fit@coefficients,
            convergence = fit@convergence,
            m2 = fit@m2,
            m3 = fit@m3,
            m4 = fit@m4,
            residuals = fit@residuals,
            error = NULL,
            success = TRUE
        )
    }, error = function(e) {
        list(
            coefficients = NA,
            convergence = FALSE,
            m2 = NA, m3 = NA, m4 = NA,
            residuals = NA,
            error = e$message,
            success = FALSE
        )
    })
    
    return(result)
}

# ==============================================================================
# Metrics Computation
# ==============================================================================

#' Compute estimation metrics
#' @param estimates Vector of parameter estimates
#' @param true_value True parameter value
#' @return List with bias, RMSE, variance, MAE
compute_metrics <- function(estimates, true_value) {
    valid_estimates <- estimates[!is.na(estimates)]
    
    if (length(valid_estimates) == 0) {
        return(list(bias = NA, rmse = NA, variance = NA, mae = NA, n_valid = 0))
    }
    
    bias <- mean(valid_estimates) - true_value
    rmse <- sqrt(mean((valid_estimates - true_value)^2))
    variance <- var(valid_estimates)
    mae <- mean(abs(valid_estimates - true_value))
    
    list(
        bias = bias,
        rmse = rmse,
        variance = variance,
        mae = mae,
        n_valid = length(valid_estimates)
    )
}

# ==============================================================================
# Scenario Definitions
# ==============================================================================

SCENARIOS <- list(
    # SAR(1,1)_12 Model
    list(
        name = "SAR_1_1_12",
        type = "SAR",
        true_params = list(ar1 = 0.5, sar1 = 0.6),
        param_names = c("ar1", "sar1"),
        generator = function(n, innov_type) {
            generate_sar_series(n, ar = 0.5, sar = 0.6, s = 12, 
                              innov_type = innov_type)
        },
        fit_func = function(y, method) {
            sar_pmm2(y, order = c(1, 1), 
                    season = list(period = 12),
                    method = method, verbose = FALSE)
        }
    ),
    
    # SMA(1)_12 Model
    list(
        name = "SMA_1_12",
        type = "SMA",
        true_params = list(sma1 = 0.6),
        param_names = c("sma1"),
        generator = function(n, innov_type) {
            generate_sma_series(n, sma = 0.6, s = 12, innov_type = innov_type)
        },
        fit_func = function(y, method) {
            sma_pmm2(y, order = 1, 
                    season = list(period = 12),
                    method = method, verbose = FALSE)
        }
    ),
    
    # SARMA(1,0)×(1,1)_12 Model
    list(
        name = "SARMA_1_0_1_1_12",
        type = "SARMA",
        true_params = list(ar1 = 0.5, sar1 = 0.6, sma1 = 0.4),
        param_names = c("ar1", "sar1", "sma1"),
        generator = function(n, innov_type) {
            generate_sarma_series(n, ar = 0.5, sar = 0.6, sma = 0.4, 
                                s = 12, innov_type = innov_type)
        },
        fit_func = function(y, method) {
            sarma_pmm2(y, order = c(1, 1, 0, 1), 
                      season = list(period = 12),
                      method = method, verbose = FALSE)
        }
    ),
    
    # SARIMA(1,1,0)×(1,1,1)_12 Model
    list(
        name = "SARIMA_1_1_0_1_1_1_12",
        type = "SARIMA",
        true_params = list(ar1 = 0.4, sar1 = 0.5, sma1 = 0.6),
        param_names = c("ar1", "sar1", "sma1"),
        generator = function(n, innov_type) {
            generate_sarima_series(n, ar = 0.4, d = 1, sar = 0.5, 
                                 D = 1, sma = 0.6, s = 12, 
                                 innov_type = innov_type)
        },
        fit_func = function(y, method) {
            sarima_pmm2(y, order = c(1, 1, 0, 1), 
                       seasonal = list(order = c(1, 1, 1), period = 12),
                       method = method, verbose = FALSE)
        }
    )
)

INNOVATION_TYPES <- c("gaussian", "gamma", "lognormal", "exponential")
METHODS <- c("pmm2", "css")  # Compare PMM2 against CSS

# ==============================================================================
# Single Replication Worker
# ==============================================================================

#' Run a single Monte Carlo replication
#' @param scenario Scenario definition
#' @param n Sample size
#' @param innov_type Innovation type
#' @param rep_id Replication ID
#' @return List with results for all methods
run_single_replication <- function(scenario, n, innov_type, rep_id) {
    # Generate data
    y <- scenario$generator(n, innov_type)
    
    # Fit with different methods
    results <- list()
    
    for (method in METHODS) {
        fit_result <- safe_fit(scenario$fit_func, y = y, method = method)
        
        # Extract g-factor for PMM2
        g_factor <- if (method == "pmm2" && fit_result$success) {
            tryCatch({
                pmm2_variance_factor(fit_result$m2, fit_result$m3, 
                                    fit_result$m4)$g
            }, error = function(e) NA)
        } else {
            NA
        }
        
        results[[method]] <- list(
            coefficients = fit_result$coefficients,
            convergence = fit_result$convergence,
            g_factor = g_factor,
            success = fit_result$success
        )
    }
    
    return(results)
}

# ==============================================================================
# Main Simulation Loop
# ==============================================================================

#' Run complete Monte Carlo simulation for one configuration
#' @param scenario Scenario definition
#' @param n Sample size
#' @param n_reps Number of replications
#' @param innov_type Innovation type
#' @return Data frame with all results
run_monte_carlo <- function(scenario, n, n_reps, innov_type) {
    cat(sprintf("  Running %s, n=%d, reps=%d, innov=%s\n", 
                scenario$name, n, n_reps, innov_type))
    
    # Run replications in parallel
    results_list <- mclapply(1:n_reps, function(i) {
        run_single_replication(scenario, n, innov_type, i)
    }, mc.cores = CONFIG$n_cores)
    
    # Organize results by method and parameter
    n_params <- length(scenario$param_names)
    
    organized_results <- list()
    
    for (method in METHODS) {
        method_coefs <- matrix(NA, nrow = n_reps, ncol = n_params)
        method_converged <- logical(n_reps)
        method_g <- numeric(n_reps)
        
        for (i in 1:n_reps) {
            res <- results_list[[i]][[method]]
            if (res$success && length(res$coefficients) == n_params) {
                method_coefs[i, ] <- res$coefficients
                method_converged[i] <- res$convergence
                method_g[i] <- res$g_factor
            }
        }
        
        # Compute metrics for each parameter
        param_metrics <- list()
        for (j in 1:n_params) {
            param_name <- scenario$param_names[j]
            true_val <- scenario$true_params[[param_name]]
            param_metrics[[param_name]] <- compute_metrics(method_coefs[, j], true_val)
        }
        
        organized_results[[method]] <- list(
            metrics = param_metrics,
            convergence_rate = mean(method_converged),
            mean_g = if (method == "pmm2") mean(method_g, na.rm = TRUE) else NA,
            coefficients = method_coefs
        )
    }
    
    return(organized_results)
}

# ==============================================================================
# Execute All Simulations
# ==============================================================================

all_results <- list()
timestamp_start <- Sys.time()

for (scenario in SCENARIOS) {
    cat("\n")
    cat("==============================================================================\n")
    cat("Scenario:", scenario$name, "\n")
    cat("==============================================================================\n")
    
    scenario_results <- list()
    
    for (innov_type in INNOVATION_TYPES) {
        cat("\nInnovation type:", innov_type, "\n")
        
        innov_results <- list()
        
        for (n in CONFIG$sample_sizes) {
            for (n_reps in CONFIG$n_replications) {
                config_name <- sprintf("n%d_reps%d", n, n_reps)
                
                mc_result <- run_monte_carlo(scenario, n, n_reps, innov_type)
                
                innov_results[[config_name]] <- mc_result
                
                # Print summary
                cat(sprintf("\n  Results for n=%d, reps=%d:\n", n, n_reps))
                for (method in METHODS) {
                    cat(sprintf("    %s:\n", toupper(method)))
                    res <- mc_result[[method]]
                    cat(sprintf("      Convergence rate: %.3f\n", 
                              res$convergence_rate))
                    if (method == "pmm2") {
                        cat(sprintf("      Mean g-factor: %.4f\n", res$mean_g))
                    }
                    
                    for (param_name in scenario$param_names) {
                        metrics <- res$metrics[[param_name]]
                        cat(sprintf("      %s - RMSE: %.4f, Bias: %.4f, Var: %.4f\n",
                                  param_name, metrics$rmse, metrics$bias, 
                                  metrics$variance))
                    }
                }
            }
        }
        
        scenario_results[[innov_type]] <- innov_results
    }
    
    all_results[[scenario$name]] <- scenario_results
}

timestamp_end <- Sys.time()
elapsed_time <- difftime(timestamp_end, timestamp_start, units = "mins")

# ==============================================================================
# Save Results
# ==============================================================================

cat("\n")
cat("==============================================================================\n")
cat("Saving Results\n")
cat("==============================================================================\n")

# Save full results
results_file <- file.path(CONFIG$output_dir, "comprehensive_mc_results.rds")
saveRDS(list(
    results = all_results,
    config = CONFIG,
    scenarios = SCENARIOS,
    timestamp_start = timestamp_start,
    timestamp_end = timestamp_end,
    elapsed_time = elapsed_time
), results_file)

cat("Full results saved to:", results_file, "\n")

# ==============================================================================
# Generate Summary Report
# ==============================================================================

cat("\n")
cat("==============================================================================\n")
cat("SUMMARY REPORT\n")
cat("==============================================================================\n")
cat("Execution time:", round(elapsed_time, 2), "minutes\n")
cat("Total scenarios:", length(SCENARIOS), "\n")
cat("Innovation types:", length(INNOVATION_TYPES), "\n")
cat("Sample sizes:", length(CONFIG$sample_sizes), "\n")
cat("Replication counts:", length(CONFIG$n_replications), "\n\n")

# Create summary table
summary_data <- list()

for (scenario_name in names(all_results)) {
    scenario_res <- all_results[[scenario_name]]
    
    for (innov_type in names(scenario_res)) {
        innov_res <- scenario_res[[innov_type]]
        
        for (config_name in names(innov_res)) {
            config_res <- innov_res[[config_name]]
            
            # Extract n and n_reps from config_name
            parts <- strsplit(config_name, "_")[[1]]
            n <- as.numeric(gsub("n", "", parts[1]))
            n_reps <- as.numeric(gsub("reps", "", parts[2]))
            
            # Compare PMM2 vs CSS
            pmm2_res <- config_res[["pmm2"]]
            css_res <- config_res[["css"]]
            
            # Get first parameter for comparison
            first_param <- names(pmm2_res$metrics)[1]
            
            pmm2_rmse <- pmm2_res$metrics[[first_param]]$rmse
            css_rmse <- css_res$metrics[[first_param]]$rmse
            rmse_improvement <- (css_rmse - pmm2_rmse) / css_rmse * 100
            
            pmm2_var <- pmm2_res$metrics[[first_param]]$variance
            css_var <- css_res$metrics[[first_param]]$variance
            var_reduction <- (css_var - pmm2_var) / css_var * 100
            
            summary_data[[length(summary_data) + 1]] <- list(
                Scenario = scenario_name,
                Innovation = innov_type,
                n = n,
                Reps = n_reps,
                Parameter = first_param,
                PMM2_RMSE = pmm2_rmse,
                CSS_RMSE = css_rmse,
                RMSE_Improvement = rmse_improvement,
                Var_Reduction = var_reduction,
                PMM2_Conv = pmm2_res$convergence_rate,
                CSS_Conv = css_res$convergence_rate,
                Mean_g = pmm2_res$mean_g
            )
        }
    }
}

# Convert to data frame
summary_df <- do.call(rbind, lapply(summary_data, as.data.frame))

# Save summary
summary_file <- file.path(CONFIG$output_dir, "summary_report.csv")
write.csv(summary_df, summary_file, row.names = FALSE)
cat("\nSummary table saved to:", summary_file, "\n")

# Print key findings
cat("\n")
cat("==============================================================================\n")
cat("KEY FINDINGS\n")
cat("==============================================================================\n\n")

# Average improvements by scenario
for (scenario_name in unique(summary_df$Scenario)) {
    scenario_subset <- summary_df[summary_df$Scenario == scenario_name, ]
    avg_rmse_imp <- mean(scenario_subset$RMSE_Improvement, na.rm = TRUE)
    avg_var_red <- mean(scenario_subset$Var_Reduction, na.rm = TRUE)
    
    cat(sprintf("%s:\n", scenario_name))
    cat(sprintf("  Average RMSE improvement: %.2f%%\n", avg_rmse_imp))
    cat(sprintf("  Average variance reduction: %.2f%%\n", avg_var_red))
    cat(sprintf("  PMM2 convergence: %.3f\n", 
                mean(scenario_subset$PMM2_Conv, na.rm = TRUE)))
    cat(sprintf("  CSS convergence: %.3f\n", 
                mean(scenario_subset$CSS_Conv, na.rm = TRUE)))
    cat("\n")
}

cat("==============================================================================\n")
cat("Testing complete!\n")
cat("==============================================================================\n")
