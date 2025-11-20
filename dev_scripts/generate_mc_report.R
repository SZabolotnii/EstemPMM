#!/usr/bin/env Rscript

# ==============================================================================
# Generate Monte Carlo Report
# Analyzes results from run_comprehensive_mc.R
# ==============================================================================

suppressPackageStartupMessages({
    library(methods)
})

# Configuration
RESULTS_DIR <- "mc_results_comprehensive"
OUTPUT_FILE <- "monte_carlo_report.md"

# Helper function to compute metrics
compute_metrics <- function(estimates, true_val) {
    estimates <- estimates[!is.na(estimates)]
    if (length(estimates) == 0) {
        return(list(bias = NA, mse = NA, var = NA))
    }

    bias <- mean(estimates - true_val)
    var_est <- var(estimates)
    mse <- bias^2 + var_est

    list(bias = bias, mse = mse, var = var_est)
}

# Main processing
cat("Generating Monte Carlo Report...\n")

files <- list.files(RESULTS_DIR, pattern = "\\.rds$", full.names = TRUE)
if (length(files) == 0) {
    stop("No result files found in ", RESULTS_DIR)
}

# Initialize report
report_lines <- c(
    "# Comprehensive Monte Carlo Simulation Report",
    "",
    sprintf("Generated: %s", Sys.time()),
    "",
    "## Summary of Variance Reduction (PMM2 vs CSS)",
    "",
    "| Scenario | Innovation | Sample Size | Parameter | Var Red (%) | MSE Red (%) | Convergence (PMM2) |",
    "|---|---|---|---|---|---|---|"
)

# Process each file
for (file in files) {
    # Parse filename: Scenario_Innovation_nSize.rds
    # Example: SAR_1_1_gamma_n100.rds
    basename <- tools::file_path_sans_ext(basename(file))
    parts <- strsplit(basename, "_")[[1]]

    n_part <- parts[length(parts)]
    n <- as.integer(sub("n", "", n_part))
    innov <- parts[length(parts) - 1]
    scenario_name <- paste(parts[1:(length(parts) - 2)], collapse = "_")

    results <- readRDS(file)

    # Extract coefficients
    # Structure: list of lists (rep_id, pmm2, css)

    # Determine true parameters based on scenario name (hardcoded for now based on runner)
    true_params <- switch(scenario_name,
        "SAR_1_1" = c(0.5, 0.6), # AR, SAR
        "SMA_1" = c(0.6), # SMA
        "SARMA_1_0_1_1" = c(0.5, 0.6, 0.4), # AR, SAR, SMA
        "SARIMA_1_1_0_0_1_1" = c(0.4, 0.6) # AR, SMA
    )

    param_names <- switch(scenario_name,
        "SAR_1_1" = c("AR", "SAR"),
        "SMA_1" = c("SMA"),
        "SARMA_1_0_1_1" = c("AR", "SAR", "SMA"),
        "SARIMA_1_1_0_0_1_1" = c("AR", "SMA")
    )

    n_params <- length(true_params)

    # Collect estimates
    pmm2_ests <- matrix(NA, nrow = length(results), ncol = n_params)
    css_ests <- matrix(NA, nrow = length(results), ncol = n_params)
    converged_pmm2 <- logical(length(results))

    for (i in seq_along(results)) {
        res <- results[[i]]
        if (!is.null(res$pmm2$coef) && length(res$pmm2$coef) == n_params) {
            pmm2_ests[i, ] <- res$pmm2$coef
            converged_pmm2[i] <- res$pmm2$converged
        }
        if (!is.null(res$css$coef) && length(res$css$coef) == n_params) {
            css_ests[i, ] <- res$css$coef
        }
    }

    # Compute metrics for each parameter
    for (j in 1:n_params) {
        m_pmm2 <- compute_metrics(pmm2_ests[, j], true_params[j])
        m_css <- compute_metrics(css_ests[, j], true_params[j])

        var_red <- (1 - m_pmm2$var / m_css$var) * 100
        mse_red <- (1 - m_pmm2$mse / m_css$mse) * 100
        conv_rate <- mean(converged_pmm2, na.rm = TRUE) * 100

        report_lines <- c(report_lines, sprintf(
            "| %s | %s | %d | %s | %.2f | %.2f | %.1f%% |",
            scenario_name, innov, n, param_names[j], var_red, mse_red, conv_rate
        ))
    }
}

# Write report
writeLines(report_lines, OUTPUT_FILE)
cat("Report generated: ", OUTPUT_FILE, "\n")
