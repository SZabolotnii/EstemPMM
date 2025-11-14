#!/usr/bin/env Rscript

# ==============================================================================
# Focused Monte Carlo study for complex SARMA/SARIMA models
# Innovations: Gaussian and centered Gamma
# Sample sizes: 100, 200, 500; Replications: 1000
# ==============================================================================

suppressPackageStartupMessages({
  library(EstemPMM)
})

set.seed(20250115)

N_REPLICATIONS <- 1000
SAMPLE_SIZES <- c(100, 200, 500)
SEASONAL_PERIOD <- 12
INNOVATION_TYPES <- c("gaussian", "gamma_centered")
OUTPUT_DIR <- "test_results"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("==============================================================================\n")
cat("Monte Carlo Study: SARMA/SARIMA under non-Gaussian innovations\n")
cat("==============================================================================\n")
cat("Replications:", N_REPLICATIONS, "\n")
cat("Sample sizes:", paste(SAMPLE_SIZES, collapse = ", "), "\n")
cat("Seasonal period:", SEASONAL_PERIOD, "\n")
cat("Innovation types:", paste(INNOVATION_TYPES, collapse = ", "), "\n\n")

# ------------------------------------------------------------------------------
# Helper utilities
# ------------------------------------------------------------------------------

generate_innovations <- function(type, n) {
  switch(type,
         "gaussian" = rnorm(n, mean = 0, sd = 1),
         "gamma_centered" = rgamma(n, shape = 2, scale = 1) - 2,
         stop("Unknown innovation type: ", type))
}

generate_sarma_series <- function(n, ar, sar, sma,
                                  s = 12, innovation_type = "gaussian",
                                  burn_in = 200) {
  innovations <- generate_innovations(innovation_type, n + burn_in)

  model_list <- list()
  if (!is.null(ar)) model_list$ar <- ar
  if (!is.null(sar) || !is.null(sma)) {
    seasonal <- list(period = s)
    if (!is.null(sar)) seasonal$sar <- sar
    if (!is.null(sma)) seasonal$sma <- sma
    model_list$seasonal <- seasonal
  }

  series <- arima.sim(n = n + burn_in, model = model_list, innov = innovations)
  as.numeric(tail(series, n))
}

generate_sarima_series <- function(n, ar, sar, sma,
                                   s = 12, innovation_type = "gaussian",
                                   burn_in = 200) {
  innovations <- generate_innovations(innovation_type, n + burn_in)

  series <- arima.sim(
    n = n + burn_in,
    list(order = c(1, 1, 0), ar = ar,
         seasonal = list(order = c(1, 1, 1), period = s,
                         sar = sar, sma = sma)),
    innov = innovations
  )
  as.numeric(tail(series, n))
}

compute_metrics <- function(estimates, true_value) {
  bias <- mean(estimates - true_value)
  rmse <- sqrt(mean((estimates - true_value)^2))
  variance <- var(estimates)
  mae <- mean(abs(estimates - true_value))
  list(bias = bias, rmse = rmse, variance = variance, mae = mae)
}

safe_fit <- function(fit_func, ...) {
  tryCatch({
    fit <- fit_func(...)
    list(
      coefficients = fit@coefficients,
      convergence = fit@convergence,
      m2 = fit@m2,
      m3 = fit@m3,
      m4 = fit@m4
    )
  }, error = function(e) {
    list(coefficients = NA, convergence = FALSE, m2 = NA, m3 = NA, m4 = NA)
  })
}

# ------------------------------------------------------------------------------
# Scenario definitions
# ------------------------------------------------------------------------------

true_sarma <- list(
  name = "SARMA(1,0)×(1,1)_12",
  params = c(AR1 = 0.5, SAR1 = 0.6, SMA1 = 0.4),
  fit_fun = function(y, method) {
    sarma_pmm2(
      y,
      order = c(1, 1, 0, 1),
      season = list(period = SEASONAL_PERIOD),
      method = method,
      verbose = FALSE
    )
  },
  generator = function(n, innovation_type) {
    generate_sarma_series(
      n,
      ar = 0.5,
      sar = 0.6,
      sma = 0.4,
      s = SEASONAL_PERIOD,
      innovation_type = innovation_type
    )
  },
  param_names = c("AR(1)", "SAR(1)", "SMA(1)")
)

true_sarima <- list(
  name = "SARIMA(1,1,0)×(1,1,1)_12",
  params = c(AR1 = 0.4, SAR1 = 0.5, SMA1 = 0.6),
  fit_fun = function(y, method) {
    sarima_pmm2(
      y,
      order = c(1, 1, 0, 1),
      seasonal = list(order = c(1, 1), period = SEASONAL_PERIOD),
      method = method,
      verbose = FALSE
    )
  },
  generator = function(n, innovation_type) {
    generate_sarima_series(
      n,
      ar = 0.4,
      sar = 0.5,
      sma = 0.6,
      s = SEASONAL_PERIOD,
      innovation_type = innovation_type
    )
  },
  param_names = c("AR(1)", "SAR(1)", "SMA(1)")
)

scenarios <- list(sarma = true_sarma, sarima = true_sarima)

all_results <- list()
summary_rows <- list()

for (innov_type in INNOVATION_TYPES) {
  cat("-----------------------------------------------------------------------------\n")
  cat("Innovation type:", innov_type, "\n")
  cat("-----------------------------------------------------------------------------\n")

  scenario_results <- list()

  for (scenario_name in names(scenarios)) {
    scenario <- scenarios[[scenario_name]]
    cat("\nScenario:", scenario$name, "\n")

    scenario_store <- list()

    for (n in SAMPLE_SIZES) {
      cat("  Sample size:", n, "\n")

      pmm2_coef <- matrix(NA_real_, N_REPLICATIONS, length(scenario$params))
      css_coef <- matrix(NA_real_, N_REPLICATIONS, length(scenario$params))
      pmm2_g <- numeric(N_REPLICATIONS)
      pmm2_converged <- 0
      css_converged <- 0

      pb <- txtProgressBar(min = 0, max = N_REPLICATIONS, style = 3)
      for (i in seq_len(N_REPLICATIONS)) {
        series <- scenario$generator(n, innovation_type = innov_type)

        fit_pmm2 <- safe_fit(scenario$fit_fun, series, method = "pmm2")
        if (length(fit_pmm2$coefficients) == length(scenario$params) &&
            !any(is.na(fit_pmm2$coefficients))) {
          pmm2_coef[i, ] <- fit_pmm2$coefficients
          if (!any(is.na(c(fit_pmm2$m2, fit_pmm2$m3, fit_pmm2$m4)))) {
            pmm2_g[i] <- pmm2_variance_factor(
              fit_pmm2$m2, fit_pmm2$m3, fit_pmm2$m4
            )$g
          } else {
            pmm2_g[i] <- NA
          }
          if (isTRUE(fit_pmm2$convergence)) {
            pmm2_converged <- pmm2_converged + 1
          }
        }

        fit_css <- safe_fit(scenario$fit_fun, series, method = "css")
        if (length(fit_css$coefficients) == length(scenario$params) &&
            !any(is.na(fit_css$coefficients))) {
          css_coef[i, ] <- fit_css$coefficients
          if (isTRUE(fit_css$convergence)) {
            css_converged <- css_converged + 1
          }
        }

        setTxtProgressBar(pb, i)
      }
      close(pb)

      valid_pmm2 <- stats::complete.cases(pmm2_coef)
      valid_css <- stats::complete.cases(css_coef)

      metrics_pmm2 <- lapply(seq_along(scenario$params), function(j) {
        compute_metrics(pmm2_coef[valid_pmm2, j], scenario$params[j])
      })
      metrics_css <- lapply(seq_along(scenario$params), function(j) {
        compute_metrics(css_coef[valid_css, j], scenario$params[j])
      })

      scenario_store[[paste0("n", n)]] <- list(
        pmm2 = list(
          coefficients = metrics_pmm2,
          convergence_rate = pmm2_converged / N_REPLICATIONS,
          mean_g = mean(pmm2_g[valid_pmm2], na.rm = TRUE)
        ),
        css = list(
          coefficients = metrics_css,
          convergence_rate = css_converged / N_REPLICATIONS
        )
      )

      for (j in seq_along(scenario$params)) {
        metric_names <- c("bias", "rmse", "variance", "mae")
        p_vals <- metrics_pmm2[[j]]
        c_vals <- metrics_css[[j]]
        var_reduction <- 1 - p_vals$variance / c_vals$variance

        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          innovation = innov_type,
          scenario = scenario$name,
          sample_size = n,
          parameter = scenario$param_names[j],
          method = "PMM2",
          bias = p_vals$bias,
          rmse = p_vals$rmse,
          variance = p_vals$variance,
          mae = p_vals$mae,
          convergence = scenario_store[[paste0("n", n)]]$pmm2$convergence_rate,
          mean_g = scenario_store[[paste0("n", n)]]$pmm2$mean_g,
          variance_reduction = var_reduction
        )

        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          innovation = innov_type,
          scenario = scenario$name,
          sample_size = n,
          parameter = scenario$param_names[j],
          method = "CSS",
          bias = c_vals$bias,
          rmse = c_vals$rmse,
          variance = c_vals$variance,
          mae = c_vals$mae,
          convergence = scenario_store[[paste0("n", n)]]$css$convergence_rate,
          mean_g = NA,
          variance_reduction = NA
        )
      }

      cat("    PMM2 convergence rate:",
          round(scenario_store[[paste0("n", n)]]$pmm2$convergence_rate * 100, 1),
          "%\n")
      cat("    CSS convergence rate:",
          round(scenario_store[[paste0("n", n)]]$css$convergence_rate * 100, 1),
          "%\n")
    }

    scenario_results[[scenario_name]] <- scenario_store
  }

  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  save_path <- file.path(
    OUTPUT_DIR,
    paste0("sarma_sarima_mc_", innov_type, "_", timestamp, ".rds")
  )
  saveRDS(
    list(
      innovation = innov_type,
      scenarios = scenario_results,
      metadata = list(
        replications = N_REPLICATIONS,
        sample_sizes = SAMPLE_SIZES,
        seasonal_period = SEASONAL_PERIOD,
        timestamp = Sys.time()
      )
    ),
    file = save_path
  )
  cat("\nStored detailed results for", innov_type, "at", save_path, "\n\n")

  all_results[[innov_type]] <- scenario_results
}

summary_df <- do.call(rbind, summary_rows)
summary_path <- file.path(
  OUTPUT_DIR,
  paste0("sarma_sarima_mc_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
)
utils::write.csv(summary_df, summary_path, row.names = FALSE)

cat("==============================================================================\n")
cat("Simulation complete. Summary table saved to:", summary_path, "\n")
cat("==============================================================================\n")

# ------------------------------------------------------------------------------
# Markdown report generation
# ------------------------------------------------------------------------------

format_metric <- function(x, digits = 4) {
  ifelse(is.na(x), "NA", sprintf(paste0("%.", digits, "f"), x))
}

format_percent <- function(x) {
  ifelse(is.na(x), "NA", sprintf("%.1f", 100 * x))
}

pm_summary <- subset(summary_df, method == "PMM2")
md_lines <- c(
  "# SARMA/SARIMA Monte Carlo Report",
  "",
  sprintf("- Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("- Replications per setting: %s", N_REPLICATIONS),
  sprintf("- Sample sizes: %s", paste(SAMPLE_SIZES, collapse = ", ")),
  sprintf("- Innovation types: %s", paste(INNOVATION_TYPES, collapse = ", ")),
  sprintf("- Seasonal period: %s", SEASONAL_PERIOD),
  "",
  "## Summary by scenario and innovation"
)

for (innov in unique(pm_summary$innovation)) {
  md_lines <- c(md_lines, "", sprintf("### Innovation: %s", innov))
  innov_df <- pm_summary[pm_summary$innovation == innov, ]

  for (scenario in unique(innov_df$scenario)) {
    scenario_df <- innov_df[innov_df$scenario == scenario, ]
    scenario_df <- scenario_df[order(scenario_df$sample_size,
                                     scenario_df$parameter), ]

    md_lines <- c(
      md_lines,
      "",
      sprintf("#### %s", scenario),
      "",
      "| n | Parameter | Bias | RMSE | Variance | MAE | Var Red (%) | Convergence | mean g |",
      "|---|-----------|------|------|----------|-----|-------------|-------------|--------|"
    )

    for (row_idx in seq_len(nrow(scenario_df))) {
      row <- scenario_df[row_idx, ]
      md_lines <- c(
        md_lines,
        sprintf("| %d | %s | %s | %s | %s | %s | %s | %s | %s |",
                row$sample_size,
                row$parameter,
                format_metric(row$bias),
                format_metric(row$rmse),
                format_metric(row$variance),
                format_metric(row$mae),
                format_percent(row$variance_reduction),
                format_percent(row$convergence),
                format_metric(row$mean_g))
      )
    }
  }
}

md_lines <- c(
  md_lines,
  "",
  sprintf("Raw summary CSV: `%s`", summary_path),
  sprintf("Detailed RDS artifacts saved in `%s`.", OUTPUT_DIR)
)

report_path <- file.path(
  OUTPUT_DIR,
  paste0("sarma_sarima_mc_report_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".md")
)
writeLines(md_lines, con = report_path)

cat("Markdown report saved to:", report_path, "\n")
cat("==============================================================================\n")
