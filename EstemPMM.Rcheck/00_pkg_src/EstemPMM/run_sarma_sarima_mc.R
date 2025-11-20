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
METHODS <- c("pmm2", "css", "mle")

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

      coef_store <- lapply(METHODS, function(m) {
        matrix(NA_real_, N_REPLICATIONS, length(scenario$params))
      })
      names(coef_store) <- METHODS

      convergence_counts <- setNames(numeric(length(METHODS)), METHODS)
      pmm2_g <- numeric(N_REPLICATIONS)

      pb <- txtProgressBar(min = 0, max = N_REPLICATIONS, style = 3)
      for (i in seq_len(N_REPLICATIONS)) {
        series <- scenario$generator(n, innovation_type = innov_type)

        for (method in METHODS) {
          fit <- safe_fit(scenario$fit_fun, series, method = method)
          if (length(fit$coefficients) == length(scenario$params) &&
              !any(is.na(fit$coefficients))) {
            coef_store[[method]][i, ] <- fit$coefficients
            if (isTRUE(fit$convergence)) {
              convergence_counts[method] <- convergence_counts[method] + 1
            }
            if (method == "pmm2" &&
                !any(is.na(c(fit$m2, fit$m3, fit$m4)))) {
              pmm2_g[i] <- pmm2_variance_factor(
                fit$m2, fit$m3, fit$m4
              )$g
            }
          }
        }

        setTxtProgressBar(pb, i)
      }
      close(pb)

valid_mask <- lapply(coef_store, stats::complete.cases)
names(valid_mask) <- names(coef_store)

      metrics <- lapply(METHODS, function(method) {
        lapply(seq_along(scenario$params), function(j) {
          compute_metrics(coef_store[[method]][valid_mask[[method]], j],
                          scenario$params[j])
        })
      })
      names(metrics) <- METHODS

      css_variances <- sapply(seq_along(scenario$params), function(j) {
        metrics$css[[j]]$variance
      })

      scenario_store[[paste0("n", n)]] <- list(
        pmm2 = list(
          coefficients = metrics$pmm2,
          convergence_rate = convergence_counts["pmm2"] / N_REPLICATIONS,
          mean_g = mean(pmm2_g[valid_mask$pmm2], na.rm = TRUE)
        ),
        css = list(
          coefficients = metrics$css,
          convergence_rate = convergence_counts["css"] / N_REPLICATIONS
        ),
        mle = list(
          coefficients = metrics$mle,
          convergence_rate = convergence_counts["mle"] / N_REPLICATIONS
        )
      )

      for (method in METHODS) {
        method_pretty <- toupper(method)
        method_pretty <- if (method_pretty == "PMM2") "PMM2" else method_pretty

        for (j in seq_along(scenario$params)) {
          vals <- metrics[[method]][[j]]
          var_reduction <- if (method == "css") {
            NA
          } else {
            1 - vals$variance / css_variances[j]
          }

          summary_rows[[length(summary_rows) + 1]] <- data.frame(
            innovation = innov_type,
            scenario = scenario$name,
            sample_size = n,
            parameter = scenario$param_names[j],
            method = method_pretty,
            bias = vals$bias,
            rmse = vals$rmse,
            variance = vals$variance,
            mae = vals$mae,
            convergence = scenario_store[[paste0("n", n)]][[method]]$convergence_rate,
            mean_g = if (method == "pmm2") scenario_store[[paste0("n", n)]]$pmm2$mean_g else NA,
            variance_reduction = var_reduction
          )
        }
      }

      for (method in METHODS) {
        cat("    ", toupper(method), " convergence rate:",
            round(scenario_store[[paste0("n", n)]][[method]]$convergence_rate * 100, 1),
            "%\n", sep = "")
      }
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

md_data <- summary_df
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

for (innov in unique(md_data$innovation)) {
  md_lines <- c(md_lines, "", sprintf("### Innovation: %s", innov))
  innov_df <- md_data[md_data$innovation == innov, ]

  for (scenario in unique(innov_df$scenario)) {
    scenario_df <- innov_df[innov_df$scenario == scenario, ]
    scenario_df <- scenario_df[order(scenario_df$sample_size,
                                     scenario_df$parameter,
                                     scenario_df$method), ]

    md_lines <- c(
      md_lines,
      "",
      sprintf("#### %s", scenario),
      "",
      "| n | Parameter | Method | Bias | RMSE | Variance | MAE | Var Red (%) | Convergence | mean g |",
      "|---|-----------|--------|------|------|----------|-----|-------------|-------------|--------|"
    )

    for (row_idx in seq_len(nrow(scenario_df))) {
      row <- scenario_df[row_idx, ]
      md_lines <- c(
        md_lines,
        sprintf("| %d | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
                row$sample_size,
                row$parameter,
                row$method,
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
