#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(EstemPMM))
suppressPackageStartupMessages(library(MASS))

get_script_dir <- function() {
  frame_files <- vapply(sys.frames(), function(f) {
    if (!is.null(f$ofile)) f$ofile else NA_character_
  }, character(1), USE.NAMES = FALSE)
  frame_files <- frame_files[!is.na(frame_files)]
  if (length(frame_files) > 0) {
    return(dirname(normalizePath(frame_files[1])))
  }
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- cmd_args[grepl("^--file=", cmd_args)]
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  normalizePath(getwd())
}

script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, ".."))
original_wd <- getwd()
on.exit(setwd(original_wd), add = TRUE)
setwd(project_root)

args <- commandArgs(trailingOnly = TRUE)
reps <- 2000L
seed <- 12345L
standardize_innov <- TRUE
css_method <- "CSS"
allow_m_est <- TRUE
for (arg in args) {
  if (grepl("^--reps=", arg)) {
    reps <- as.integer(sub("^--reps=", "", arg))
  } else if (grepl("^--seed=", arg)) {
    seed <- as.integer(sub("^--seed=", "", arg))
  } else if (grepl("^--standardize-innov=", arg)) {
    value <- tolower(sub("^--standardize-innov=", "", arg))
    standardize_innov <- value %in% c("1", "true", "t", "yes", "y")
  } else if (grepl("^--css-method=", arg)) {
    css_method <- toupper(sub("^--css-method=", "", arg))
  } else if (grepl("^--m-est=", arg)) {
    value <- tolower(sub("^--m-est=", "", arg))
    allow_m_est <- value %in% c("1", "true", "t", "yes", "y")
  }
}
if (!is.finite(reps) || reps < 10L) {
  stop("Argument --reps must be an integer >= 10")
}
if (!is.finite(seed)) {
  stop("Argument --seed must be an integer")
}
if (!css_method %in% c("CSS", "CSS-ML", "ML")) {
  stop("Argument --css-method must be one of CSS, CSS-ML, ML")
}
set.seed(seed)

gaussian_generator <- function(n) stats::rnorm(n)
gamma_raw_generator <- function(n) stats::rgamma(n, shape = 2, scale = 1)
lognormal_sigma <- 0.5513963
lognormal_raw_generator <- function(n) stats::rlnorm(n, meanlog = 0, sdlog = lognormal_sigma)
chisq_raw_generator <- function(n) stats::rchisq(n, df = 3)

lognormal_mean <- exp(lognormal_sigma^2 / 2)
lognormal_sd <- sqrt((exp(lognormal_sigma^2) - 1) * exp(lognormal_sigma^2))

distribution_configs <- list(
  gaussian = list(
    label = "Gaussian",
    raw_generator = gaussian_generator,
    mean_raw = 0,
    sd_raw = 1,
    skewness = 0,
    kurtosis = 0
  ),
  gamma = list(
    label = "Gamma",
    raw_generator = gamma_raw_generator,
    mean_raw = 2,
    sd_raw = sqrt(2),
    skewness = 2 / sqrt(2),
    kurtosis = 6 / 2
  ),
  lognormal = list(
    label = "Lognormal",
    raw_generator = lognormal_raw_generator,
    mean_raw = lognormal_mean,
    sd_raw = lognormal_sd,
    skewness = (exp(lognormal_sigma^2) + 2) * sqrt(exp(lognormal_sigma^2) - 1),
    kurtosis = exp(4 * lognormal_sigma^2) + 2 * exp(3 * lognormal_sigma^2) +
      3 * exp(2 * lognormal_sigma^2) - 6
  ),
  chisq = list(
    label = "Chi-squared",
    raw_generator = chisq_raw_generator,
    mean_raw = 3,
    sd_raw = sqrt(6),
    skewness = sqrt(8 / 3),
    kurtosis = 12 / 3
  )
)

for (dist_name in names(distribution_configs)) {
  cfg <- distribution_configs[[dist_name]]
  cfg$generator <- local({
    raw_fun <- cfg$raw_generator
    mean_raw <- cfg$mean_raw
    sd_raw <- cfg$sd_raw
    function(n) {
      vals <- raw_fun(n)
      if (standardize_innov) {
        if (sd_raw > 0) {
          vals <- (vals - mean_raw) / sd_raw
        } else {
          vals <- vals - mean_raw
        }
      }
      vals
    }
  })
  distribution_configs[[dist_name]] <- cfg
}

monte_carlo_configs <- list(
  list(
    id = "ARIMA(1,1,0)",
    order = c(1, 1, 0),
    ar = c(0.7),
    ma = numeric(0),
    sample_sizes = c(100L, 200L, 500L, 1000L),
    distributions = c("gaussian", "gamma", "lognormal", "chisq"),
    methods = c("CSS", "OLS", "PMM2"),
    true_params = c(phi1 = 0.7),
    baseline_re = "CSS",
    baseline_vr = "OLS"
  ),
  list(
    id = "ARIMA(0,1,1)",
    order = c(0, 1, 1),
    ar = numeric(0),
    ma = c(-0.5),
    sample_sizes = c(500L),
    distributions = c("gaussian", "gamma", "lognormal", "chisq"),
    methods = c("CSS", "PMM2"),
    true_params = c(theta1 = -0.5),
    baseline_re = "CSS",
    baseline_vr = "CSS"
  ),
  list(
    id = "ARIMA(1,1,1)",
    order = c(1, 1, 1),
    ar = c(0.6),
    ma = c(-0.4),
    sample_sizes = c(500L),
    distributions = c("gamma", "lognormal", "chisq"),
    methods = c("CSS", "PMM2"),
    true_params = c(phi1 = 0.6, theta1 = -0.4),
    baseline_re = "CSS",
    baseline_vr = "CSS"
  ),
  list(
    id = "ARIMA(2,1,0)",
    order = c(2, 1, 0),
    ar = c(0.5, 0.3),
    ma = numeric(0),
    sample_sizes = c(500L),
    distributions = c("gamma"),
    methods = c("CSS", "PMM2"),
    true_params = c(phi1 = 0.5, phi2 = 0.3),
    baseline_re = "CSS",
    baseline_vr = "CSS"
  )
)

if (allow_m_est) {
  for (i in seq_along(monte_carlo_configs)) {
    cfg <- monte_carlo_configs[[i]]
    if (cfg$order[1] > 0) {
      cfg$methods <- unique(c(cfg$methods, "M-EST"))
    }
    monte_carlo_configs[[i]] <- cfg
  }
}

calc_skewness <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 3) return(NA_real_)
  m <- mean(x)
  m2 <- mean((x - m)^2)
  m3 <- mean((x - m)^3)
  if (!is.finite(m2) || m2 <= 0) return(NA_real_)
  m3 / (m2^(3 / 2))
}

calc_excess_kurtosis <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 4) return(NA_real_)
  m <- mean(x)
  m2 <- mean((x - m)^2)
  m4 <- mean((x - m)^4)
  if (!is.finite(m2) || m2 <= 0) return(NA_real_)
  m4 / (m2^2) - 3
}

simulate_series <- function(order, ar, ma, n, generator, burnin = 200L) {
  ar_coefs <- if (length(ar)) ar else NULL
  ma_coefs <- if (length(ma)) ma else NULL
  sim <- stats::arima.sim(
    model = list(order = order, ar = ar_coefs, ma = ma_coefs),
    n = n + burnin,
    rand.gen = generator
  )
  as.numeric(tail(sim, n))
}

extract_css <- function(fit, param_names) {
  if (is.null(fit) || !is.list(fit) || is.null(fit$coef)) {
    return(rep(NA_real_, length(param_names)))
  }
  css_names <- names(fit$coef)
  res <- rep(NA_real_, length(param_names))
  names(res) <- param_names
  for (nm in param_names) {
    idx <- as.integer(sub("(phi|theta)", "", nm))
    if (grepl("^phi", nm)) {
      target <- sprintf("ar%d", idx)
    } else {
      target <- sprintf("ma%d", idx)
    }
    pos <- match(target, css_names)
    if (!is.na(pos)) {
      res[nm] <- fit$coef[pos]
    }
  }
  res
}

fit_ols <- function(series, order, param_names) {
  if (order[1] == 0) {
    return(rep(NA_real_, length(param_names)))
  }
  diff_series <- diff(series, differences = order[2])
  embed_mat <- tryCatch(stats::embed(diff_series, order[1] + 1), error = function(e) NULL)
  if (is.null(embed_mat)) return(rep(NA_real_, length(param_names)))
  y <- embed_mat[, 1]
  X <- embed_mat[, -1, drop = FALSE]
  if (nrow(X) <= ncol(X)) return(rep(NA_real_, length(param_names)))
  fit <- tryCatch(stats::lm.fit(x = X, y = y), error = function(e) NULL)
  if (is.null(fit) || any(!is.finite(fit$coefficients))) {
    return(rep(NA_real_, length(param_names)))
  }
  coefs <- as.numeric(fit$coefficients)
  res <- rep(NA_real_, length(param_names))
  res[grepl("^phi", param_names)] <- coefs
  res
}

fit_pmm2 <- function(series, order, param_names) {
  fit <- tryCatch(
    EstemPMM::arima_pmm2(series, order = order, include.mean = FALSE, verbose = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit) || !methods::is(fit, "ARIMAPMM2")) {
    return(list(coefs = rep(NA_real_, length(param_names)), residuals = NA_real_))
  }
  coefs <- as.numeric(fit@coefficients)
  res <- rep(NA_real_, length(param_names))
  names(res) <- param_names
  phi_idx <- seq_len(order[1])
  theta_idx <- seq_len(order[3])
  if (length(phi_idx)) {
    res[seq_along(phi_idx)] <- coefs[phi_idx]
  }
  if (length(theta_idx)) {
    res[length(phi_idx) + seq_along(theta_idx)] <- coefs[length(phi_idx) + seq_along(theta_idx)]
  }
  list(coefs = res, residuals = fit@residuals)
}

fit_m_est <- function(series, order, param_names) {
  if (order[1] == 0) {
    return(rep(NA_real_, length(param_names)))
  }
  if (order[3] > 0) {
    # Robust M-estimation not implemented for MA components
    return(rep(NA_real_, length(param_names)))
  }
  diff_series <- diff(series, differences = order[2])
  embed_mat <- tryCatch(stats::embed(diff_series, order[1] + 1), error = function(e) NULL)
  if (is.null(embed_mat)) return(rep(NA_real_, length(param_names)))
  y <- embed_mat[, 1]
  X <- embed_mat[, -1, drop = FALSE]
  if (nrow(X) <= ncol(X) || any(!is.finite(y)) || any(!is.finite(X))) {
    return(rep(NA_real_, length(param_names)))
  }
  fit <- tryCatch(MASS::rlm(x = X, y = y, intercept = FALSE, psi = MASS::psi.huber,
                            maxit = 100),
                  error = function(e) NULL)
  res <- rep(NA_real_, length(param_names))
  names(res) <- param_names
  if (!is.null(fit) && all(is.finite(fit$coefficients))) {
    res[seq_along(fit$coefficients)] <- as.numeric(fit$coefficients)
  }
  res
}

calc_metrics <- function(values, true_value) {
  finite_vals <- values[is.finite(values)]
  successes <- length(finite_vals)
  if (!successes) {
    return(list(
      bias = NA_real_,
      variance = NA_real_,
      mse = NA_real_,
      rmse = NA_real_,
      successes = 0L
    ))
  }
  diff_vals <- finite_vals - true_value
  bias <- mean(diff_vals)
  variance <- stats::var(finite_vals)
  mse <- mean(diff_vals^2)
  rmse <- sqrt(mse)
  list(
    bias = bias,
    variance = variance,
    mse = mse,
    rmse = rmse,
    successes = successes
  )
}

results_dir <- file.path("results", "monte_carlo")
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

metrics_records <- list()
raw_estimates_store <- list()  # For bootstrap CI computation
residual_store <- lapply(distribution_configs, function(x) list(skew = numeric(0), kurtosis = numeric(0)))

total_tasks <- sum(vapply(monte_carlo_configs, function(cfg) {
  length(cfg$sample_sizes) * length(cfg$distributions)
}, integer(1)))
task_counter <- 0L

for (cfg in monte_carlo_configs) {
  param_names <- names(cfg$true_params)
  for (dist_name in cfg$distributions) {
    dist_cfg <- distribution_configs[[dist_name]]
    for (n in cfg$sample_sizes) {
      task_counter <- task_counter + 1L
      cat(sprintf(
        "[%d/%d] Simulating %s | %s | N = %d | reps = %d\n",
        task_counter, total_tasks, cfg$id, dist_cfg$label, n, reps
      ))
      estimates <- lapply(cfg$methods, function(m) {
        mat <- matrix(NA_real_, nrow = reps, ncol = length(param_names))
        colnames(mat) <- param_names
        mat
      })
      names(estimates) <- cfg$methods
      pb <- utils::txtProgressBar(min = 0, max = reps, style = 3)
      for (r in seq_len(reps)) {
        series <- simulate_series(cfg$order, cfg$ar, cfg$ma, n, dist_cfg$generator)
        if ("CSS" %in% cfg$methods) {
          css_fit <- tryCatch(
            stats::arima(series, order = cfg$order, method = css_method, include.mean = FALSE),
            error = function(e) NULL
          )
          estimates[["CSS"]][r, ] <- extract_css(css_fit, param_names)
        }
        if ("OLS" %in% cfg$methods) {
          estimates[["OLS"]][r, ] <- fit_ols(series, cfg$order, param_names)
        }
        if ("PMM2" %in% cfg$methods) {
          pm_fit <- fit_pmm2(series, cfg$order, param_names)
          estimates[["PMM2"]][r, ] <- pm_fit$coefs
          if (!all(is.na(pm_fit$coefs)) &&
              identical(cfg$id, "ARIMA(1,1,0)") &&
              n == 500L) {
            res_vals <- pm_fit$residuals
            res_vals <- res_vals[is.finite(res_vals)]
            if (length(res_vals) > 10L) {
              s_val <- calc_skewness(res_vals)
              k_val <- calc_excess_kurtosis(res_vals)
              if (is.finite(s_val)) {
                residual_store[[dist_name]]$skew <- c(residual_store[[dist_name]]$skew, s_val)
              }
              if (is.finite(k_val)) {
                residual_store[[dist_name]]$kurtosis <- c(residual_store[[dist_name]]$kurtosis, k_val)
              }
            }
          }
        }
        if ("M-EST" %in% cfg$methods) {
          estimates[["M-EST"]][r, ] <- fit_m_est(series, cfg$order, param_names)
        }
        utils::setTxtProgressBar(pb, r)
      }
      close(pb)

      combo_records <- list()
      for (method in cfg$methods) {
        mat <- estimates[[method]]
        for (param in param_names) {
          stats <- calc_metrics(mat[, param], cfg$true_params[[param]])
          combo_records[[length(combo_records) + 1L]] <- data.frame(
            model = cfg$id,
            distribution = dist_cfg$label,
            sample_size = n,
            method = method,
            parameter = param,
            bias = stats$bias,
            variance = stats$variance,
            mse = stats$mse,
            rmse = stats$rmse,
            successes = stats$successes,
            replicates = reps,
            stringsAsFactors = FALSE
          )
        }
      }
      combo_df <- do.call(rbind, combo_records)
      combo_df$relative_efficiency <- NA_real_
      combo_df$variance_reduction <- NA_real_
      for (param in param_names) {
        baseline_re_val <- combo_df$mse[
          combo_df$parameter == param & combo_df$method == cfg$baseline_re
        ]
        baseline_vr_val <- combo_df$variance[
          combo_df$parameter == param & combo_df$method == cfg$baseline_vr
        ]
        for (i in seq_len(nrow(combo_df))) {
          if (!identical(combo_df$parameter[i], param)) next
          if (!is.na(baseline_re_val) &&
              baseline_re_val > 0 &&
              !is.na(combo_df$mse[i]) &&
              combo_df$mse[i] > 0) {
            combo_df$relative_efficiency[i] <- baseline_re_val / combo_df$mse[i]
          } else {
            combo_df$relative_efficiency[i] <- NA_real_
          }
          if (combo_df$method[i] == cfg$baseline_re) {
            combo_df$relative_efficiency[i] <- 1
          }
          if (combo_df$method[i] == cfg$baseline_vr) {
            combo_df$variance_reduction[i] <- 0
          } else if (!is.na(baseline_vr_val) &&
                     baseline_vr_val > 0 &&
                     !is.na(combo_df$variance[i])) {
            combo_df$variance_reduction[i] <- 100 * (baseline_vr_val - combo_df$variance[i]) / baseline_vr_val
          } else {
            combo_df$variance_reduction[i] <- NA_real_
          }
          if (cfg$baseline_vr == "OLS" && combo_df$method[i] == "CSS") {
            combo_df$variance_reduction[i] <- NA_real_
          }
        }
      }
      metrics_records[[length(metrics_records) + 1L]] <- combo_df

      # Store raw estimates for bootstrap CI computation
      for (method in cfg$methods) {
        mat <- estimates[[method]]
        for (param in param_names) {
          key <- paste(cfg$id, dist_cfg$label, n, method, param, sep = "|")
          raw_estimates_store[[key]] <- mat[, param]
        }
      }
    }
  }
}

metrics_df <- do.call(rbind, metrics_records)
metrics_path <- file.path(results_dir, "monte_carlo_metrics.csv")
utils::write.csv(metrics_df, metrics_path, row.names = FALSE)

# Save raw estimates and configuration for bootstrap CI computation
raw_estimates_path <- file.path(results_dir, "raw_estimates.rds")
saveRDS(list(
  estimates = raw_estimates_store,
  configs = monte_carlo_configs
), raw_estimates_path)
cat(sprintf("Raw estimates saved to: %s\n", raw_estimates_path))

arima110_summary <- metrics_df[
  metrics_df$model == "ARIMA(1,1,0)" & metrics_df$parameter == "phi1",
]
utils::write.csv(
  arima110_summary,
  file.path(results_dir, "arima110_summary.csv"),
  row.names = FALSE
)

re_table <- subset(
  arima110_summary,
  method == "PMM2",
  select = c("distribution", "sample_size", "relative_efficiency")
)
re_wide <- reshape(
  re_table,
  idvar = "distribution",
  timevar = "sample_size",
  direction = "wide"
)
names(re_wide) <- sub("relative_efficiency\\.", "RE_N", names(re_wide))
utils::write.csv(
  re_wide,
  file.path(results_dir, "arima110_re_vs_sample_size.csv"),
  row.names = FALSE
)

arima011_summary <- metrics_df[
  metrics_df$model == "ARIMA(0,1,1)" &
    metrics_df$sample_size == 500L &
    metrics_df$parameter == "theta1",
]
utils::write.csv(
  arima011_summary,
  file.path(results_dir, "arima011_summary.csv"),
  row.names = FALSE
)

arima111_summary <- metrics_df[
  metrics_df$model == "ARIMA(1,1,1)" &
    metrics_df$sample_size == 500L,
]
utils::write.csv(
  arima111_summary,
  file.path(results_dir, "arima111_summary.csv"),
  row.names = FALSE
)

arima210_summary <- metrics_df[
  metrics_df$model == "ARIMA(2,1,0)" &
    metrics_df$sample_size == 500L,
]
utils::write.csv(
  arima210_summary,
  file.path(results_dir, "arima210_summary.csv"),
  row.names = FALSE
)

residual_rows <- list()
for (dist_name in names(residual_store)) {
  dist_cfg <- distribution_configs[[dist_name]]
  skew_vals <- residual_store[[dist_name]]$skew
  kurt_vals <- residual_store[[dist_name]]$kurtosis
  if (!length(skew_vals)) next
  residual_rows[[length(residual_rows) + 1L]] <- data.frame(
    distribution = dist_cfg$label,
    true_skewness = dist_cfg$skewness,
    residual_skewness_mean = mean(skew_vals),
    residual_skewness_sd = stats::sd(skew_vals),
    residual_kurtosis_mean = mean(kurt_vals),
    residual_kurtosis_sd = stats::sd(kurt_vals),
    replicates = length(skew_vals),
    stringsAsFactors = FALSE
  )
}
if (length(residual_rows)) {
  residual_df <- do.call(rbind, residual_rows)
  utils::write.csv(
    residual_df,
    file.path(results_dir, "arima110_residual_cumulants.csv"),
    row.names = FALSE
  )
}

article_reference <- read.csv(text = "
distribution,sample_size,method,target_bias,target_variance,target_mse,target_rmse,target_re,target_vr
Gaussian,100,CSS,-0.0012,0.0180,0.0180,0.1342,1.00,
Gaussian,100,OLS,-0.0008,0.0182,0.0182,0.1349,1.00,0.0
Gaussian,100,PMM2,-0.0010,0.0181,0.0181,0.1346,1.00,0.5
Gaussian,200,CSS,-0.0005,0.0088,0.0088,0.0938,1.00,
Gaussian,200,OLS,-0.0003,0.0089,0.0089,0.0943,0.99,0.0
Gaussian,200,PMM2,-0.0004,0.0088,0.0088,0.0938,1.01,1.1
Gaussian,500,CSS,-0.0001,0.0034,0.0034,0.0583,1.00,
Gaussian,500,OLS,-0.0002,0.0035,0.0035,0.0592,0.97,0.0
Gaussian,500,PMM2,-0.0001,0.0034,0.0034,0.0583,1.03,2.9
Gaussian,1000,CSS,0.0000,0.0017,0.0017,0.0412,1.00,
Gaussian,1000,OLS,0.0000,0.0017,0.0017,0.0412,1.00,0.0
Gaussian,1000,PMM2,0.0000,0.0017,0.0017,0.0412,1.00,0.0
Gamma,100,CSS,-0.0015,0.0195,0.0195,0.1396,1.00,
Gamma,100,OLS,-0.0012,0.0198,0.0198,0.1407,0.98,0.0
Gamma,100,PMM2,-0.0009,0.0142,0.0142,0.1192,1.39,28.3
Gamma,200,CSS,-0.0007,0.0096,0.0096,0.0980,1.00,
Gamma,200,OLS,-0.0005,0.0097,0.0097,0.0985,0.99,0.0
Gamma,200,PMM2,-0.0004,0.0067,0.0067,0.0819,1.45,30.9
Gamma,500,CSS,-0.0002,0.0038,0.0038,0.0616,1.00,
Gamma,500,OLS,-0.0002,0.0038,0.0038,0.0616,1.00,0.0
Gamma,500,PMM2,-0.0001,0.0024,0.0024,0.0490,1.58,36.8
Gamma,1000,CSS,0.0000,0.0019,0.0019,0.0436,1.00,
Gamma,1000,OLS,0.0000,0.0019,0.0019,0.0436,1.00,0.0
Gamma,1000,PMM2,0.0000,0.0012,0.0012,0.0346,1.58,36.8
Lognormal,100,CSS,-0.0018,0.0210,0.0210,0.1449,1.00,
Lognormal,100,OLS,-0.0015,0.0213,0.0213,0.1460,0.99,0.0
Lognormal,100,PMM2,-0.0010,0.0138,0.0138,0.1175,1.54,35.2
Lognormal,200,CSS,-0.0008,0.0103,0.0103,0.1015,1.00,
Lognormal,200,OLS,-0.0006,0.0104,0.0104,0.1020,0.99,0.0
Lognormal,200,PMM2,-0.0004,0.0065,0.0065,0.0806,1.60,37.5
Lognormal,500,CSS,-0.0002,0.0040,0.0040,0.0632,1.00,
Lognormal,500,OLS,-0.0002,0.0041,0.0041,0.0640,0.98,0.0
Lognormal,500,PMM2,-0.0001,0.0024,0.0024,0.0490,1.71,41.5
Lognormal,1000,CSS,0.0000,0.0020,0.0020,0.0447,1.00,
Lognormal,1000,OLS,0.0000,0.0020,0.0020,0.0447,1.00,0.0
Lognormal,1000,PMM2,0.0000,0.0012,0.0012,0.0346,1.67,40.0
Chi-squared,100,CSS,-0.0016,0.0202,0.0202,0.1421,1.00,
Chi-squared,100,OLS,-0.0013,0.0205,0.0205,0.1432,0.99,0.0
Chi-squared,100,PMM2,-0.0008,0.0130,0.0130,0.1140,1.58,36.6
Chi-squared,200,CSS,-0.0007,0.0099,0.0099,0.0995,1.00,
Chi-squared,200,OLS,-0.0005,0.0100,0.0100,0.1000,0.99,0.0
Chi-squared,200,PMM2,-0.0003,0.0058,0.0058,0.0762,1.72,42.0
Chi-squared,500,CSS,-0.0002,0.0039,0.0039,0.0625,1.00,
Chi-squared,500,OLS,-0.0002,0.0040,0.0040,0.0632,0.98,0.0
Chi-squared,500,PMM2,-0.0001,0.0021,0.0021,0.0458,1.90,47.5
Chi-squared,1000,CSS,0.0000,0.0019,0.0019,0.0436,1.00,
Chi-squared,1000,OLS,0.0000,0.0020,0.0020,0.0447,0.95,0.0
Chi-squared,1000,PMM2,0.0000,0.0011,0.0011,0.0332,1.82,45.0
", stringsAsFactors = FALSE)

sim_subset <- subset(
  arima110_summary,
  method %in% c("CSS", "OLS", "PMM2", "M-EST"),
  select = c("distribution", "sample_size", "method", "bias", "variance", "mse", "rmse", "relative_efficiency")
)
names(sim_subset)[names(sim_subset) == "bias"] <- "bias_sim"
names(sim_subset)[names(sim_subset) == "variance"] <- "variance_sim"
names(sim_subset)[names(sim_subset) == "mse"] <- "mse_sim"
names(sim_subset)[names(sim_subset) == "rmse"] <- "rmse_sim"
names(sim_subset)[names(sim_subset) == "relative_efficiency"] <- "re_sim"

article_comparison <- merge(
  sim_subset,
  article_reference,
  by = c("distribution", "sample_size", "method"),
  all = TRUE,
  suffixes = c("_sim", "_article")
)

if (nrow(article_comparison) > 0) {
  article_comparison$bias_diff <- article_comparison$bias_sim - article_comparison$target_bias
  article_comparison$variance_ratio <- article_comparison$variance_sim / article_comparison$target_variance
  article_comparison$mse_ratio <- article_comparison$mse_sim / article_comparison$target_mse
  article_comparison$rmse_ratio <- article_comparison$rmse_sim / article_comparison$target_rmse
}

utils::write.csv(
  article_comparison,
  file.path(results_dir, "article_comparison.csv"),
  row.names = FALSE
)

cat("\nMonte Carlo simulation completed.\n")
cat(sprintf("Metrics saved to: %s\n", metrics_path))
