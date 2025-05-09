################################################################
## 1. Необхідні пакети й допоміжні функції
################################################################


library(EstemPMM)
library(dplyr)
library(moments)     # для skewness(), kurtosis()
library(ggplot2)   # опційно, якщо захочете ggplot-візуалізації

# (Припускаємо, що у вас є функції ar_pmm2, ma_pmm2, arma_pmm2, arima_pmm2)

# Генерація AR(p), MA(q), ARMA, ARIMA рядів із гамма-інноваціями
generate_ar <- function(n, ar_coefs, shape_gamma = 2, burn_in = 50) {
  sim_data <- arima.sim(
    model = list(ar = ar_coefs),
    n = n + burn_in,
    rand.gen = function(nn) {
      (rgamma(nn, shape = shape_gamma) - shape_gamma)/sqrt(shape_gamma)
    }
  )
  sim_data[(burn_in+1):(n+burn_in)]
}


generate_ma <- function(n, ma_coefs, shape_gamma = 2, burn_in = 50) {
  sim_data <- arima.sim(
    model = list(ma = ma_coefs),
    n = n + burn_in,
    rand.gen = function(nn) {
      (rgamma(nn, shape = shape_gamma) - shape_gamma)/sqrt(shape_gamma)
    }
  )
  sim_data[(burn_in+1):(n+burn_in)]
}

generate_arma <- function(n, ar_coefs, ma_coefs, shape_gamma = 2, burn_in = 50) {
  sim_data <- arima.sim(
    model = list(ar = ar_coefs, ma = ma_coefs),
    n = n + burn_in,
    rand.gen = function(nn) {
      (rgamma(nn, shape = shape_gamma) - shape_gamma)/sqrt(shape_gamma)
    }
  )
  sim_data[(burn_in+1):(n+burn_in)]
}

generate_arima <- function(n, ar_coefs, d, ma_coefs, shape_gamma = 2, burn_in = 50) {
  # Генеруємо стаціонарний ARMA(p,q)
  sim_arma <- arima.sim(
    model = list(ar = ar_coefs, ma = ma_coefs),
    n = n + burn_in,
    rand.gen = function(nn) {
      (rgamma(nn, shape = shape_gamma) - shape_gamma)/sqrt(shape_gamma)
    }
  )
  sim_arma_clean <- sim_arma[(burn_in+1):(n+burn_in)]
  # Робимо d-разів інтеграцію
  if(d > 0) {
    for(i in seq_len(d)) {
      sim_arma_clean <- cumsum(sim_arma_clean)
    }
  }
  sim_arma_clean
}

################################################################
## 2. Універсальна функція fit_model: підгонка і повернення (коефіцієнти, залишки)
################################################################

fit_model <- function(series, model_type, method,
                      ar_coefs=NULL, ma_coefs=NULL, d=0, include.mean=FALSE) {
  # Повертає list(coefs=..., res=..., c3=..., c4=..., g2=...)

  if(model_type=="AR") {
    p <- length(ar_coefs)
    # Оцінка:
    if(method=="MLE" || method=="ML") {
      fit <- arima(series, order = c(p,0,0), method = "ML", include.mean=include.mean)
      coefs_est <- unname(fit$coef[1:p])
      res <- residuals(fit)
    } else if(method=="OLS") {
      # Проста регресія
      n <- length(series)
      X <- matrix(0, nrow=n-p, ncol=p)
      for(i in seq_len(p)) {
        X[,i] <- series[(p - i + 1):(n - i)]
      }
      y <- series[(p+1):n]
      if(nrow(X) == length(y)) {
        b_ols <- lm.fit(X,y)$coefficients
        coefs_est <- as.numeric(b_ols)
        # Залишки:
        res <- y - X %*% b_ols
      } else {
        coefs_est <- rep(NA, p)
        res <- rep(NA, length(series))
      }
    } else if(method=="YW") {
      f_yw <- ar.yw(series, aic=FALSE, order.max=p)
      coefs_est <- f_yw$ar
      # Залишки:
      # ar.yw не зберігає напряму *коректні* res, але "f_yw$resid" є
      # Використаємо:
      res <- f_yw$resid
    } else if(method=="PMM2") {
      f_pmm2 <- ar_pmm2(series, order=p, method="pmm2", include.mean=include.mean)
      coefs_est <- f_pmm2@coefficients
      res <- f_pmm2@residuals
    } else {
      stop("Unsupported AR method: ", method)
    }

  } else if(model_type=="MA") {
    q <- length(ma_coefs)
    if(method %in% c("CSS","ML")) {
      fit <- arima(series, order=c(0,0,q), method=method, include.mean=include.mean)
      coefs_est <- unname(fit$coef[1:q])
      res <- residuals(fit)
    } else if(method=="PMM2") {
      f_pmm2 <- ma_pmm2(series, order=q, method="pmm2", include.mean=include.mean)
      coefs_est <- f_pmm2@coefficients
      res <- f_pmm2@residuals
    } else {
      stop("Unsupported MA method: ", method)
    }

  } else if(model_type=="ARMA") {
    p <- length(ar_coefs)
    q <- length(ma_coefs)
    if(method %in% c("CSS","ML")) {
      fit <- arima(series, order=c(p,0,q), method=method, include.mean=include.mean)
      coefs_est <- unname(fit$coef[1:(p+q)])
      res <- residuals(fit)
    } else if(method=="PMM2") {
      f_pmm2 <- arma_pmm2(series, order=c(p,q), method="pmm2", include.mean=include.mean)
      coefs_est <- f_pmm2@coefficients
      res <- f_pmm2@residuals
    } else {
      stop("Unsupported ARMA method: ", method)
    }

  } else if(model_type=="ARIMA") {
    p <- length(ar_coefs)
    q <- length(ma_coefs)
    if(method %in% c("CSS","ML")) {
      fit <- arima(series, order=c(p,d,q), method=method, include.mean=include.mean)
      coefs_est <- unname(fit$coef[1:(p+q)])
      res <- residuals(fit)
    } else if(method=="PMM2") {
      f_pmm2 <- arima_pmm2(series, order=c(p,d,q), method="pmm2", include.mean=include.mean)
      coefs_est <- f_pmm2@coefficients
      res <- f_pmm2@residuals
    } else {
      stop("Unsupported ARIMA method: ", method)
    }

  } else {
    stop("Unknown model_type in fit_model")
  }

  # Обчислення c3, c4 і g2:
  # skewness() і kurtosis() з пакета moments
  # kurtosis(...) дає повну 4-ту моментну характеристику => ексцес = kurtosis(...) - 3
  c3 <- if(all(is.finite(res))) skewness(res, na.rm=TRUE) else NA
  # c4 = excess kurtosis
  c4 <- if(all(is.finite(res))) kurtosis(res, na.rm=TRUE) - 3 else NA
  g2 <- if(!is.na(c3) && !is.na(c4)) 1 - c3^2 / (2 + c4) else NA

  list(coefs = coefs_est, res = res, c3 = c3, c4 = c4, g2 = g2)
}


################################################################
## 3. Основний цикл Монте-Карло
################################################################

W <- 1000 # ількість експериментів
N <- 100 # Довжина часового ряду

monte_carlo_comparison <- function(R = W,
                                   model_list = list(
                                     list(type="AR",    ar=c(0.7, -0.3), ma=NULL, d=0, n=N),
                                     list(type="MA",    ar=NULL, ma=c(0.6, -0.1), d=0, n=N),
                                     list(type="ARMA",  ar=0.7,  ma=-0.4,  d=0, n=N),
                                     list(type="ARIMA", ar=0.7,  ma=-0.4,  d=1, n=N)
                                   ),
                                   shape_gamma = 2,
                                   methods_ar   = c("MLE","OLS","PMM2","YW"),
                                   methods_ma   = c("CSS","ML","PMM2"),
                                   methods_arma = c("CSS","ML","PMM2"),
                                   methods_arima= c("CSS","ML","PMM2")
){
  results <- data.frame()

  for(mdl in model_list) {
    mtype <- toupper(mdl$type)  # "AR", "MA", "ARMA", "ARIMA"
    ar_params <- mdl$ar
    ma_params <- mdl$ma
    d         <- mdl$d
    n         <- mdl$n

    # Істинні коефіцієнти та назви
    if(mtype=="AR") {
      true_coefs  <- ar_params
      param_names <- paste0("ar", seq_along(ar_params))
      method_set  <- methods_ar
      gen_fun <- function() generate_ar(n, ar_params, shape_gamma = shape_gamma)
    } else if(mtype=="MA") {
      true_coefs  <- ma_params
      param_names <- paste0("ma", seq_along(ma_params))
      method_set  <- methods_ma
      gen_fun <- function() generate_ma(n, ma_params, shape_gamma = shape_gamma)
    } else if(mtype=="ARMA") {
      p <- length(ar_params)
      q <- length(ma_params)
      true_coefs  <- c(ar_params, ma_params)
      param_names <- c(paste0("ar", seq_len(p)), paste0("ma", seq_len(q)))
      method_set  <- methods_arma
      gen_fun <- function() generate_arma(n, ar_params, ma_params, shape_gamma = shape_gamma)
    } else if(mtype=="ARIMA") {
      true_coefs  <- c(ar_params, ma_params)
      param_names <- c("ar1", "ma1")
      method_set  <- methods_arima
      gen_fun <- function() generate_arima(n, ar_params, d, ma_params, shape_gamma=shape_gamma)
    } else {
      stop("Unknown model type in model_list: ", mtype)
    }

    # Цикл по R симуляціях
    for(rep_i in seq_len(R)) {
      series <- gen_fun()

      for(mtd in method_set) {
        # Викликаємо fit_model
        fobj <- fit_model(series, model_type = mtype, method=mtd,
                          ar_coefs=ar_params, ma_coefs=ma_params, d=d,
                          include.mean=FALSE)
        est <- fobj$coefs
        # Розміри можуть не співпадати, але припустимо, що length(est) = length(true_coefs)

        # Зберігаємо в data.frame - один рядок на кожний параметр,
        # але skewness, kurtosis, g2 однакові для всіх параметрів => дублюємо
        for(i in seq_along(est)) {
          tmp <- data.frame(
            rep       = rep_i,
            model     = mtype,
            method    = mtd,
            param     = param_names[i],
            estimate  = est[i],
            true_value= true_coefs[i],
            c3        = fobj$c3,
            c4        = fobj$c4,
            g2        = fobj$g2
          )
          results <- rbind(results, tmp)
        }
      }
    }
  }
  results
}

################################################################
## 4. Запуск, статистики та візуалізація
################################################################

set.seed(123)
mc_results <- monte_carlo_comparison(R = 200)  # Змінюйте R за потреби

# Обчислимо Bias, MSE
stats_summary <- mc_results %>%
  group_by(model, method, param) %>%
  summarise(
    mean_est = mean(estimate, na.rm=TRUE),
    bias     = mean(estimate - true_value, na.rm=TRUE),
    MSE      = mean((estimate - true_value)^2, na.rm=TRUE),
    # Середні c3, c4, g2 (по всіх симуляціях)
    mean_c3  = mean(c3, na.rm=TRUE),
    mean_c4  = mean(c4, na.rm=TRUE),
    mean_g2  = mean(g2, na.rm=TRUE),
    .groups  = "drop"
  )

cat("\n===== Summary of Bias, MSE, and mean(g2) =====\n")
print(stats_summary, n=50)

# 4.1. Порівнюємо з PMM2, але тепер ratio = MSE_PMM2 / MSE
pm2_only <- stats_summary %>%
  filter(method == "PMM2") %>%
  rename(MSE_PMM2 = MSE, bias_PMM2 = bias) %>%
  select(model, param, MSE_PMM2, bias_PMM2)

stats_rel <- stats_summary %>%
  left_join(pm2_only, by=c("model","param")) %>%
  mutate(
    MSE_ratio = MSE_PMM2 / MSE,    # Якщо >1 => у PMM2 гірша MSE
    # (або краще?? Залежить як читаєте)
    MSE_diff  = MSE - MSE_PMM2,
    bias_diff = bias - bias_PMM2
  )

cat("\n===== Relative comparison to PMM2 =====\n")
print(stats_rel, n=50)

# 5. Бокс-плоти оцінок, без «outliers»
unique_modpar <- unique(paste(mc_results$model, mc_results$param, sep="_"))

par(mfrow = c(2,2))
for(mmp in unique_modpar) {
  sub_df <- subset(mc_results, paste(model, param, sep="_") == mmp)
  true_val <- unique(sub_df$true_value)
  main_title <- paste0("Model=", unique(sub_df$model),
                       " Param=", unique(sub_df$param))

  boxplot(estimate ~ method, data=sub_df,
          main = main_title, ylab = "Estimate", xlab = "Method",
          col = "lightblue", border = "darkblue",
          outline=FALSE  # ось головний ключ, що не показує «outliers»
  )
  abline(h = true_val, col="red", lty=2)
}

# Якщо хочете подивитись також на g2 (чи c3, c4) - можна побудувати аналогічні бокс-плоти:
# boxplot(g2 ~ method, data=subset(mc_results, model=="MA"), outline=FALSE)

################################################################
## Кінець скрипта
################################################################
