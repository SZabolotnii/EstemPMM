pkgname <- "EstemPMM"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "EstemPMM-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('EstemPMM')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("DCOILWTICO")
### * DCOILWTICO

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: DCOILWTICO
### Title: WTI Crude Oil Prices
### Aliases: DCOILWTICO
### Keywords: datasets

### ** Examples

data(DCOILWTICO)
head(DCOILWTICO)
summary(DCOILWTICO$DCOILWTICO)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("DCOILWTICO", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ar_pmm2")
### * ar_pmm2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ar_pmm2
### Title: Fit an AR model using PMM2 (wrapper)
### Aliases: ar_pmm2

### ** Examples

## No test: 
# Fit AR(2) model with default variant
x <- arima.sim(n = 200, list(ar = c(0.7, -0.3)))
fit1 <- ar_pmm2(x, order = 2)
coef(fit1)

# Compare variants
fit2 <- ar_pmm2(x, order = 2, pmm2_variant = "unified_iterative")
fit3 <- ar_pmm2(x, order = 2, pmm2_variant = "linearized")
## End(No test)
  



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ar_pmm2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ar_pmm3")
### * ar_pmm3

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ar_pmm3
### Title: Fit an AR model using PMM3
### Aliases: ar_pmm3

### ** Examples

## No test: 
set.seed(42)
x <- arima.sim(n = 200, list(ar = 0.7),
               innov = runif(200, -sqrt(3), sqrt(3)))
fit <- ar_pmm3(x, order = 1)
coef(fit)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ar_pmm3", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("arima_pmm2")
### * arima_pmm2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: arima_pmm2
### Title: Fit an ARIMA model using PMM2 (wrapper)
### Aliases: arima_pmm2

### ** Examples

## No test: 
# Fit ARIMA(1,1,1) model to non-stationary series
x <- cumsum(arima.sim(n = 200, list(ar = 0.6, ma = -0.4)))
fit1 <- arima_pmm2(x, order = c(1, 1, 1))
coef(fit1)

# ARIMA(2,1,0) - random walk with AR(2) innovations
x2 <- cumsum(arima.sim(n = 250, list(ar = c(0.7, -0.3))))
fit2 <- arima_pmm2(x2, order = c(2, 1, 0), pmm2_variant = "unified_global")

# ARIMA(0,2,2) - double differencing with MA(2)
x3 <- cumsum(cumsum(arima.sim(n = 300, list(ma = c(0.5, 0.3)))))
fit3 <- arima_pmm2(x3, order = c(0, 2, 2), pmm2_variant = "unified_iterative")
## End(No test)
  



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("arima_pmm2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("arima_pmm3")
### * arima_pmm3

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: arima_pmm3
### Title: Fit an ARIMA model using PMM3
### Aliases: arima_pmm3

### ** Examples

## No test: 
set.seed(42)
x <- cumsum(arima.sim(n = 200, list(ar = 0.6),
            innov = runif(200, -sqrt(3), sqrt(3))))
fit <- arima_pmm3(x, order = c(1, 1, 0))
coef(fit)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("arima_pmm3", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("arma_pmm2")
### * arma_pmm2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: arma_pmm2
### Title: Fit an ARMA model using PMM2 (wrapper)
### Aliases: arma_pmm2

### ** Examples

## No test: 
# Fit ARMA(2,1) model
x <- arima.sim(n = 250, list(ar = c(0.7, -0.3), ma = 0.5))
fit1 <- arma_pmm2(x, order = c(2, 1))
coef(fit1)

# Try iterative variant for better accuracy
fit2 <- arma_pmm2(x, order = c(2, 1), pmm2_variant = "unified_iterative")

# Higher-order ARMA
x2 <- arima.sim(n = 300, list(ar = c(0.6, -0.2), ma = c(0.4, 0.3)))
fit3 <- arma_pmm2(x2, order = c(2, 2), pmm2_variant = "unified_global")
## End(No test)
  



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("arma_pmm2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("arma_pmm3")
### * arma_pmm3

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: arma_pmm3
### Title: Fit an ARMA model using PMM3
### Aliases: arma_pmm3

### ** Examples

## No test: 
set.seed(42)
x <- arima.sim(n = 250, list(ar = 0.7, ma = -0.3),
               innov = runif(250, -sqrt(3), sqrt(3)))
fit <- arma_pmm3(x, order = c(1, 1))
coef(fit)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("arma_pmm3", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("auto_mpg")
### * auto_mpg

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: auto_mpg
### Title: Auto MPG Dataset
### Aliases: auto_mpg
### Keywords: datasets

### ** Examples

data(auto_mpg)
# PMM2 example: MPG vs Acceleration (asymmetric residuals)
fit_ols <- lm(mpg ~ acceleration, data = auto_mpg)
pmm_skewness(residuals(fit_ols))  # gamma3 ~ 0.5 -> PMM2
pmm_dispatch(residuals(fit_ols))
fit_pmm2 <- lm_pmm2(mpg ~ acceleration, data = auto_mpg, na.action = na.omit)
coef(fit_pmm2)  # compare with coef(fit_ols)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("auto_mpg", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("compare_sar_methods")
### * compare_sar_methods

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: compare_sar_methods
### Title: Compare SAR model estimation methods
### Aliases: compare_sar_methods

### ** Examples

## No test: 
set.seed(42)
y <- arima.sim(n = 120,
  model = list(order = c(1, 0, 0), ar = 0.7,
    seasonal = list(order = c(1, 0, 0), ar = 0.5, period = 12)))
compare_sar_methods(y, order = c(1, 1), period = 12)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compare_sar_methods", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("compare_with_ols")
### * compare_with_ols

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: compare_with_ols
### Title: Compare PMM2 with OLS
### Aliases: compare_with_ols

### ** Examples

## No test: 
result <- compare_with_ols(mpg ~ wt, data = mtcars)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compare_with_ols", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("create_sar_matrix")
### * create_sar_matrix

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: create_sar_matrix
### Title: Create design matrix for seasonal AR model
### Aliases: create_sar_matrix

### ** Examples

## No test: 
# Simple SAR(1)_12 model
x <- rnorm(120)
X <- create_sar_matrix(x, p = 0, P = 1, s = 12)

# AR(1) + SAR(1)_12 additive model
X <- create_sar_matrix(x, p = 1, P = 1, s = 12)

# AR(1) x SAR(1)_12 multiplicative model
X <- create_sar_matrix(x, p = 1, P = 1, s = 12, multiplicative = TRUE)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("create_sar_matrix", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("create_sarma_matrix")
### * create_sarma_matrix

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: create_sarma_matrix
### Title: Create design matrix for seasonal ARMA model
### Aliases: create_sarma_matrix

### ** Examples

## No test: 
# Simple SARMA(1,0) x (1,0)_12 model (AR+SAR, no MA)
x <- rnorm(120)
residuals <- rnorm(120)
X <- create_sarma_matrix(x, residuals, p = 1, P = 1, q = 0, Q = 0, s = 12)

# Full SARMA(1,1) x (1,1)_12 model
X <- create_sarma_matrix(x, residuals, p = 1, P = 1, q = 1, Q = 1, s = 12)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("create_sarma_matrix", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("djia2002")
### * djia2002

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: djia2002
### Title: Dow Jones Industrial Average Daily Data (July-December 2002)
### Aliases: djia2002
### Keywords: datasets

### ** Examples

data(djia2002)
# AR(1) with PMM2
changes <- na.omit(djia2002$change)
pmm_skewness(changes)  # positive skewness -> PMM2
fit <- ar_pmm2(changes, order = 1)
summary(fit)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("djia2002", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lm_pmm2")
### * lm_pmm2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lm_pmm2
### Title: PMM2: Main function for PMM2 (S=2)
### Aliases: lm_pmm2

### ** Examples

set.seed(123)
n <- 80
x <- rnorm(n)
y <- 2 + 3 * x + rt(n, df = 3)
dat <- data.frame(y = y, x = x)

fit <- lm_pmm2(y ~ x, data = dat)
summary(fit, formula = y ~ x, data = dat)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lm_pmm2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lm_pmm3")
### * lm_pmm3

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lm_pmm3
### Title: PMM3: Fit linear model using Polynomial Maximization Method
###   (S=3)
### Aliases: lm_pmm3

### ** Examples

set.seed(123)
n <- 100
x <- rnorm(n)
y <- 2 + 3 * x + runif(n, -sqrt(3), sqrt(3))
dat <- data.frame(y = y, x = x)

fit <- lm_pmm3(y ~ x, data = dat)
summary(fit)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lm_pmm3", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ma_pmm2")
### * ma_pmm2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ma_pmm2
### Title: Fit an MA model using PMM2 (wrapper)
### Aliases: ma_pmm2

### ** Examples

## No test: 
# Fit MA(1) model with linearized variant (recommended)
x <- arima.sim(n = 200, list(ma = 0.6))
fit1 <- ma_pmm2(x, order = 1, pmm2_variant = "linearized")
coef(fit1)

# Compare with unified_global (best accuracy)
fit2 <- ma_pmm2(x, order = 1, pmm2_variant = "unified_global")

# Higher-order MA
x2 <- arima.sim(n = 300, list(ma = c(0.7, -0.4, 0.2)))
fit3 <- ma_pmm2(x2, order = 3, pmm2_variant = "linearized")
## End(No test)
  



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ma_pmm2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ma_pmm3")
### * ma_pmm3

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ma_pmm3
### Title: Fit an MA model using PMM3
### Aliases: ma_pmm3

### ** Examples

## No test: 
set.seed(42)
x <- arima.sim(n = 200, list(ma = 0.6),
               innov = runif(200, -sqrt(3), sqrt(3)))
fit <- ma_pmm3(x, order = 1)
coef(fit)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ma_pmm3", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pmm_dispatch")
### * pmm_dispatch

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pmm_dispatch
### Title: Automatic PMM method selection
### Aliases: pmm_dispatch

### ** Examples

set.seed(42)
x <- rnorm(200); eps <- runif(200, -1, 1)
y <- 1 + 2 * x + eps
fit_ols <- lm(y ~ x)
pmm_dispatch(residuals(fit_ols))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pmm_dispatch", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("predict-PMM2fit-method")
### * predict-PMM2fit-method

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: predict,PMM2fit-method
### Title: Prediction method for PMM2fit objects
### Aliases: predict,PMM2fit-method

### ** Examples

## No test: 
# Fit model
fit <- lm_pmm2(mpg ~ wt + hp, data = mtcars)

# Predict on new data
newdata <- data.frame(wt = c(2.5, 3.0), hp = c(100, 150))
predictions <- predict(fit, newdata = newdata)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("predict-PMM2fit-method", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sar_pmm2")
### * sar_pmm2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sar_pmm2
### Title: Fit Seasonal AR model using PMM2 method
### Aliases: sar_pmm2

### ** Examples

## No test: 
# Generate synthetic seasonal data
n <- 120
y <- arima.sim(n = n, list(ar = 0.7, seasonal = list(sar = 0.5, period = 12)))

# Fit SAR(1,1)_12 model with PMM2
fit <- sar_pmm2(y, order = c(1, 1), season = list(period = 12))
summary(fit)

# Simple seasonal model (no non-seasonal component)
fit_pure_sar <- sar_pmm2(y, order = c(0, 1), season = list(period = 12))

# Compare with OLS
fit_ols <- sar_pmm2(y, order = c(1, 1), season = list(period = 12), method = "ols")
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sar_pmm2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sarima_pmm2")
### * sarima_pmm2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sarima_pmm2
### Title: Fit a Seasonal ARIMA model using PMM2 method
### Aliases: sarima_pmm2

### ** Examples

## No test: 
set.seed(123)
n <- 200
y <- arima.sim(n = n,
  model = list(order = c(1, 0, 1), ar = 0.5, ma = 0.3,
    seasonal = list(order = c(1, 0, 0), ar = 0.4, period = 12)))
fit <- sarima_pmm2(y,
  order = c(1, 0, 1, 0),
  seasonal = list(order = c(1, 0), period = 12))
summary(fit)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sarima_pmm2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sarma_pmm2")
### * sarma_pmm2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sarma_pmm2
### Title: Fit a Seasonal ARMA model using PMM2 method
### Aliases: sarma_pmm2

### ** Examples

## No test: 
# Generate synthetic seasonal data with SARMA structure
set.seed(123)
n <- 200
y <- arima.sim(n = n, list(
  ar = 0.5, ma = 0.3,
  seasonal = list(sar = 0.6, sma = 0.4, period = 12)
))

# Fit SARMA(1,1,1,1)_12 model with PMM2
fit <- sarma_pmm2(y, order = c(1, 1, 1, 1), season = list(period = 12))
summary(fit)

# Pure seasonal model (no non-seasonal components)
fit_pure <- sarma_pmm2(y, order = c(0, 1, 0, 1), season = list(period = 12))
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sarma_pmm2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sma_pmm2")
### * sma_pmm2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sma_pmm2
### Title: Fit a Seasonal MA model using PMM2
### Aliases: sma_pmm2

### ** Examples

## No test: 
# Generate synthetic seasonal data
set.seed(123)
n <- 120
s <- 12
theta <- 0.6

# Gamma innovations (asymmetric)
innov <- rgamma(n, shape = 2, scale = 1) - 2
y <- numeric(n)
for (t in 1:n) {
  ma_term <- if (t > s) theta * innov[t - s] else 0
  y[t] <- innov[t] + ma_term
}

# Fit SMA(1)_12 model with PMM2
fit <- sma_pmm2(y, order = 1, season = list(period = 12))
summary(fit)

# Compare with CSS
fit_css <- sma_pmm2(y, order = 1, season = list(period = 12), method = "css")
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sma_pmm2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
