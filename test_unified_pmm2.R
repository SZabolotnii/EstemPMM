source("R/unified_pmm2.R")

# --- Toy Example: Nonlinear Regression y = a * exp(b * x) + epsilon ---
# True parameters
a_true <- 2.0
b_true <- 0.5
n <- 200
set.seed(123)

x <- seq(0, 5, length.out = n)
y_true <- a_true * exp(b_true * x)

# Generate asymmetric errors (Exponential Power Distribution or similar)
# Let's just use Chi-squared centered to have skewness
errors <- (rchisq(n, df = 3) - 3) * 0.5
y <- y_true + errors

# Define Model Functions
# Model: f(theta, x) = theta[1] * exp(theta[2] * x)
# Residuals: e = y - f(theta, x)

fn_residuals_toy <- function(theta) {
    y - (theta[1] * exp(theta[2] * x))
}

fn_jacobian_toy <- function(theta) {
    # J[i, j] = d(f_i)/d(theta_j)
    # d(f)/da = exp(b*x)
    # d(f)/db = a * x * exp(b*x)

    J <- matrix(0, nrow = n, ncol = 2)
    exp_bx <- exp(theta[2] * x)

    J[, 1] <- exp_bx # df/da
    J[, 2] <- theta[1] * x * exp_bx # df/db

    J
}

# --- Test 1: Classical NLS (Least Squares) ---
cat("\n--- Classical NLS ---\n")
nls_fit <- nls(y ~ a * exp(b * x), start = list(a = 1, b = 0.1))
theta_nls <- coef(nls_fit)
print(theta_nls)

# --- Test 2: One-step PMM2 ---
cat("\n--- One-step PMM2 ---\n")
res_onestep <- pmm2_nonlinear_onestep(theta_nls, fn_residuals_toy, fn_jacobian_toy, verbose = TRUE)
print(res_onestep$coefficients)

# --- Test 3: Iterative PMM2 ---
cat("\n--- Iterative PMM2 ---\n")
# Start from a slightly perturbed point to see convergence
theta_start <- theta_nls + c(0.1, -0.05)
res_iter <- pmm2_nonlinear_iterative(theta_start, fn_residuals_toy, fn_jacobian_toy, verbose = TRUE)
print(res_iter$coefficients)

cat("\nTrue Parameters: a =", a_true, ", b =", b_true, "\n")
