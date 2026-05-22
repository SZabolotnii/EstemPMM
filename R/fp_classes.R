# fp_classes.R - S4 class for PMM-FP (Fractional Polynomial PMM) estimator
#
# Part of the PMM-FP extension. See:
#   /Users/serhiizabolotnii/Science/PMM-FP/Spec_PMM_FP.md
#   /Users/serhiizabolotnii/Science/PMM-FP/paper/decisions.md (D1v2, D2)

#' S4 class for storing PMM-FP (Fractional Polynomial) model results
#'
#' Inherits from `PMM2fit`. Adds slots specific to the fractional-polynomial
#' basis search (Royston-Altman 1994; Zabolotnii 2026, PMM-FP article) and
#' the four numerical safeguards applied to the Newton refit (paper §3.3-3.5).
#'
#' @slot selected_powers numeric vector of fractional powers retained in the
#'   final model (subset of `P_a` or `P_b` depending on `track`)
#' @slot track character either `"pos"` (basis `P_a = {0, 0.5, 1, 2, 3}`,
#'   K = 10) or `"full"` (basis `P_b = {-2, -1, -0.5, 0, 0.5, 1, 2, 3}`, K = 16)
#' @slot all_models_ic data.frame with one row per candidate model, holding the
#'   powers, BIC and AIC values; used for diagnostics
#' @slot criterion character information criterion used for selection
#'   (`"BIC"` or `"AIC"`)
#' @slot shift_eps numeric shift applied to `x` when non-positive values were
#'   present (0 if no shift was needed)
#' @slot nr_diagnostics named numeric vector with Newton-Raphson convergence
#'   diagnostics for the selected model: `k_star` (iterations at convergence),
#'   `max_step` (largest parameter change), `final_grad_norm` (||nabla L||),
#'   `kappa_F` (condition number of last Hessian/correlant), and flags
#'   `flag_step_halving`, `flag_tikhonov`, `flag_reinit`, `flag_fallback`
#'   (each 0 or 1 indicating whether the respective safeguard was activated)
#' @slot nr_diverged logical: `TRUE` when both the primary Newton run and the
#'   reinitialised run failed to converge, so the OLS-FP fallback was used
#'
#' @include pmm2_classes.R
#' @exportClass PMM_FP_fit
setClass("PMM_FP_fit",
         contains = "PMM2fit",
         slots = c(selected_powers = "numeric",
                   track           = "character",
                   all_models_ic   = "data.frame",
                   criterion       = "character",
                   shift_eps       = "numeric",
                   nr_diagnostics  = "numeric",
                   nr_diverged     = "logical"))
