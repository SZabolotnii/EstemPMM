# fp_classes.R - S4 class for PMM-FP (Fractional Polynomial PMM) estimator
#
# Part of the PMM-FP extension. See:
#   /Users/serhiizabolotnii/Science/PMM-FP/Spec_PMM_FP.md
#   /Users/serhiizabolotnii/Science/PMM-FP/paper/decisions.md (D1v2, D2)

#' S4 class for storing PMM-FP (Fractional Polynomial) model results
#'
#' Inherits from `PMM2fit`. Adds slots specific to the fractional-polynomial
#' basis search (Royston-Altman 1994; Zabolotnii 2026, PMM-FP article).
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
#'
#' @include pmm2_classes.R
#' @exportClass PMM_FP_fit
setClass("PMM_FP_fit",
         contains = "PMM2fit",
         slots = c(selected_powers = "numeric",
                   track           = "character",
                   all_models_ic   = "data.frame",
                   criterion       = "character",
                   shift_eps       = "numeric"))
