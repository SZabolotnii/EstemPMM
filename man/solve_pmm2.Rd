% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pmm_utils.R
\name{solve_pmm2}
\alias{solve_pmm2}
\title{Universal solver for PMM2 system of equations}
\usage{
solve_pmm2(
  b_init,
  X,
  y,
  m2,
  m3,
  m4,
  max_iter = 1000,
  tol = 1e-05,
  regularize = TRUE,
  reg_lambda = 1e-08,
  verbose = FALSE
)
}
\arguments{
\item{b_init}{Initial parameter estimates (usually from OLS or MLE)}

\item{X}{Design matrix (including intercept and all predictors)}

\item{y}{Response vector}

\item{m2}{Second central moment of residuals}

\item{m3}{Third central moment of residuals}

\item{m4}{Fourth central moment of residuals}

\item{max_iter}{Maximum number of iterations}

\item{tol}{Convergence tolerance}

\item{regularize}{Whether to add regularization to Jacobian matrix}

\item{reg_lambda}{Regularization parameter}

\item{verbose}{Print progress information}
}
\value{
Vector of PMM2 parameter estimates
}
\description{
Universal solver for PMM2 system of equations
}
\keyword{internal}
