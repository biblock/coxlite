#' coxlite: A Lightweight Re-Implementation of the Cox Proportional Hazards
#' Model
#'
#' Fits Cox proportional hazards models with its own Newton-Raphson code
#' instead of calling [survival::coxph()]. Start with [cox_ph()], and see
#' `vignette("coxlite")` for an intro and `vignette("benchmarks")` for the
#' comparison with \pkg{survival}.
#'
#' @keywords internal
#' @useDynLib coxlite, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats vcov logLik nobs predict
"_PACKAGE"
