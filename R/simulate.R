#' Simulate survival data from a Cox model
#'
#' Simulates right-censored data with a Weibull baseline hazard
#' \eqn{h_0(t) = \lambda \gamma t^{\gamma - 1}}. Event times come from
#' inverting the survival function, and censoring times are exponential
#' (plus an optional cutoff `max_time`). Use `round_to` to get ties.
#'
#' @param n Number of subjects.
#' @param beta True coefficients (one per covariate).
#' @param x Covariate matrix. If `NULL`, covariates are standard normal.
#' @param shape Weibull shape \eqn{\gamma}. 1 gives exponential times.
#' @param scale Weibull scale \eqn{\lambda}.
#' @param cens_rate Rate for exponential censoring (0 = none).
#' @param max_time End of follow-up.
#' @param n_strata Number of strata. Stratum \eqn{s} gets scale
#'   \eqn{s \lambda}.
#' @param round_to Round times up to a multiple of this (to make ties).
#'
#' @return Data frame with `time`, `status`, `x1`, `x2`, ..., and `stratum`
#'   if `n_strata > 1`.
#' @export
#' @examples
#' set.seed(42)
#' d <- sim_cox(500, beta = c(1, -0.5, 0.25), cens_rate = 0.1, round_to = 1)
#' head(d)
#' mean(d$status)   # proportion with an event
#' cox_ph(Surv(time, status) ~ x1 + x2 + x3, data = d)
sim_cox <- function(n, beta, x = NULL, shape = 1, scale = 0.1,
                    cens_rate = 0.05, max_time = Inf, n_strata = 1,
                    round_to = NULL) {
  p <- length(beta)
  if (length(n) != 1L || n < 1) {
    stop("'n' must be a positive integer.", call. = FALSE)
  }
  if (p < 1L) stop("'beta' must have at least one element.", call. = FALSE)
  if (shape <= 0 || scale <= 0) {
    stop("'shape' and 'scale' must be positive.", call. = FALSE)
  }
  if (is.null(x)) {
    x <- matrix(stats::rnorm(n * p), n, p)
  } else {
    x <- as.matrix(x)
    if (nrow(x) != n || ncol(x) != p) {
      stop("'x' must be an n x length(beta) matrix.", call. = FALSE)
    }
  }
  colnames(x) <- paste0("x", seq_len(p))

  stratum <- sample.int(n_strata, n, replace = TRUE)
  rate <- scale * stratum * exp(drop(x %*% beta))
  event_time <- (-log(stats::runif(n)) / rate)^(1 / shape)
  if (cens_rate > 0) {
    cens_time <- stats::rexp(n, cens_rate)
  } else {
    cens_time <- rep(Inf, n)
  }
  cens_time <- pmin(cens_time, max_time)

  time <- pmin(event_time, cens_time)
  status <- as.integer(event_time <= cens_time)
  if (!is.null(round_to)) time <- ceiling(time / round_to) * round_to

  out <- data.frame(time = time, status = status, x)
  if (n_strata > 1) out$stratum <- factor(stratum)
  return(out)
}
