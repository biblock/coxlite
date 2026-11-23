#' Cumulative baseline hazard
#'
#' Breslow estimate of the cumulative baseline hazard \eqn{H_0(t)}:
#' \deqn{\hat H_0(t) = \sum_{t_j \le t} \frac{d_j}{\sum_{i \in R(t_j)}
#'   \exp(x_i^\top \hat\beta)},}
#' where \eqn{d_j} is the number of events at \eqn{t_j} and \eqn{R(t_j)} is
#' the risk set. Efron fits use the Efron version of the increment at tied
#' times, like [survival::basehaz()] does. Stratified fits get one curve per
#' stratum.
#'
#' @param fit A fitted `"coxlite"` object.
#' @param centered `TRUE` (default) gives the hazard at the mean covariate
#'   values, `FALSE` gives it at all covariates = 0.
#'
#' @return Data frame with `hazard` and `time` (and `strata` if the model has
#'   strata), same as [survival::basehaz()].
#' @export
#' @examples
#' fit <- cox_ph(Surv(time, status) ~ trt + age, data = sim_trial)
#' H0 <- cox_basehaz(fit)
#' head(H0)
#'
#' # survival curve for an average subject
#' plot(H0$time, exp(-H0$hazard), type = "s",
#'      xlab = "Months", ylab = "Survival probability")
cox_basehaz <- function(fit, centered = TRUE) {
  if (!inherits(fit, "coxlite")) {
    stop("'fit' must be a coxlite model.", call. = FALSE)
  }
  n <- fit$n
  strata_id <- if (is.null(fit$strata)) rep(1L, n) else as.integer(fit$strata)
  ord <- order(strata_id, fit$time)
  time <- fit$time[ord]
  status <- fit$status[ord]
  st <- strata_id[ord]
  risk <- exp(fit$linear.predictors[ord])

  # blocks = rows with the same stratum and time
  new_block <- c(TRUE, time[-1L] != time[-n] | st[-1L] != st[-n])
  block <- cumsum(new_block)
  first <- which(new_block)
  nd <- as.vector(rowsum(status, block))
  d0 <- as.vector(rowsum(risk * status, block))

  # risk set total = sum of risk from this row to the end of the stratum
  at_risk <- stats::ave(risk, st, FUN = function(r) rev(cumsum(rev(r))))
  s0 <- at_risk[first]

  increment <- numeric(length(first))
  for (j in which(nd > 0)) {
    if (fit$ties == "efron") {
      frac <- (seq_len(nd[j]) - 1) / nd[j]
      increment[j] <- sum(1 / (s0[j] - frac * d0[j]))
    } else {
      increment[j] <- nd[j] / s0[j]
    }
  }

  block_strata <- st[first]
  hazard <- stats::ave(increment, block_strata, FUN = cumsum)
  if (!centered) hazard <- hazard * exp(-sum(fit$means * fit$coefficients))

  out <- data.frame(hazard = hazard, time = time[first])
  if (!is.null(fit$strata)) {
    out$strata <- factor(levels(fit$strata)[block_strata],
                         levels = levels(fit$strata))
  }
  return(out)
}
