#' Control parameters for the Newton-Raphson algorithm
#'
#' @param eps Stop when the relative change in the log partial likelihood is
#'   below this (same rule as [survival::coxph.control()]).
#' @param iter.max Max number of Newton-Raphson iterations.
#' @param max.halving Max number of times a step can be halved in one
#'   iteration.
#'
#' @return A list of control parameters.
#' @export
#' @examples
#' cox_control(eps = 1e-10, iter.max = 50)
cox_control <- function(eps = 1e-9, iter.max = 20, max.halving = 20) {
  if (!is.numeric(eps) || length(eps) != 1L || eps <= 0) {
    stop("'eps' must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(iter.max) || length(iter.max) != 1L || iter.max < 1) {
    stop("'iter.max' must be a single number >= 1.", call. = FALSE)
  }
  if (!is.numeric(max.halving) || length(max.halving) != 1L ||
      max.halving < 0) {
    stop("'max.halving' must be a single number >= 0.", call. = FALSE)
  }
  return(list(eps = eps, iter.max = as.integer(iter.max),
              max.halving = as.integer(max.halving)))
}

#' Fit a Cox model from a covariate matrix
#'
#' Fits the model from a numeric covariate matrix instead of a formula.
#' This is what [cox_ph()] calls, and it's useful in simulations.
#'
#' @details
#' Covariates are centred at their means before fitting (except columns that
#' only take -1, 0 and 1, like dummy variables, same as [survival::coxph()]).
#' This doesn't change \eqn{\hat\beta}.
#'
#' Each Newton-Raphson iteration computes the log partial likelihood, the
#' score \eqn{U(\beta)} and the information \eqn{I(\beta)}, then updates
#' \deqn{\beta \leftarrow \beta + I(\beta)^{-1} U(\beta).}
#' If the log-likelihood goes down, the step is halved. The variance estimate
#' is \eqn{I(\hat\beta)^{-1}}. With `strata`, risk sets are formed within
#' each stratum.
#'
#' @param x Numeric covariate matrix, one row per subject. Needs full column
#'   rank.
#' @param time Numeric vector of follow-up times.
#' @param status 1 for an event, 0 for censored (logical is fine too).
#' @param strata Optional vector of strata.
#' @param ties `"efron"` (default) or `"breslow"`.
#' @param engine `"cpp"` (default) or `"R"`. They give the same answer, the
#'   C++ one is just faster.
#' @param init Starting values. Defaults to all zeros.
#' @param control A list from [cox_control()].
#'
#' @return A `"coxlite"` object, which is a list with
#' \describe{
#'   \item{coefficients}{Estimated regression coefficients.}
#'   \item{var}{Variance matrix of the coefficients.}
#'   \item{loglik}{Log partial likelihood at `init` and at the estimate.}
#'   \item{score}{Score test statistic.}
#'   \item{wald.test}{Wald test statistic.}
#'   \item{iter, converged}{Number of iterations, and whether it converged.}
#'   \item{linear.predictors}{Centred linear predictors, in the original
#'     row order.}
#'   \item{means}{Means used for centring.}
#'   \item{n, nevent}{Number of observations and events.}
#'   \item{time, status, strata}{The data, kept for [cox_basehaz()].}
#'   \item{ties, engine}{Which methods were used.}
#' }
#'
#' @seealso [cox_ph()] for the formula interface.
#' @export
#' @examples
#' set.seed(1)
#' d <- sim_cox(200, beta = c(0.5, -1))
#' fit <- cox_fit(as.matrix(d[, c("x1", "x2")]), d$time, d$status)
#' fit
#'
#' # same fit with the R engine
#' fit_r <- cox_fit(as.matrix(d[, c("x1", "x2")]), d$time, d$status,
#'                  engine = "R")
#' all.equal(coef(fit), coef(fit_r))
cox_fit <- function(x, time, status, strata = NULL,
                    ties = c("efron", "breslow"),
                    engine = c("cpp", "R"),
                    init = NULL, control = cox_control()) {
  ties <- match.arg(ties)
  engine <- match.arg(engine)

  # check inputs
  x <- as.matrix(x)
  if (!is.numeric(x)) stop("'x' must be numeric.", call. = FALSE)
  storage.mode(x) <- "double"
  n <- nrow(x)
  p <- ncol(x)
  if (p < 1L) stop("The model needs at least one covariate.", call. = FALSE)
  if (is.null(colnames(x))) colnames(x) <- paste0("x", seq_len(p))
  if (length(time) != n || length(status) != n) {
    stop("'time' and 'status' must have one entry per row of 'x'.",
         call. = FALSE)
  }
  if (!is.numeric(time) || any(!is.finite(time)) || any(time < 0)) {
    stop("'time' must contain finite, non-negative numbers.", call. = FALSE)
  }
  if (is.logical(status)) status <- as.integer(status)
  if (!all(status %in% c(0, 1))) {
    stop("'status' must be coded 0 (censored) / 1 (event).", call. = FALSE)
  }
  status <- as.integer(status)
  if (sum(status) == 0L) stop("There are no events in the data.", call. = FALSE)
  if (any(!is.finite(x))) {
    stop("'x' contains missing or infinite values.", call. = FALSE)
  }

  if (is.null(strata)) {
    strata_f <- NULL
    strata_id <- rep(1L, n)
  } else {
    if (length(strata) != n) {
      stop("'strata' must have one entry per row of 'x'.", call. = FALSE)
    }
    strata_f <- droplevels(as.factor(strata))
    strata_id <- as.integer(strata_f)
  }

  beta <- if (is.null(init)) rep(0, p) else as.numeric(init)
  if (length(beta) != p) stop("'init' must have length ", p, ".", call. = FALSE)

  # centre covariates (but not 0/1 dummies, same as coxph)
  means <- colMeans(x)
  means[apply(x, 2L, function(z) all(z %in% c(-1, 0, 1)))] <- 0
  xc <- sweep(x, 2L, means)

  qx <- qr(xc)
  if (qx$rank < p) {
    aliased <- colnames(x)[qx$pivot[(qx$rank + 1L):p]]
    stop("The design matrix is rank deficient; remove the collinear ",
         "column(s): ", paste(aliased, collapse = ", "), call. = FALSE)
  }

  # engines need rows sorted by stratum, then decreasing time
  ord <- order(strata_id, -time)
  xs <- xc[ord, , drop = FALSE]
  ts <- time[ord]
  ds <- status[ord]
  ss <- strata_id[ord]
  efron <- ties == "efron"
  if (engine == "cpp") {
    derivs <- function(b) cox_derivs_cpp(ts, ds, xs, ss, b, efron)
  } else {
    derivs <- function(b) cox_derivs_r(ts, ds, xs, ss, b, efron)
  }

  nr <- newton_raphson(derivs, beta, control)
  if (!nr$converged) {
    warning("Newton-Raphson did not converge in ", control$iter.max,
            " iterations; consider increasing 'iter.max'.", call. = FALSE)
  }

  coef <- stats::setNames(nr$beta, colnames(x))
  vcov <- chol_inverse(nr$info)
  dimnames(vcov) <- list(colnames(x), colnames(x))

  fit <- list(
    coefficients = coef,
    var = vcov,
    loglik = c(nr$loglik_init, nr$loglik),
    score = nr$score_test,
    wald.test = sum(coef * (nr$info %*% coef)),
    iter = nr$iter,
    converged = nr$converged,
    linear.predictors = drop(xc %*% coef),
    means = means,
    n = n,
    nevent = sum(status),
    time = time,
    status = status,
    strata = strata_f,
    ties = ties,
    engine = engine
  )
  class(fit) <- "coxlite"
  return(fit)
}

# derivs(beta) returns list(loglik, score, info)
newton_raphson <- function(derivs, beta, control) {
  cur <- derivs(beta)
  if (!is.finite(cur$loglik)) {
    stop("The log partial likelihood is not finite at the initial values.",
         call. = FALSE)
  }
  loglik_init <- cur$loglik
  score_test <- sum(cur$score * chol_solve(cur$info, cur$score))

  converged <- FALSE
  iter <- 0L
  while (iter < control$iter.max) {
    iter <- iter + 1L
    step <- chol_solve(cur$info, cur$score)
    new <- derivs(beta + step)

    # step halving
    halvings <- 0L
    slack <- control$eps * abs(cur$loglik)
    while ((!is.finite(new$loglik) || new$loglik < cur$loglik - slack) &&
           halvings < control$max.halving) {
      halvings <- halvings + 1L
      step <- step / 2
      new <- derivs(beta + step)
    }

    change <- abs(new$loglik - cur$loglik)
    beta <- beta + step
    cur <- new
    if (change <= control$eps * abs(cur$loglik)) {
      converged <- TRUE
      break
    }
  }

  return(list(beta = beta, loglik = cur$loglik, info = cur$info,
              loglik_init = loglik_init, score_test = score_test,
              iter = iter, converged = converged))
}

# cholesky helpers for the information matrix
chol_solve <- function(a, b) {
  r <- safe_chol(a)
  return(backsolve(r, forwardsolve(t(r), b)))
}

chol_inverse <- function(a) {
  return(chol2inv(safe_chol(a)))
}

safe_chol <- function(a) {
  tryCatch(chol(a), error = function(e) {
    stop("The information matrix is not positive definite. The model may ",
         "be over-parameterised, or a coefficient may be infinite ",
         "(e.g. a covariate perfectly separates events).", call. = FALSE)
  })
}
