#' Predictions from a fitted Cox model
#'
#' Linear predictor or relative risk for the original data or new data.
#'
#' @param object A `"coxlite"` fit.
#' @param newdata Data frame for [cox_ph()] fits, or a matrix with the same
#'   columns for [cox_fit()] fits. Leave it out to predict on the fitted data.
#' @param type `"lp"` (linear predictor) or `"risk"` (\eqn{\exp} of it).
#' @param reference `"sample"` (default) centres at the covariate means from
#'   the fit, `"zero"` doesn't. Same as in [survival::predict.coxph()].
#' @param se.fit Return standard errors too? Only works with `newdata`. For
#'   `type = "risk"` this uses the delta method, so it won't match
#'   \pkg{survival}, which does it differently.
#' @param ... Not used.
#'
#' @return A vector of predictions, or a list with `fit` and `se.fit`.
#' @export
#' @examples
#' fit <- cox_ph(Surv(time, status) ~ trt + age, data = sim_trial)
#' head(predict(fit))
#'
#' new_patients <- data.frame(trt = c("control", "treatment"), age = 60)
#' predict(fit, new_patients, type = "risk", se.fit = TRUE)
predict.coxlite <- function(object, newdata, type = c("lp", "risk"),
                            reference = c("sample", "zero"),
                            se.fit = FALSE, ...) {
  type <- match.arg(type)
  reference <- match.arg(reference)
  beta <- object$coefficients
  shift <- sum(object$means * beta)

  if (missing(newdata)) {
    x <- NULL
    lp <- object$linear.predictors
    if (reference == "zero") lp <- lp + shift
  } else {
    x <- new_design(object, newdata)
    lp <- drop(x %*% beta)
    if (reference == "sample") lp <- lp - shift
    names(lp) <- rownames(x)
  }

  fit <- if (type == "lp") lp else exp(lp)
  if (!se.fit) return(fit)

  if (is.null(x)) {
    stop("se.fit = TRUE requires 'newdata' (the fitted design matrix is ",
         "not stored).", call. = FALSE)
  }
  if (reference == "sample") x <- sweep(x, 2L, object$means)
  se <- sqrt(rowSums((x %*% object$var) * x))
  if (type == "risk") se <- se * fit
  return(list(fit = fit, se.fit = se))
}

# design matrix for newdata
new_design <- function(object, newdata) {
  # cox_fit() models don't have terms, newdata should already be a matrix
  if (is.null(object$terms)) {
    x <- as.matrix(newdata)
    if (ncol(x) != length(object$coefficients)) {
      stop("'newdata' must have ", length(object$coefficients), " columns.",
           call. = FALSE)
    }
    return(x)
  }
  tt <- stats::delete.response(object$terms)
  mf <- stats::model.frame(tt, newdata, xlev = object$xlevels,
                           na.action = stats::na.pass)
  x <- stats::model.matrix(tt, mf, contrasts.arg = object$contrasts)
  return(x[, colnames(x) != "(Intercept)", drop = FALSE])
}
