#' Methods for fitted coxlite models
#'
#' Print, summary and extractor methods for fits from [cox_ph()] and
#' [cox_fit()]. `confint()`, `AIC()` and `BIC()` also work through the
#' default methods.
#'
#' @param object,x A `"coxlite"` fit (or its summary).
#' @param conf.int Confidence level for the hazard ratio intervals.
#' @param digits Number of digits to print.
#' @param ... Not used.
#'
#' @return `vcov()` gives the variance matrix, `logLik()` the log partial
#'   likelihood, `nobs()` the number of events (like \pkg{survival}), and
#'   `summary()` a `"summary.coxlite"` object.
#'
#' @name coxlite-methods
#' @examples
#' fit <- cox_ph(Surv(time, status) ~ trt + age, data = sim_trial)
#' coef(fit)
#' vcov(fit)
#' exp(confint(fit))   # CIs for the hazard ratios
#' logLik(fit)
#' AIC(fit)
#' summary(fit, conf.int = 0.9)
NULL

#' @rdname coxlite-methods
#' @export
vcov.coxlite <- function(object, ...) {
  return(object$var)
}

#' @rdname coxlite-methods
#' @export
logLik.coxlite <- function(object, ...) {
  ll <- object$loglik[2L]
  attr(ll, "df") <- length(object$coefficients)
  attr(ll, "nobs") <- object$nevent
  class(ll) <- "logLik"
  return(ll)
}

#' @rdname coxlite-methods
#' @export
nobs.coxlite <- function(object, ...) {
  return(object$nevent)
}

#' @rdname coxlite-methods
#' @export
print.coxlite <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  print_call(x$call)
  print(coef_table(x), digits = digits)
  lr <- lr_test(x)
  cat("\nLikelihood ratio test = ", format(round(lr$stat, 2)),
      " on ", lr$df, " df, p = ", format.pval(lr$p, digits = digits),
      "\n", sep = "")
  cat("n = ", x$n, ", number of events = ", x$nevent, "\n", sep = "")
  invisible(x)
}

#' @rdname coxlite-methods
#' @export
summary.coxlite <- function(object, conf.int = 0.95, ...) {
  beta <- object$coefficients
  se <- sqrt(diag(object$var))
  z <- stats::qnorm((1 + conf.int) / 2)

  # hazard ratios + CIs
  ci <- cbind(exp(beta), exp(-beta), exp(beta - z * se), exp(beta + z * se))
  pct <- format(100 * conf.int, trim = TRUE)
  dimnames(ci) <- list(names(beta), c("exp(coef)", "exp(-coef)",
                                      paste0("lower .", pct),
                                      paste0("upper .", pct)))

  # LR, Wald and score tests (same as coxph)
  df <- length(beta)
  wald_p <- stats::pchisq(object$wald.test, df, lower.tail = FALSE)
  score_p <- stats::pchisq(object$score, df, lower.tail = FALSE)
  tests <- rbind(
    "Likelihood ratio test" = unlist(lr_test(object)),
    "Wald test" = c(object$wald.test, df, wald_p),
    "Score (logrank) test" = c(object$score, df, score_p)
  )
  colnames(tests) <- c("statistic", "df", "p")

  out <- list(call = object$call,
              n = object$n,
              nevent = object$nevent,
              coefficients = coef_table(object),
              conf.int = ci,
              tests = tests,
              ties = object$ties,
              stratified = !is.null(object$strata))
  class(out) <- "summary.coxlite"
  return(out)
}

#' @rdname coxlite-methods
#' @export
print.summary.coxlite <- function(x,
                                  digits = max(3L, getOption("digits") - 3L),
                                  ...) {
  print_call(x$call)
  cat("  n = ", x$n, ", number of events = ", x$nevent,
      "  (ties: ", x$ties, if (x$stratified) ", stratified", ")\n\n",
      sep = "")
  stats::printCoefmat(x$coefficients, digits = digits, P.values = TRUE,
                      has.Pvalue = TRUE, signif.legend = TRUE)
  cat("\n")
  print(x$conf.int, digits = digits)
  cat("\n")
  for (i in seq_len(nrow(x$tests))) {
    cat(formatC(rownames(x$tests)[i], width = -22), "= ",
        format(round(x$tests[i, "statistic"], 2)), " on ",
        x$tests[i, "df"], " df,   p = ",
        format.pval(x$tests[i, "p"], digits = digits), "\n", sep = "")
  }
  invisible(x)
}

print_call <- function(call) {
  if (!is.null(call)) {
    cat("Call:\n")
    print(call)
    cat("\n")
  }
}

coef_table <- function(object) {
  beta <- object$coefficients
  se <- sqrt(diag(object$var))
  z <- beta / se
  return(cbind(coef = beta, "exp(coef)" = exp(beta), "se(coef)" = se, z = z,
               "Pr(>|z|)" = 2 * stats::pnorm(-abs(z))))
}

# LR test vs. the starting values
lr_test <- function(object) {
  stat <- 2 * (object$loglik[2L] - object$loglik[1L])
  df <- length(object$coefficients)
  return(list(statistic = stat, df = df,
              p = stats::pchisq(stat, df, lower.tail = FALSE)))
}
