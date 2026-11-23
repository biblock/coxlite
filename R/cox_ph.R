#' Fit a Cox proportional hazards model
#'
#' Fits the Cox model
#' \deqn{h(t \mid x) = h_0(t) \exp(x^\top \beta)}
#' to right-censored data using a formula, like [survival::coxph()]. The
#' actual fitting is done by [cox_fit()].
#'
#' @param formula Formula with `Surv(time, status)` on the left. Use
#'   `strata(z)` on the right to give each level of `z` its own baseline
#'   hazard.
#' @param data Data frame with the variables in `formula`.
#' @param ties `"efron"` (default) or `"breslow"`.
#' @param engine `"cpp"` (default) or `"R"`; see [cox_fit()].
#' @param init Starting values for the coefficients.
#' @param control A list from [cox_control()].
#' @param na.action What to do with missing values (default
#'   [stats::na.omit()]).
#'
#' @return A `"coxlite"` object (see [cox_fit()]). `print()`, `summary()`,
#'   `coef()`, `vcov()`, `confint()`, `logLik()`, `AIC()` and `predict()` all
#'   work on it.
#'
#' @seealso [cox_fit()], [predict.coxlite()], [cox_basehaz()]
#' @export
#' @examples
#' # example data included in the package
#' fit <- cox_ph(Surv(time, status) ~ trt + age + biomarker, data = sim_trial)
#' fit
#' summary(fit)
#'
#' # stratify by centre and use Breslow ties
#' fit_s <- cox_ph(Surv(time, status) ~ trt + age + strata(center),
#'                 data = sim_trial, ties = "breslow")
#' coef(fit_s)
#'
#' # compare with survival
#' fit_ref <- survival::coxph(Surv(time, status) ~ trt + age + biomarker,
#'                            data = sim_trial)
#' all.equal(coef(fit), coef(fit_ref))
cox_ph <- function(formula, data, ties = c("efron", "breslow"),
                   engine = c("cpp", "R"), init = NULL,
                   control = cox_control(), na.action = stats::na.omit) {
  call <- match.call()
  ties <- match.arg(ties)
  engine <- match.arg(engine)
  if (missing(data)) data <- environment(formula)

  terms_all <- stats::terms(formula, specials = "strata", data = data)
  mf <- stats::model.frame(terms_all, data = data, na.action = na.action)
  y <- stats::model.response(mf)
  if (!inherits(y, "Surv")) {
    stop("The response must be a survival object created by Surv().",
         call. = FALSE)
  }
  if (attr(y, "type") != "right") {
    stop("Only right-censored data, Surv(time, status), are supported.",
         call. = FALSE)
  }

  # pull strata() terms out of the formula
  strata_vars <- attr(terms_all, "specials")$strata
  terms_x <- terms_all
  strata <- NULL
  if (length(strata_vars) > 0) {
    factors <- attr(terms_all, "factors")
    strata_terms <- which(colSums(factors[strata_vars, , drop = FALSE]) > 0)
    if (any(colSums(factors[-strata_vars, strata_terms, drop = FALSE]) > 0)) {
      stop("Interactions with strata() terms are not supported.", call. = FALSE)
    }
    if (length(strata_terms) == ncol(factors)) {
      stop("The model needs at least one covariate besides strata().",
           call. = FALSE)
    }
    terms_x <- stats::drop.terms(terms_all, strata_terms, keep.response = TRUE)
    if (length(strata_vars) == 1L) {
      strata <- mf[[strata_vars]]
    } else {
      strata <- interaction(mf[strata_vars], drop = TRUE, sep = ", ")
    }
  }

  # keep the intercept so factors get the usual contrasts, then drop it
  x <- stats::model.matrix(terms_x, mf)
  contrasts <- attr(x, "contrasts")
  x <- x[, colnames(x) != "(Intercept)", drop = FALSE]
  if (ncol(x) == 0L) {
    stop("The model needs at least one covariate.", call. = FALSE)
  }

  fit <- cox_fit(x, time = y[, "time"], status = y[, "status"],
                 strata = strata, ties = ties, engine = engine,
                 init = init, control = control)

  fit$call <- call
  fit$terms <- terms_x
  fit$xlevels <- stats::.getXlevels(terms_x, mf)
  fit$contrasts <- contrasts
  fit$na.action <- attr(mf, "na.action")
  return(fit)
}

#' @importFrom survival Surv
#' @export
survival::Surv

#' @importFrom survival strata
#' @export
survival::strata
