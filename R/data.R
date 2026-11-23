#' Simulated randomised trial with time-to-event outcome
#'
#' Fake two-arm trial with 400 patients at four centres, simulated from a
#' Cox model (code in `data-raw/sim_trial.R`). Times are rounded to the
#' nearest half month so there are lots of ties.
#'
#' True log hazard ratios: -0.5 for treatment, 0.03 per year of age, 0.4 for
#' biomarker and 0.2 for male. The baseline hazard is Weibull and differs by
#' centre.
#'
#' @format A data frame with 400 rows and 7 variables:
#' \describe{
#'   \item{time}{Follow-up time in months.}
#'   \item{status}{1 = death, 0 = censored.}
#'   \item{trt}{`control` or `treatment`.}
#'   \item{age}{Age in years.}
#'   \item{biomarker}{Standardised biomarker.}
#'   \item{sex}{`female` or `male`.}
#'   \item{center}{Centre, `A` to `D`.}
#' }
#' @source Simulated, see `data-raw/sim_trial.R`.
#' @examples
#' str(sim_trial)
#' cox_ph(Surv(time, status) ~ trt + age + biomarker + sex + strata(center),
#'        data = sim_trial)
"sim_trial"
