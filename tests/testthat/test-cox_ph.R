test_that("matches coxph with Efron and Breslow ties", {
  f <- Surv(time, status) ~ trt + age + biomarker + sex
  for (ties in c("efron", "breslow")) {
    fit <- cox_ph(f, sim_trial, ties = ties)
    ref <- survival::coxph(f, sim_trial, ties = ties)
    expect_same_fit(fit, ref)
  }
})

test_that("matches coxph with strata", {
  f <- Surv(time, status) ~ trt + age + strata(center)
  expect_same_fit(cox_ph(f, sim_trial), survival::coxph(f, sim_trial))
})

test_that("handles factors and missing values like coxph", {
  lung <- survival::lung
  lung$status <- lung$status - 1
  f <- Surv(time, status) ~ age + sex + factor(ph.ecog)
  fit <- cox_ph(f, lung)
  ref <- survival::coxph(f, lung)
  expect_same_fit(fit, ref)
  expect_equal(fit$n, ref$n)
})

test_that("bad formulas give errors", {
  expect_error(cox_ph(time ~ age, sim_trial), "Surv")
  expect_error(cox_ph(Surv(time, status) ~ strata(center), sim_trial),
               "at least one covariate")
})
