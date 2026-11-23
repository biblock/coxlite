test_that("cox_basehaz matches survival::basehaz", {
  f <- Surv(time, status) ~ trt + age + strata(center)
  fit <- cox_ph(f, sim_trial)
  ref <- survival::coxph(f, sim_trial)
  h <- cox_basehaz(fit)
  hr <- survival::basehaz(ref)
  expect_equal(h$hazard, hr$hazard)
  expect_equal(h$time, hr$time)
})
