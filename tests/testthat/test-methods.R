test_that("summary and logLik match coxph", {
  f <- Surv(time, status) ~ trt + age + biomarker
  fit <- cox_ph(f, sim_trial)
  ref <- survival::coxph(f, sim_trial)
  expect_equal(summary(fit)$conf.int, summary(ref)$conf.int,
               ignore_attr = TRUE)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)))
  expect_output(print(summary(fit)), "Likelihood ratio test")
})
