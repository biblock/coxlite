test_that("predict matches coxph on new data", {
  f <- Surv(time, status) ~ trt + age + biomarker + sex
  fit <- cox_ph(f, sim_trial)
  ref <- survival::coxph(f, sim_trial)
  nd <- sim_trial[1:5, ]
  p <- predict(fit, nd, se.fit = TRUE)
  pr <- predict(ref, nd, se.fit = TRUE, reference = "sample")
  expect_equal(p$fit, pr$fit, ignore_attr = TRUE)
  expect_equal(p$se.fit, pr$se.fit, ignore_attr = TRUE)
})
