test_that("R and C++ engines agree", {
  f <- Surv(time, status) ~ trt + age + biomarker + strata(center)
  for (ties in c("efron", "breslow")) {
    a <- cox_ph(f, sim_trial, ties = ties)
    b <- cox_ph(f, sim_trial, ties = ties, engine = "R")
    expect_equal(a[c("coefficients", "var", "loglik")],
                 b[c("coefficients", "var", "loglik")])
  }
})
