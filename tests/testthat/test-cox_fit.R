test_that("cox_fit catches bad input", {
  x <- matrix(rnorm(20), 10, 2)
  tm <- 1:10
  st <- rep(c(0, 1), 5)
  expect_error(cox_fit(x, tm, st + 1), "coded 0")
  expect_error(cox_fit(x, tm, rep(0, 10)), "no events")
  expect_error(cox_fit(cbind(x, x[, 1]), tm, st), "rank deficient")
})
