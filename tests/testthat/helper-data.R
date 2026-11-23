# Check the main pieces of a coxlite fit against a coxph fit
expect_same_fit <- function(fit, ref) {
  expect_equal(coef(fit), coef(ref))
  expect_equal(vcov(fit), vcov(ref), ignore_attr = TRUE)
  expect_equal(fit$loglik, ref$loglik)
}
