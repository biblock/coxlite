test_that("sim_cox output looks right", {
  set.seed(1)
  d <- sim_cox(200, beta = c(1, -1), n_strata = 2, round_to = 1)
  expect_named(d, c("time", "status", "x1", "x2", "stratum"))
  expect_true(all(d$status %in% 0:1))
  expect_true(anyDuplicated(d$time) > 0)
})
