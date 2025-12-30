test_that("fit_pooling runs on minimal data", {

  skip_if_not_installed("rstan")
  skip_on_cran()

  fit <- fit.pooling(data = stan_data, chains = 1, iter = 200, warmup = 100, check = FALSE)
  expect_true(methods::is(fit, "stanfit"))
})


