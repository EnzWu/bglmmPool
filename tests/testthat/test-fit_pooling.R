test_that("fit_pooling runs on minimal data", {

  skip_if_not_installed("rstan")
  skip_on_cran()

  fit <- fit.pooling(data = stan_data, chains = 1, iter = 200, warmup = 100, check = FALSE)
  expect_true(methods::is(fit, "stanfit"))
  s1 <- rstan::extract(fit, pars = c("RD_overall_all"), permuted = TRUE)
  s1$RD_overall_all
})

test_that("Stan file is shipped with the package", {
  f <- system.file("stan", "Pooling.stan", package = "bglmmPool")
  expect_true(nzchar(f))
  expect_true(file.exists(f))
})

