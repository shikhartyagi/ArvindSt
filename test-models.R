test_that("fit_rw1 returns ArvindFit object", {
  dat <- simulate_arvind_data(n = 50, seed = 1)
  m <- fit_rw1(Y ~ X1 + X2 + X3, dat, seed = 42)

  expect_s3_class(m, "ArvindFit")
  expect_equal(m$model_type, "RW1-approx")
  expect_length(m$fitted, m$n)
  expect_length(m$residuals, m$n)
  expect_true(is.numeric(m$theta))
  expect_true(m$theta > 0)
})

test_that("diagnostics_arvind returns a data frame with 21 metrics", {
  dat <- simulate_arvind_data(n = 50, seed = 1)
  m <- fit_rw1(Y ~ X1 + X2 + X3, dat, seed = 42)
  d <- diagnostics_arvind(m)

  expect_s3_class(d, "data.frame")
  expect_equal(nrow(d), 1)
  expect_true("MSE" %in% names(d))
  expect_true("RMSE" %in% names(d))
  expect_true("R2" %in% names(d))
  expect_true("AIC" %in% names(d))
  expect_true("KS_stat" %in% names(d))
  expect_true("Theta" %in% names(d))
})

test_that("cv_arvind returns correct structure", {
  dat <- simulate_arvind_data(n = 50, seed = 1)
  m <- fit_rw1(Y ~ X1 + X2 + X3, dat, seed = 42)
  cv <- cv_arvind(m, k_folds = 3, rolling = FALSE, seed = 42)

  expect_true(is.list(cv))
  expect_length(cv$cv_rmse, 3)
  expect_true(cv$mean_cv_rmse > 0)
})

test_that("summary_arvind works with single model", {
  dat <- simulate_arvind_data(n = 50, seed = 1)
  m <- fit_rw1(Y ~ X1 + X2 + X3, dat, seed = 42)
  comp <- summary_arvind(m, comparison_plots = FALSE)

  expect_s3_class(comp, "data.frame")
  expect_equal(nrow(comp), 1)
})
