test_that("fit_arvind_mle recovers theta from simulated data", {
  set.seed(42)
  true_theta <- 2
  x <- rarvind(500, true_theta)
  fit <- fit_arvind_mle(x)
  expect_false(is.na(fit$theta))
  expect_equal(fit$theta, true_theta, tolerance = 0.5)
  expect_true(is.finite(fit$negloglik))
})

test_that("fit_arvind_mle returns NA for insufficient data", {
  fit <- fit_arvind_mle(c(1, 2))  # only 2 observations
  expect_true(is.na(fit$theta))
})

test_that("fit_arvind_mle returns NA for non-finite data", {
  fit <- fit_arvind_mle(c(1, 2, Inf, 4))
  expect_true(is.na(fit$theta))
})

test_that("arvind_mean_fn returns positive values", {
  for (th in c(0.5, 1, 2, 5)) {
    mu <- arvind_mean_fn(th)
    expect_true(!is.na(mu))
    expect_true(mu > 0)
  }
})

test_that("arvind_var_fn returns positive values", {
  for (th in c(0.5, 1, 2, 5)) {
    v <- arvind_var_fn(th)
    expect_true(!is.na(v))
    expect_true(v > 0)
  }
})

test_that("arvind_mean_fn decreases with theta", {
  means <- sapply(c(0.5, 1, 2, 5), arvind_mean_fn)
  # Mean should generally decrease as theta increases
  expect_true(means[1] > means[4])
})
