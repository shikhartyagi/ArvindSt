test_that("rarvind generates n positive values", {
  set.seed(42)
  x <- rarvind(100, theta = 2)
  expect_length(x, 100)
  expect_true(all(x > 0))
})

test_that("rarvind sample mean is close to theoretical mean", {
  set.seed(42)
  theta <- 2
  x <- rarvind(5000, theta)
  mu_th <- arvind_mean_fn(theta)
  expect_equal(mean(x), mu_th, tolerance = 0.1)
})

test_that("rarvind_centred has mean approximately 0", {
  set.seed(42)
  eps <- rarvind_centred(5000, theta = 2)
  expect_length(eps, 5000)
  expect_equal(mean(eps), 0, tolerance = 0.1)
})

test_that("rarvind works for different theta values", {
  set.seed(42)
  for (th in c(0.5, 1, 5, 10)) {
    x <- rarvind(50, th)
    expect_length(x, 50)
    expect_true(all(x > 0))
  }
})
