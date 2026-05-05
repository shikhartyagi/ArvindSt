test_that("darvind handles vector theta edge case", {
  # Single theta applied to vector x
  vals <- darvind(c(0.5, 1, 2), theta = 1)
  expect_length(vals, 3)
  expect_true(all(vals > 0))
})

test_that("parvind returns values in [0, 1]", {
  x_seq <- c(-1, 0, 0.001, 0.5, 1, 5, 50)
  vals <- parvind(x_seq, theta = 1)
  expect_true(all(vals >= 0))
  expect_true(all(vals <= 1))
})

test_that("Arvind survival function equals 1 - CDF", {
  x <- 1.5
  theta <- 3
  cdf <- parvind(x, theta)
  surv <- exp(-theta * x^2) / (1 + theta * x)
  expect_equal(cdf, 1 - surv, tolerance = 1e-10)
})

test_that("Hazard rate is positive", {
  theta <- 2
  x_seq <- seq(0.01, 2, length.out = 100)
  # h(x) = f(x) / S(x)
  f_vals <- darvind(x_seq, theta)
  s_vals <- parvind(x_seq, theta, lower.tail = FALSE)
  h_vals <- f_vals / s_vals
  
  expect_true(all(h_vals > 0))
  expect_true(all(!is.na(h_vals)))
})
