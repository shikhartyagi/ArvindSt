test_that("parvind returns 0 at x = 0", {
  expect_equal(parvind(0, theta = 1), 0)
  expect_equal(parvind(0, theta = 5), 0)
})

test_that("parvind approaches 1 for large x", {
  val <- parvind(100, theta = 1)
  expect_true(val > 0.999)
})

test_that("parvind is monotone increasing", {
  x_seq <- seq(0.01, 10, length.out = 100)
  cdf_vals <- parvind(x_seq, theta = 2)
  diffs <- diff(cdf_vals)
  expect_true(all(diffs >= 0))
})

test_that("parvind lower.tail=FALSE is complement", {
  x <- 1.5
  theta <- 2
  lower <- parvind(x, theta, lower.tail = TRUE)
  upper <- parvind(x, theta, lower.tail = FALSE)
  expect_equal(lower + upper, 1, tolerance = 1e-12)
})

test_that("parvind is consistent with darvind via numerical integration", {
  theta <- 2
  x <- 1.5
  cdf_numerical <- integrate(function(t) darvind(t, theta),
                              lower = 0, upper = x)$value
  cdf_analytical <- parvind(x, theta)
  expect_equal(cdf_numerical, cdf_analytical, tolerance = 1e-6)
})

test_that("qarvind inverts parvind", {
  theta <- 2
  for (p in c(0.1, 0.25, 0.5, 0.75, 0.9)) {
    q <- qarvind(p, theta)
    p_back <- parvind(q, theta)
    expect_equal(p_back, p, tolerance = 1e-6)
  }
})

test_that("qarvind handles boundary values", {
  expect_equal(qarvind(0, theta = 1), 0)
  expect_equal(qarvind(1, theta = 1), Inf)
})
