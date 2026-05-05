test_that("darvind returns correct values for positive x and theta", {
  # PDF should be positive for x > 0
  val <- darvind(1, theta = 1)
  expect_true(val > 0)

  # Vectorised input
  vals <- darvind(c(0.1, 0.5, 1, 2, 5), theta = 2)
  expect_length(vals, 5)
  expect_true(all(vals > 0))
})

test_that("darvind returns 0 for non-positive x", {
  expect_equal(darvind(0, theta = 1), 0)
  expect_equal(darvind(-1, theta = 1), 0)
  expect_equal(darvind(-5, theta = 2), 0)
})

test_that("darvind log=TRUE gives log of density", {
  x <- 1.5
  theta <- 2
  d_val <- darvind(x, theta, log = FALSE)
  ld_val <- darvind(x, theta, log = TRUE)
  expect_equal(log(d_val), ld_val, tolerance = 1e-12)
})

test_that("darvind PDF integrates to 1", {
  for (th in c(0.5, 1, 2, 5)) {
    integral <- integrate(function(x) darvind(x, th),
                          lower = 0, upper = Inf,
                          subdivisions = 2000)$value
    expect_equal(integral, 1, tolerance = 1e-6)
  }
})

test_that("darvind returns valid probabilities for a range of x", {
  theta <- 2
  x_seq <- seq(0, 5, length.out = 500)
  d_vals <- darvind(x_seq, theta)
  
  expect_true(all(!is.na(d_vals)))
  expect_true(all(d_vals >= 0))
})
