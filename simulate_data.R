#' Generate Simulated Data for Examples
#'
#' Creates a small simulated dataset that mimics the structure needed for
#' demonstrating the \pkg{ArvindSt} model-fitting functions. Useful for
#' examples and testing.
#'
#' @param n integer; number of observations to generate (default: 60).
#' @param seed integer; random seed for reproducibility (default: 42).
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{Y}{numeric; simulated response variable.}
#'   \item{X1}{numeric; first covariate.}
#'   \item{X2}{numeric; second covariate.}
#'   \item{X3}{numeric; third covariate.}
#'   \item{Group}{factor; grouping variable with 4 levels.}
#' }
#'
#' @examples
#' dat <- simulate_arvind_data(n = 50, seed = 1)
#' head(dat)
#'
#' @export
simulate_arvind_data <- function(n = 60, seed = 42) {
  set.seed(seed)
  X1 <- rnorm(n, 50, 10)
  X2 <- rnorm(n, 30, 5)
  X3 <- rnorm(n, 20, 8)
  eps <- rarvind(n, theta = 2)
  mu_eps <- tryCatch(
    integrate(function(x) x * darvind(x, 2), 0, Inf)$value,
    error = function(e) mean(eps))
  Y <- 100 + 2 * X1 - 1.5 * X2 + 0.8 * X3 + 10 * (eps - mu_eps)
  Group <- factor(rep(c("A", "B", "C", "D"), length.out = n))
  data.frame(Y = Y, X1 = X1, X2 = X2, X3 = X3, Group = Group)
}
