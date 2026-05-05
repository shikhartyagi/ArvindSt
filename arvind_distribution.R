#' Arvind Distribution Density Function
#'
#' Computes the probability density function (PDF) of the Arvind distribution.
#'
#' The Arvind distribution with parameter \eqn{\theta > 0} has PDF
#' \deqn{f(x; \theta) = \frac{\theta(1 + 2x + 2\theta x^2)}{(1 + \theta x)^2}
#'   \exp(-\theta x^2), \quad x > 0.}
#'
#' @param x numeric vector of quantiles.
#' @param theta positive numeric scalar; the distribution parameter.
#' @param log logical; if \code{TRUE}, log-density is returned. Default
#'   \code{FALSE}.
#'
#' @return A numeric vector of density values (or log-density values when
#'   \code{log = TRUE}).
#'
#' @examples
#' # Evaluate the PDF at several points
#' darvind(c(0.5, 1, 2), theta = 1)
#'
#' # Log-density
#' darvind(1, theta = 2, log = TRUE)
#'
#' # Returns 0 for x <= 0
#' darvind(-1, theta = 1)
#'
#' @export
darvind <- function(x, theta, log = FALSE) {
  val <- rep(-Inf, length(x))
  idx <- (x > 0) & (theta > 0)
  if (any(idx)) {
    xx <- x[idx]
    th <- theta
    num <- th * (1 + 2 * xx + 2 * th * xx^2)
    den <- (1 + th * xx)^2
    lv  <- log(num) - log(den) - th * xx^2
    val[idx] <- lv
  }
  if (log) return(val) else return(exp(val))
}

#' Arvind Distribution Function (CDF)
#'
#' Computes the cumulative distribution function (CDF) of the Arvind
#' distribution.
#'
#' The CDF is given by
#' \deqn{F(x; \theta) = 1 - \frac{1}{1 + \theta x} \exp(-\theta x^2),
#'   \quad x > 0.}
#'
#' @param q numeric vector of quantiles.
#' @param theta positive numeric scalar; the distribution parameter.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#'   \eqn{P(X \le q)}{P(X <= q)}; otherwise \eqn{P(X > q)}.
#'
#' @return A numeric vector of probabilities.
#'
#' @examples
#' parvind(1, theta = 1)
#' parvind(c(0.5, 1, 2), theta = 2)
#' parvind(1, theta = 1, lower.tail = FALSE)
#'
#' @export
parvind <- function(q, theta, lower.tail = TRUE) {
  val <- rep(0, length(q))
  idx <- (q > 0) & (theta > 0)
  if (any(idx)) {
    qq <- q[idx]
    th <- theta
    surv <- exp(-th * qq^2) / (1 + th * qq)
    val[idx] <- 1 - surv
  }
  if (lower.tail) return(val) else return(1 - val)
}

#' Arvind Distribution Quantile Function
#'
#' Computes quantiles of the Arvind distribution by numerical inversion
#' of the CDF using \code{\link[stats]{uniroot}}.
#'
#' @param p numeric vector of probabilities (\eqn{0 \le p \le 1}).
#' @param theta positive numeric scalar; the distribution parameter.
#'
#' @return A numeric vector of quantiles.
#'
#' @examples
#' qarvind(0.5, theta = 1)
#' qarvind(c(0.25, 0.5, 0.75), theta = 2)
#'
#' @export
qarvind <- function(p, theta) {
  sapply(p, function(pp) {
    if (pp <= 0) return(0)
    if (pp >= 1) return(Inf)
    fn <- function(x) parvind(x, theta) - pp
    tryCatch(
      uniroot(fn, lower = 1e-10,
              upper = 100 / theta, tol = 1e-12)$root,
      error = function(e) {
        tryCatch(
          uniroot(fn, lower = 1e-10,
                  upper = 1000 / theta,
                  tol = 1e-10)$root,
          error = function(e2) NA_real_)
      })
  })
}

#' Random Generation from the Arvind Distribution
#'
#' Generates random variates from the Arvind distribution using a rejection
#' sampling algorithm with a half-normal proposal distribution.
#'
#' @param n positive integer; number of random variates to generate.
#' @param theta positive numeric scalar; the distribution parameter.
#'
#' @return A numeric vector of length \code{n} containing positive random
#'   variates.
#'
#' @examples
#' set.seed(42)
#' x <- rarvind(100, theta = 1)
#' summary(x)
#'
#' @export
rarvind <- function(n, theta) {
  samples <- numeric(0)
  mode_x <- tryCatch(
    optimize(function(x) darvind(x, theta),
             interval = c(0, 10 / sqrt(theta)),
             maximum = TRUE)$maximum,
    error = function(e) 0.5 / sqrt(theta))
  f_max   <- darvind(mode_x, theta) * 1.5
  prop_sd <- max(1 / sqrt(theta), 0.1)
  while (length(samples) < n) {
    needed     <- (n - length(samples)) * 3
    candidates <- abs(rnorm(needed, 0, prop_sd))
    f_vals <- darvind(candidates, theta)
    g_vals <- f_max * dnorm(candidates, 0, prop_sd) *
              2 / dnorm(0, 0, prop_sd)
    g_vals <- pmax(g_vals,
                   f_max * exp(-theta * candidates^2))
    u      <- runif(needed)
    accept <- u <= f_vals / pmax(g_vals, f_vals * 1.01)
    samples <- c(samples, candidates[accept])
  }
  return(samples[1:n])
}

#' Centred Random Generation from the Arvind Distribution
#'
#' Generates centred Arvind variates with approximately zero mean, suitable
#' for use as error terms and innovation terms in stochastic regression models.
#'
#' The centred variate is computed as \eqn{\tilde{\varepsilon} =
#' \varepsilon - \mu_A(\theta)}, where \eqn{\varepsilon \sim
#' \mathrm{Arvind}(\theta)} and \eqn{\mu_A(\theta)} is the mean of the
#' Arvind distribution.
#'
#' @param n positive integer; number of random variates to generate.
#' @param theta positive numeric scalar; the distribution parameter.
#'
#' @return A numeric vector of length \code{n} with approximately zero mean.
#'
#' @examples
#' set.seed(42)
#' eps <- rarvind_centred(1000, theta = 2)
#' mean(eps)  # approximately 0
#'
#' @export
rarvind_centred <- function(n, theta) {
  raw  <- rarvind(n, theta)
  mu_a <- tryCatch(
    integrate(function(x) x * darvind(x, theta),
              0, Inf)$value,
    error = function(e) mean(raw))
  return(raw - mu_a)
}
