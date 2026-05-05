#' Maximum Likelihood Estimation for the Arvind Distribution
#'
#' Fits the Arvind distribution to a vector of positive observations by
#' maximum likelihood. Optimisation is performed on the log-scale via the
#' Brent method.
#'
#' @param e_pos numeric vector of strictly positive observations.
#'
#' @return A list with components:
#' \describe{
#'   \item{theta}{numeric; the MLE of theta.}
#'   \item{negloglik}{numeric; the minimised negative log-likelihood.}
#' }
#'
#' @examples
#' set.seed(42)
#' x <- rarvind(200, theta = 2)
#' fit_arvind_mle(x)
#'
#' @export
fit_arvind_mle <- function(e_pos) {
  if (length(e_pos) < 3 || any(!is.finite(e_pos)))
    return(list(theta = NA, negloglik = NA))
  neg_ll <- function(log_th) {
    th <- exp(log_th)
    ne <- length(e_pos)
    ll <- ne * log(th) +
      sum(log(pmax(1e-300,
                   1 + 2 * e_pos + 2 * th * e_pos^2))) -
      2 * sum(log(1 + th * e_pos)) - th * sum(e_pos^2)
    if (!is.finite(ll)) return(1e15)
    return(-ll)
  }
  opt <- tryCatch(
    optim(log(1 / mean(e_pos)), neg_ll,
          method = "Brent", lower = -10, upper = 15),
    error = function(e)
      list(par = 0, value = 1e15))
  list(theta = exp(opt$par), negloglik = opt$value)
}

#' Transform Residuals for Arvind Fitting
#'
#' Transforms raw residuals to positive values suitable for fitting the Arvind
#' distribution by shifting and standardising.
#'
#' @param resid_raw numeric vector of raw residuals.
#' @param Y_ref numeric vector of observed response values (used for scaling).
#'
#' @return A list with components:
#' \describe{
#'   \item{shift}{numeric; the shift applied.}
#'   \item{sigma}{numeric; the standard deviation used for standardisation.}
#'   \item{e_pos}{numeric vector; positive standardised residuals.}
#'   \item{theta}{numeric; MLE of the Arvind parameter.}
#'   \item{negloglik}{numeric; negative log-likelihood at the MLE.}
#' }
#'
#' @keywords internal
make_arvind_resid <- function(resid_raw, Y_ref) {
  shift_val <- min(resid_raw) - 0.01 * sd(Y_ref)
  resid_pos <- resid_raw - shift_val
  sigma_s   <- sd(resid_pos)
  e_std     <- resid_pos / sigma_s
  e_pos     <- e_std[e_std > 0]
  arv       <- fit_arvind_mle(e_pos)
  list(shift = shift_val, sigma = sigma_s,
       e_pos = e_pos,
       theta = arv$theta, negloglik = arv$negloglik)
}

#' Mean of the Arvind Distribution
#'
#' Computes the theoretical mean of the Arvind distribution with parameter
#' \code{theta} by numerical integration.
#'
#' @param theta positive numeric scalar; the distribution parameter.
#'
#' @return A numeric scalar giving the theoretical mean, or \code{NA} if
#'   integration fails.
#'
#' @examples
#' arvind_mean_fn(1)
#' arvind_mean_fn(2)
#'
#' @export
arvind_mean_fn <- function(theta) {
  tryCatch(integrate(function(x) x * darvind(x, theta),
           0, Inf, subdivisions = 2000,
           rel.tol = 1e-10)$value,
           error = function(e) NA)
}

#' Variance of the Arvind Distribution
#'
#' Computes the theoretical variance of the Arvind distribution with parameter
#' \code{theta} by numerical integration.
#'
#' @param theta positive numeric scalar; the distribution parameter.
#'
#' @return A numeric scalar giving the theoretical variance, or \code{NA} if
#'   integration fails.
#'
#' @examples
#' arvind_var_fn(1)
#' arvind_var_fn(2)
#'
#' @export
arvind_var_fn <- function(theta) {
  mu  <- arvind_mean_fn(theta)
  mu2 <- tryCatch(
    integrate(function(x) x^2 * darvind(x, theta),
              0, Inf, subdivisions = 2000,
              rel.tol = 1e-10)$value,
    error = function(e) NA)
  if (is.na(mu) || is.na(mu2)) return(NA)
  return(mu2 - mu^2)
}
