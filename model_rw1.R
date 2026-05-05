#' Fit Random Walk on Coefficients Model (RW1-approx)
#'
#' Fits a stochastic regression model with time-varying coefficients evolving
#' as a random walk with Arvind-distributed innovations. The observation
#' errors also follow the Arvind distribution.
#'
#' @param formula an object of class \code{\link[stats]{formula}} specifying
#'   the model (e.g., \code{Y ~ X1 + X2}).
#' @param data a data frame containing the variables in the formula.
#' @param theta_innov positive numeric; the Arvind parameter for state
#'   innovations (default: 2.0).
#' @param rw_scale numeric; proportion of OLS coefficients used as innovation
#'   scale (default: 0.01).
#' @param seed integer; random seed for reproducibility (default: 42).
#'
#' @return An object of class \code{"ArvindFit"}, a list containing:
#' \describe{
#'   \item{model_type}{character; \code{"RW1-approx"}.}
#'   \item{fitted}{numeric vector; fitted values.}
#'   \item{residuals}{numeric vector; raw residuals.}
#'   \item{theta}{numeric; estimated Arvind parameter for residuals.}
#'   \item{sigma}{numeric; residual scale.}
#'   \item{shift}{numeric; shift applied to residuals.}
#'   \item{e_pos}{numeric vector; positive standardised residuals.}
#'   \item{negloglik}{numeric; negative log-likelihood.}
#'   \item{beta_t}{matrix; time-varying coefficient paths.}
#'   \item{beta_final}{numeric vector; final coefficient values.}
#'   \item{sigma_rw}{numeric vector; random walk innovation scales.}
#'   \item{theta_innov}{numeric; Arvind parameter used for innovations.}
#'   \item{n}{integer; number of observations.}
#'   \item{p}{integer; number of parameters.}
#'   \item{X}{matrix; design matrix.}
#'   \item{Y}{numeric vector; response variable.}
#'   \item{formula}{the model formula.}
#'   \item{data}{the input data frame.}
#' }
#'
#' @examples
#' dat <- simulate_arvind_data(n = 50, seed = 1)
#' m1 <- fit_rw1(Y ~ X1 + X2 + X3, dat, seed = 42)
#' m1$theta
#'
#' @seealso [diagnostics_arvind()], [forecast_arvind()], [cv_arvind()]
#'
#' @export
fit_rw1 <- function(formula, data,
                    theta_innov = 2.0,
                    rw_scale = 0.01,
                    seed = 42) {
  set.seed(seed)
  mf <- model.frame(formula, data)
  Y  <- model.response(mf)
  X  <- model.matrix(formula, data)
  n  <- length(Y)
  p  <- ncol(X)

  # OLS initialisation
  ols_fit  <- lm(formula, data = data)
  beta_ols <- coef(ols_fit)
  sigma_rw <- rw_scale * abs(beta_ols)

  # State propagation
  beta_t <- matrix(NA, nrow = n, ncol = p)
  Y_hat  <- numeric(n)
  beta_current <- beta_ols

  for (t in seq_len(n)) {
    if (t > 1) {
      w <- rarvind_centred(p, theta_innov)
      beta_current <- beta_current + sigma_rw * w
    }
    beta_t[t, ] <- beta_current
    Y_hat[t]    <- sum(X[t, ] * beta_current)
  }

  # Arvind residual fitting
  resid_raw <- Y - Y_hat
  arv <- make_arvind_resid(resid_raw, Y)

  result <- list(
    model_type  = "RW1-approx",
    fitted      = Y_hat,
    residuals   = resid_raw,
    theta       = arv$theta,
    sigma       = arv$sigma,
    shift       = arv$shift,
    e_pos       = arv$e_pos,
    negloglik   = arv$negloglik,
    beta_t      = beta_t,
    beta_final  = beta_current,
    sigma_rw    = sigma_rw,
    theta_innov = theta_innov,
    n = n, p = p,
    X = X, Y = Y,
    formula = formula,
    data = data
  )
  class(result) <- "ArvindFit"
  return(result)
}
