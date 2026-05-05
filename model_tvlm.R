#' Fit Time-Varying Coefficient Linear Model (tvLM)
#'
#' Fits a time-varying coefficient linear model using kernel-weighted
#' least squares (via the \pkg{tvReg} package) with Arvind-distributed
#' residuals.
#'
#' @param formula an object of class \code{\link[stats]{formula}}.
#' @param data a data frame containing the variables in the formula.
#' @param bw numeric or \code{NULL}; the bandwidth for kernel smoothing.
#'   If \code{NULL} (default), bandwidth is selected automatically via
#'   leave-one-out cross-validation.
#' @param seed integer; random seed for reproducibility (default: 42).
#'
#' @return An object of class \code{"ArvindFit"}, a list containing the
#'   same standard fields as [fit_rw1()], plus:
#' \describe{
#'   \item{tv_coefs}{matrix; time-varying coefficient estimates.}
#'   \item{tv_fit}{the fitted \code{tvReg::tvLM} object.}
#' }
#'
#' @examples
#' dat <- simulate_arvind_data(n = 50, seed = 1)
#' m2 <- fit_tvlm(Y ~ X1 + X2 + X3, dat, bw = 0.5, seed = 42)
#' m2$theta
#'
#' @seealso [diagnostics_arvind()], [forecast_arvind()], [cv_arvind()]
#'
#' @export
fit_tvlm <- function(formula, data,
                     bw = NULL, seed = 42) {
  set.seed(seed)
  mf <- model.frame(formula, data)
  Y  <- model.response(mf)
  X  <- model.matrix(formula, data)
  n  <- length(Y)
  p  <- ncol(X)

  tv_fit <- tryCatch(
    tvReg::tvLM(formula, data = data, bw = bw),
    error = function(e)
      tvReg::tvLM(formula, data = data, bw = 0.5))

  tv_coefs  <- tv_fit$coefficients
  Y_hat     <- rowSums(X * tv_coefs)
  resid_raw <- Y - Y_hat
  arv <- make_arvind_resid(resid_raw, Y)

  result <- list(
    model_type = "tvLM",
    fitted     = Y_hat,
    residuals  = resid_raw,
    theta      = arv$theta,
    sigma      = arv$sigma,
    shift      = arv$shift,
    e_pos      = arv$e_pos,
    negloglik  = arv$negloglik,
    tv_coefs   = tv_coefs,
    tv_fit     = tv_fit,
    n = n, p = p,
    X = X, Y = Y,
    formula = formula,
    data = data
  )
  class(result) <- "ArvindFit"
  return(result)
}
