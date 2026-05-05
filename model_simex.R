#' Fit Simulation-Extrapolation (SIMEX) Model
#'
#' Fits a regression model correcting for measurement error attenuation using
#' the SIMEX algorithm with Arvind-distributed measurement noise and residuals.
#'
#' @param formula an object of class \code{\link[stats]{formula}}.
#' @param data a data frame containing the variables in the formula.
#' @param me_vars character vector of covariate names measured with error.
#'   If \code{NULL} (default), the first two term labels are used.
#' @param me_frac numeric; fraction of marginal variance used as measurement
#'   error variance (default: 0.05).
#' @param lambda_grid numeric vector; SIMEX lambda grid
#'   (default: \code{c(0.5, 1, 1.5, 2)}).
#' @param n_sim integer; number of SIMEX simulation replicates (default: 100).
#' @param theta_me positive numeric; Arvind parameter for measurement error
#'   (default: 2.0).
#' @param seed integer; random seed for reproducibility (default: 123).
#'
#' @return An object of class \code{"ArvindFit"}, a list containing the
#'   same standard fields as [fit_rw1()], plus:
#' \describe{
#'   \item{beta}{numeric vector; SIMEX-corrected coefficient estimates.}
#'   \item{simex_coefs}{matrix; coefficient estimates at each lambda level.}
#'   \item{lambda_grid}{numeric vector; the SIMEX lambda grid used.}
#'   \item{me_vars}{character vector; covariate names with measurement error.}
#'   \item{sigma2_me}{named numeric vector; measurement error variances.}
#' }
#'
#' @examples
#' dat <- simulate_arvind_data(n = 50, seed = 1)
#' m3 <- fit_simex(Y ~ X1 + X2 + X3, dat,
#'                 me_vars = c("X1", "X2"),
#'                 n_sim = 20, seed = 123)
#' m3$beta
#' m3$beta
#'
#' @seealso [diagnostics_arvind()], [forecast_arvind()], [cv_arvind()]
#'
#' @export
fit_simex <- function(formula, data,
                      me_vars = NULL,
                      me_frac = 0.05,
                      lambda_grid = c(0.5, 1, 1.5, 2),
                      n_sim = 100,
                      theta_me = 2.0,
                      seed = 123) {
  set.seed(seed)
  mf <- model.frame(formula, data)
  Y  <- model.response(mf)
  X  <- model.matrix(formula, data)
  n  <- length(Y)
  p  <- ncol(X)

  if (is.null(me_vars)) {
    cov_names <- attr(terms(formula), "term.labels")
    me_vars   <- cov_names[seq_len(min(2, length(cov_names)))]
  }
  sigma2_me <- sapply(me_vars,
    function(v) var(data[[v]]) * me_frac)

  # Simulation step
  simex_coefs <- matrix(NA, nrow = length(lambda_grid), ncol = p)
  colnames(simex_coefs) <- colnames(X)

  for (l_idx in seq_along(lambda_grid)) {
    lam <- lambda_grid[l_idx]
    coef_sims <- matrix(0, nrow = n_sim, ncol = p)
    for (s in seq_len(n_sim)) {
      data_noisy <- data
      for (v in me_vars) {
        noise <- rarvind_centred(n, theta_me) *
                 sqrt(lam * sigma2_me[v])
        data_noisy[[v]] <- data[[v]] + noise
      }
      fit_s <- lm(formula, data = data_noisy)
      coef_sims[s, ] <- coef(fit_s)
    }
    simex_coefs[l_idx, ] <- colMeans(coef_sims)
  }

  # Extrapolation step
  beta_simex <- numeric(p)
  for (j in seq_len(p)) {
    qfit <- lm(simex_coefs[, j] ~
      lambda_grid + I(lambda_grid^2))
    beta_simex[j] <- predict(qfit,
      newdata = data.frame(lambda_grid = -1))
  }
  names(beta_simex) <- colnames(X)

  Y_hat <- as.numeric(X %*% beta_simex)
  resid_raw <- Y - Y_hat
  arv <- make_arvind_resid(resid_raw, Y)

  result <- list(
    model_type   = "SIMEX",
    fitted       = Y_hat,
    residuals    = resid_raw,
    theta        = arv$theta,
    sigma        = arv$sigma,
    shift        = arv$shift,
    e_pos        = arv$e_pos,
    negloglik    = arv$negloglik,
    beta         = beta_simex,
    simex_coefs  = simex_coefs,
    lambda_grid  = lambda_grid,
    me_vars      = me_vars,
    sigma2_me    = sigma2_me,
    n = n, p = p,
    X = X, Y = Y,
    formula = formula,
    data = data
  )
  class(result) <- "ArvindFit"
  return(result)
}
