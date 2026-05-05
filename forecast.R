#' Monte Carlo Forecasting for Arvind Models
#'
#' Generates Monte Carlo forecasts with 80 percent and 95 percent prediction
#' intervals for any fitted \code{ArvindFit} model. Covariates are forecast
#' using SARIMA models (via the \pkg{forecast} package) if not supplied.
#'
#' @param fit an object of class \code{"ArvindFit"}.
#' @param newdata_sims optional named list of pre-computed covariate
#'   simulation matrices, each of dimension \code{h x nsim}.
#' @param h integer; forecast horizon in time steps (default: 120).
#' @param nsim integer; number of Monte Carlo replicates (default: 5000).
#' @param covariate_models optional list of fitted SARIMA models for
#'   covariates (auto-fitted if \code{NULL}).
#' @param seed integer; random seed for reproducibility (default: 123).
#'
#' @return A list with components:
#' \describe{
#'   \item{sims}{matrix (\code{h x nsim}); full simulation matrix.}
#'   \item{mean}{numeric vector length \code{h}; mean forecast.}
#'   \item{median}{numeric vector length \code{h}; median forecast.}
#'   \item{lo80}{numeric vector; lower 80 percent prediction interval.}
#'   \item{hi80}{numeric vector; upper 80 percent prediction interval.}
#'   \item{lo95}{numeric vector; lower 95 percent prediction interval.}
#'   \item{hi95}{numeric vector; upper 95 percent prediction interval.}
#' }
#'
#' @examples
#' \donttest{
#' dat <- simulate_arvind_data(n = 50, seed = 1)
#' m1 <- fit_rw1(Y ~ X1 + X2 + X3, dat, seed = 42)
#' fc <- forecast_arvind(m1, h = 12, nsim = 100, seed = 42)
#' head(fc$mean)
#' }
#'
#' @seealso [fit_rw1()], [diagnostics_arvind()], [cv_arvind()]
#'
#' @export
forecast_arvind <- function(fit, newdata_sims = NULL,
                            h = 120, nsim = 5000,
                            covariate_models = NULL,
                            seed = 123) {
  set.seed(seed)
  X <- fit$X
  Y <- fit$Y
  n <- fit$n
  p <- fit$p
  theta_m <- fit$theta
  sigma_m <- fit$sigma
  shift_m <- fit$shift

  # Auto-fit SARIMA for covariates if needed
  if (is.null(newdata_sims)) {
    covars <- colnames(X)[-1]
    cov_sims <- list()
    for (col in covars) {
      ts_col <- ts(fit$data[[col]], frequency = 12)
      fit_ar <- forecast::auto.arima(ts_col, seasonal = TRUE)
      sims_mat <- matrix(NA, nrow = h, ncol = nsim)
      for (s in seq_len(nsim))
        sims_mat[, s] <- simulate(fit_ar, nsim = h, future = TRUE)
      cov_sims[[col]] <- sims_mat
    }
    newdata_sims <- cov_sims
  }

  # Monte Carlo forecasting
  Y_fc_sims <- matrix(NA, nrow = h, ncol = nsim)
  covars <- names(newdata_sims)

  for (s in seq_len(nsim)) {
    X_fut <- cbind(1, sapply(covars,
      function(col) newdata_sims[[col]][, s]))
    eps_s <- rarvind(h, theta_m)

    if (fit$model_type == "RW1-approx") {
      beta_fc <- fit$beta_final
      for (tt in seq_len(h)) {
        w_fc <- rarvind_centred(p, fit$theta_innov)
        beta_fc <- beta_fc + fit$sigma_rw * w_fc
        Y_fc_sims[tt, s] <- sum(X_fut[tt, ] * beta_fc) +
          shift_m + sigma_m * eps_s[tt]
      }
    } else if (fit$model_type == "tvLM") {
      beta_last <- fit$tv_coefs[n, ]
      Y_fc_sims[, s] <- as.numeric(X_fut %*% beta_last) +
        shift_m + sigma_m * eps_s

    } else if (fit$model_type == "SIMEX") {
      Y_fc_sims[, s] <- as.numeric(X_fut %*% fit$beta) +
        shift_m + sigma_m * eps_s

    } else if (fit$model_type == "Mixed-Effects") {
      beta_fe <- lme4::fixef(fit$lme_model)
      Y_fc_sims[, s] <- as.numeric(X_fut %*% beta_fe) +
        shift_m + sigma_m * eps_s

    } else if (fit$model_type == "Regime-Switching") {
      tp <- fit$trans_probs
      curr_st <- tail(fit$states, 1)
      for (tt in seq_len(h)) {
        curr_st <- sample(seq_len(fit$nstates), 1,
                          prob = tp[curr_st, ])
        beta_fc <- fit$state_betas[[curr_st]]
        Y_fc_sims[tt, s] <- sum(X_fut[tt, ] * beta_fc) +
          shift_m + sigma_m * eps_s[tt]
      }
    }
  }

  list(
    sims   = Y_fc_sims,
    mean   = rowMeans(Y_fc_sims),
    median = apply(Y_fc_sims, 1, median),
    lo80   = apply(Y_fc_sims, 1, quantile, 0.10),
    hi80   = apply(Y_fc_sims, 1, quantile, 0.90),
    lo95   = apply(Y_fc_sims, 1, quantile, 0.025),
    hi95   = apply(Y_fc_sims, 1, quantile, 0.975)
  )
}
