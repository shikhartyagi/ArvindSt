#' K-Fold and Rolling-Window Cross-Validation
#'
#' Performs k-fold cross-validation and optionally rolling-window
#' (expanding-window) cross-validation for an \code{ArvindFit} model.
#'
#' @param fit an object of class \code{"ArvindFit"}.
#' @param k_folds integer; number of cross-validation folds (default: 5).
#' @param rolling logical; if \code{TRUE} (default), also performs
#'   rolling-window cross-validation.
#' @param n0_frac numeric; fraction of data used as initial training window
#'   for rolling CV (default: 0.5).
#' @param seed integer; random seed for reproducibility (default: 42).
#'
#' @return A list with components:
#' \describe{
#'   \item{cv_rmse}{numeric vector of length \code{k_folds}; per-fold RMSE.}
#'   \item{cv_mae}{numeric vector of length \code{k_folds}; per-fold MAE.}
#'   \item{mean_cv_rmse}{numeric; average k-fold RMSE.}
#'   \item{mean_cv_mae}{numeric; average k-fold MAE.}
#'   \item{roll_rmse}{numeric; rolling-window RMSE (or \code{NA}).}
#' }
#'
#' @examples
#' dat <- simulate_arvind_data(n = 50, seed = 1)
#' m1 <- fit_rw1(Y ~ X1 + X2 + X3, dat, seed = 42)
#' cv <- cv_arvind(m1, k_folds = 3, rolling = FALSE, seed = 42)
#' cv$mean_cv_rmse
#'
#' @seealso [diagnostics_arvind()], [forecast_arvind()]
#'
#' @export
cv_arvind <- function(fit, k_folds = 5,
                      rolling = TRUE,
                      n0_frac = 0.5,
                      seed = 42) {
  set.seed(seed)
  Y <- fit$Y
  X <- fit$X
  n <- fit$n
  fold_ids <- sample(rep(seq_len(k_folds), length.out = n))
  cv_rmse <- cv_mae <- numeric(k_folds)

  for (f in seq_len(k_folds)) {
    test_idx  <- which(fold_ids == f)
    train_idx <- which(fold_ids != f)

    if (fit$model_type == "tvLM") {
      tv_cv <- tryCatch(
        tvReg::tvLM(fit$formula,
                    data = fit$data[train_idx, ],
                    bw = 0.5),
        error = function(e) NULL)
      if (!is.null(tv_cv)) {
        ntst <- length(test_idx)
        ncoef <- min(ntst, nrow(tv_cv$coefficients))
        pred_f <- rowSums(
          X[test_idx[seq_len(ncoef)], , drop = FALSE] *
          tv_cv$coefficients[seq_len(ncoef), , drop = FALSE])
        if (ncoef < ntst)
          pred_f <- c(pred_f, rep(mean(Y[train_idx]), ntst - ncoef))
      } else {
        pred_f <- rep(mean(Y[train_idx]), length(test_idx))
      }

    } else if (fit$model_type == "Mixed-Effects") {
      me_cv <- tryCatch(
        lme4::lmer(fit$formula_full,
          data = fit$data[train_idx, ],
          REML = FALSE),
        error = function(e) NULL)
      if (!is.null(me_cv))
        pred_f <- predict(me_cv,
          newdata = fit$data[test_idx, ],
          allow.new.levels = TRUE)
      else
        pred_f <- rep(mean(Y[train_idx]), length(test_idx))

    } else {
      fit_cv <- lm(Y[train_idx] ~ X[train_idx, -1])
      pred_f <- predict(fit_cv,
        newdata = data.frame(
          X[test_idx, -1, drop = FALSE]))
    }

    cv_rmse[f] <- sqrt(mean((Y[test_idx] - pred_f)^2))
    cv_mae[f]  <- mean(abs(Y[test_idx] - pred_f))
  }

  # Rolling-window CV
  roll_rmse <- NA
  if (rolling) {
    n0 <- max(20, floor(n * n0_frac))
    if (n0 < n) {
      roll_errs <- numeric(n - n0)
      for (m_idx in n0:(n - 1)) {
        train_idx <- seq_len(m_idx)
        fit_roll <- lm(Y[train_idx] ~ X[train_idx, -1])
        pred_r <- predict(fit_roll,
          newdata = data.frame(
            t(X[m_idx + 1, -1])))
        roll_errs[m_idx - n0 + 1] <- (Y[m_idx + 1] - pred_r)^2
      }
      roll_rmse <- sqrt(mean(roll_errs))
    }
  }

  list(cv_rmse = cv_rmse, cv_mae = cv_mae,
       mean_cv_rmse = mean(cv_rmse),
       mean_cv_mae  = mean(cv_mae),
       roll_rmse = roll_rmse)
}
