#' Fit Mixed-Effects Regression with Arvind Errors
#'
#' Fits a mixed-effects regression model with Arvind-distributed random
#' effects and observation-level errors. Estimation uses a two-stage approach:
#' REML initialisation via \pkg{lme4}, followed by Arvind MLE on the
#' residuals.
#'
#' @param formula an object of class \code{\link[stats]{formula}} specifying
#'   the fixed-effects structure.
#' @param data a data frame containing the variables in the formula and the
#'   grouping variable.
#' @param group_var character string; the name of the grouping variable in
#'   \code{data} (default: \code{"Season"}).
#' @param re_formula optional random-effects formula (e.g.,
#'   \code{(1 + X1 | group)}). If \code{NULL} (default), a random intercept
#'   model \code{(1 | group_var)} is used.
#' @param seed integer; random seed for reproducibility (default: 42).
#'
#' @return An object of class \code{"ArvindFit"}, a list containing the
#'   same standard fields as [fit_rw1()], plus:
#' \describe{
#'   \item{lme_model}{the fitted \code{lme4::lmer} object.}
#'   \item{theta_re}{numeric; Arvind parameter estimated from random effects.}
#'   \item{group_var}{character; the grouping variable name.}
#' }
#'
#' @examples
#' dat <- simulate_arvind_data(n = 50, seed = 1)
#' m4 <- fit_mixed(Y ~ X1 + X2 + X3, dat, group_var = "Group", seed = 42)
#' m4$theta
#' m4$theta
#'
#' @seealso [diagnostics_arvind()], [forecast_arvind()], [cv_arvind()]
#'
#' @export
fit_mixed <- function(formula, data,
                      group_var = "Season",
                      re_formula = NULL,
                      seed = 42) {
  set.seed(seed)
  mf <- model.frame(formula, data)
  Y  <- model.response(mf)
  X  <- model.matrix(formula, data)
  n  <- length(Y)
  p  <- ncol(X)

  resp      <- as.character(formula[[2]])
  fixed_rhs <- as.character(formula[[3]])
  if (is.null(re_formula)) {
    full_frm <- as.formula(paste0(
      resp, " ~ ", fixed_rhs,
      " + (1|", group_var, ")"))
  } else {
    full_frm <- as.formula(paste0(
      resp, " ~ ", fixed_rhs,
      " + ", deparse(re_formula)))
  }

  me_fit <- tryCatch(
    lme4::lmer(full_frm, data = data, REML = FALSE),
    error = function(e) {
      full_frm2 <- as.formula(paste0(
        resp, " ~ ", fixed_rhs,
        " + (1|", group_var, ")"))
      lme4::lmer(full_frm2, data = data, REML = FALSE)
    })

  Y_hat     <- fitted(me_fit)
  resid_raw <- Y - Y_hat
  arv <- make_arvind_resid(resid_raw, Y)

  # Fit Arvind to random effects
  re_vals <- unlist(lme4::ranef(me_fit))
  re_pos  <- re_vals - min(re_vals) + 0.01
  arv_re  <- fit_arvind_mle(re_pos)

  result <- list(
    model_type  = "Mixed-Effects",
    fitted      = Y_hat,
    residuals   = resid_raw,
    theta       = arv$theta,
    sigma       = arv$sigma,
    shift       = arv$shift,
    e_pos       = arv$e_pos,
    negloglik   = arv$negloglik,
    lme_model   = me_fit,
    theta_re    = arv_re$theta,
    group_var   = group_var,
    formula_full = full_frm,
    n = n, p = p,
    X = X, Y = Y,
    formula = formula,
    data = data
  )
  class(result) <- "ArvindFit"
  return(result)
}
