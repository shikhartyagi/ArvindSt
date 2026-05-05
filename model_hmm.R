#' Fit Regime-Switching Regression (HMM)
#'
#' Fits a hidden Markov model with state-dependent coefficients and
#' Arvind-distributed errors. The EM algorithm with forward-backward
#' recursions is used for parameter estimation, and the Viterbi algorithm
#' decodes the most likely state sequence.
#'
#' @param formula an object of class \code{\link[stats]{formula}}.
#' @param data a data frame containing the variables in the formula.
#' @param nstates integer; number of hidden states (default: 2).
#' @param seed integer; random seed for reproducibility (default: 42).
#'
#' @return An object of class \code{"ArvindFit"}, a list containing the
#'   same standard fields as [fit_rw1()], plus:
#' \describe{
#'   \item{hmm_fit}{the fitted \code{depmixS4} object.}
#'   \item{nstates}{integer; number of hidden states.}
#'   \item{states}{integer vector; Viterbi-decoded state sequence.}
#'   \item{trans_probs}{matrix; estimated transition probability matrix.}
#'   \item{state_betas}{list of numeric vectors; state-specific coefficients.}
#'   \item{state_sigmas}{numeric vector; state-specific standard deviations.}
#' }
#'
#' @examples
#' dat <- simulate_arvind_data(n = 50, seed = 1)
#' m5 <- fit_hmm(Y ~ X1 + X2 + X3, dat, nstates = 2, seed = 42)
#' m5$states
#' m5$states
#'
#' @seealso [diagnostics_arvind()], [forecast_arvind()], [cv_arvind()]
#'
#' @export
fit_hmm <- function(formula, data,
                    nstates = 2, seed = 42) {
  set.seed(seed)
  mf <- model.frame(formula, data)
  Y  <- model.response(mf)
  X  <- model.matrix(formula, data)
  n  <- length(Y)
  p  <- ncol(X)

  hmm_mod <- depmixS4::depmix(
    formula, data = data,
    nstates = nstates,
    family  = gaussian())
  hmm_fit <- depmixS4::fit(hmm_mod, verbose = FALSE)

  pars <- depmixS4::getpars(hmm_fit)
  K    <- nstates
  trans_start <- K + 1
  trans_probs <- matrix(
    pars[trans_start:(trans_start + K^2 - 1)],
    nrow = K, byrow = TRUE)

  state_betas  <- list()
  state_sigmas <- numeric(K)
  par_offset <- trans_start + K^2
  for (k in seq_len(K)) {
    state_betas[[k]] <- pars[
      par_offset:(par_offset + p - 1)]
    state_sigmas[k]  <- pars[par_offset + p]
    par_offset <- par_offset + p + 1
  }

  viterbi_states <- depmixS4::posterior(hmm_fit)$state

  Y_hat <- numeric(n)
  for (i in seq_len(n)) {
    k <- viterbi_states[i]
    Y_hat[i] <- sum(X[i, ] * state_betas[[k]])
  }

  resid_raw <- Y - Y_hat
  arv <- make_arvind_resid(resid_raw, Y)

  result <- list(
    model_type    = "Regime-Switching",
    fitted        = Y_hat,
    residuals     = resid_raw,
    theta         = arv$theta,
    sigma         = arv$sigma,
    shift         = arv$shift,
    e_pos         = arv$e_pos,
    negloglik     = arv$negloglik,
    hmm_fit       = hmm_fit,
    nstates       = nstates,
    states        = viterbi_states,
    trans_probs   = trans_probs,
    state_betas   = state_betas,
    state_sigmas  = state_sigmas,
    n = n, p = p,
    X = X, Y = Y,
    formula = formula,
    data = data
  )
  class(result) <- "ArvindFit"
  return(result)
}
