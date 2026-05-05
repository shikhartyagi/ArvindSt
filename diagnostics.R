#' Goodness-of-Fit Diagnostics for Arvind Models
#'
#' Computes 21 goodness-of-fit metrics for any fitted \code{ArvindFit}
#' object, including MSE, RMSE, MAE, MAPE, R-squared, AIC, BIC,
#' Kolmogorov-Smirnov test, Anderson-Darling statistic, and more.
#'
#' @param fit an object of class \code{"ArvindFit"} returned by any of the
#'   model-fitting functions.
#'
#' @return A data frame with one row and 21 columns of diagnostics metrics.
#'   See Details for the full list.
#'
#' @details
#' The following metrics are computed:
#' \describe{
#'   \item{Model}{character; the model type.}
#'   \item{MSE}{Mean Squared Error.}
#'   \item{RMSE}{Root Mean Squared Error.}
#'   \item{MAE}{Mean Absolute Error.}
#'   \item{MAPE}{Mean Absolute Percentage Error.}
#'   \item{R2}{R-squared.}
#'   \item{AdjR2}{Adjusted R-squared.}
#'   \item{AIC}{Akaike Information Criterion.}
#'   \item{AICc}{Corrected AIC.}
#'   \item{BIC}{Bayesian Information Criterion.}
#'   \item{LogLik}{Log-likelihood at the MLE.}
#'   \item{Bias}{Mean residual.}
#'   \item{MASE}{Mean Absolute Scaled Error.}
#'   \item{DW}{Durbin-Watson statistic.}
#'   \item{LjungBox_stat}{Ljung-Box test statistic.}
#'   \item{LjungBox_p}{Ljung-Box p-value.}
#'   \item{Theta}{Estimated Arvind parameter.}
#'   \item{KS_stat}{Kolmogorov-Smirnov test statistic.}
#'   \item{KS_pvalue}{Kolmogorov-Smirnov p-value.}
#'   \item{AD_stat}{Anderson-Darling test statistic.}
#'   \item{CvM_stat}{Cramer-von Mises test statistic.}
#' }
#'
#' @examples
#' dat <- simulate_arvind_data(n = 50, seed = 1)
#' m1 <- fit_rw1(Y ~ X1 + X2 + X3, dat, seed = 42)
#' diagnostics_arvind(m1)
#'
#' @export
diagnostics_arvind <- function(fit) {
  Y <- fit$Y
  Y_hat <- fit$fitted
  resid <- fit$residuals
  theta <- fit$theta
  e_pos <- fit$e_pos
  n <- fit$n
  k <- fit$p + 2
  negloglik <- fit$negloglik

  MSE  <- mean(resid^2)
  RMSE <- sqrt(MSE)
  MAE  <- mean(abs(resid))
  MAPE <- mean(abs(resid / Y)) * 100
  SS_res <- sum(resid^2)
  SS_tot <- sum((Y - mean(Y))^2)
  R2    <- 1 - SS_res / SS_tot
  AdjR2 <- 1 - (1 - R2) * (n - 1) / (n - k - 1)
  ll    <- -negloglik
  AIC_v <- -2 * ll + 2 * k
  AICc  <- AIC_v + (2 * k * (k + 1)) / max(n - k - 1, 1)
  BIC_v <- -2 * ll + k * log(n)
  Bias  <- mean(resid)
  naive_err <- mean(abs(diff(Y)))
  MASE <- MAE / max(naive_err, 1e-10)
  DW   <- sum(diff(resid)^2) / sum(resid^2)
  LB   <- Box.test(resid,
    lag = min(10, n - 1), type = "Ljung-Box")
  KS   <- suppressWarnings(ks.test(e_pos,
    function(x) parvind(x, theta)))

  # Anderson-Darling
  u <- sort(parvind(e_pos, theta))
  m <- length(u)
  u <- pmin(pmax(u, 1e-10), 1 - 1e-10)
  AD <- -m - (1 / m) * sum(
    (2 * (seq_len(m)) - 1) * (log(u) + log(1 - rev(u))))

  # Cramer-von Mises
  CvM <- sum((u - (2 * (seq_len(m)) - 1) / (2 * m))^2) + 1 / (12 * m)

  data.frame(
    Model = fit$model_type,
    MSE = MSE, RMSE = RMSE, MAE = MAE,
    MAPE = MAPE, R2 = R2, AdjR2 = AdjR2,
    AIC = AIC_v, AICc = AICc, BIC = BIC_v,
    LogLik = ll, Bias = Bias, MASE = MASE,
    DW = DW,
    LjungBox_stat = as.numeric(LB$statistic),
    LjungBox_p = LB$p.value,
    Theta = theta,
    KS_stat = as.numeric(KS$statistic),
    KS_pvalue = KS$p.value,
    AD_stat = AD, CvM_stat = CvM,
    stringsAsFactors = FALSE)
}

#' Diagnostic Plots for Arvind Models
#'
#' Generates up to 25 diagnostic plots for a fitted \code{ArvindFit} object,
#' including observed vs fitted, residual histogram with Arvind density overlay,
#' Q-Q plot, ACF, ECDF comparison, and more.
#'
#' @param fit an object of class \code{"ArvindFit"}.
#' @param output_dir character; directory where plots are saved. Defaults to
#'   a temporary directory.
#' @param prefix character or \code{NULL}; prefix for plot filenames. If
#'   \code{NULL}, derived from the model type.
#'
#' @return The \code{fit} object is returned invisibly.
#'
#' @examples
#' \donttest{
#' dat <- simulate_arvind_data(n = 50, seed = 1)
#' m1 <- fit_rw1(Y ~ X1 + X2 + X3, dat, seed = 42)
#' plot_arvind(m1, output_dir = tempdir())
#' }
#'
#' @export
plot_arvind <- function(fit,
                        output_dir = tempdir(),
                        prefix = NULL) {
  if (is.null(prefix))
    prefix <- gsub("[^A-Za-z0-9]", "_", fit$model_type)

  Y     <- fit$Y
  Y_hat <- fit$fitted
  resid <- fit$residuals
  theta <- fit$theta
  e_pos <- fit$e_pos
  n     <- fit$n

  if (!dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE)

  theme_paper <- ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", hjust = 0.5, size = 14),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom")

  # --- Plot 1: Observed vs Fitted ---
  grDevices::png(file.path(output_dir, paste0(prefix, "_01_Obs_vs_Fit.png")),
      width = 800, height = 600)
  df1 <- data.frame(Index = seq_len(n), Observed = Y, Fitted = Y_hat)
  p1 <- ggplot2::ggplot(df1, ggplot2::aes(x = .data$Index)) +
    ggplot2::geom_line(ggplot2::aes(y = .data$Observed, colour = "Observed")) +
    ggplot2::geom_line(ggplot2::aes(y = .data$Fitted, colour = "Fitted"),
                       linetype = "dashed") +
    ggplot2::labs(title = paste(fit$model_type, "- Observed vs Fitted"),
                  x = "Time Index", y = "Value", colour = "") +
    theme_paper
  print(p1)
  grDevices::dev.off()

  # --- Plot 2: Residual Histogram with Arvind Density ---
  grDevices::png(file.path(output_dir, paste0(prefix, "_02_Resid_Hist.png")),
      width = 800, height = 600)
  x_grid <- seq(min(e_pos) * 0.8, max(e_pos) * 1.2, length.out = 200)
  arvind_dens <- darvind(x_grid, theta)
  df2 <- data.frame(e_pos = e_pos)
  p2 <- ggplot2::ggplot(df2, ggplot2::aes(x = .data$e_pos)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                            bins = 20, fill = "steelblue", alpha = 0.6) +
    ggplot2::geom_line(data = data.frame(x = x_grid, y = arvind_dens),
                       ggplot2::aes(x = .data$x, y = .data$y),
                       colour = "red", linewidth = 1) +
    ggplot2::labs(title = paste(fit$model_type, "- Residual Histogram"),
                  x = "Standardised Residual", y = "Density") +
    theme_paper
  print(p2)
  grDevices::dev.off()

  # --- Plot 3: Q-Q Plot (Arvind) ---
  grDevices::png(file.path(output_dir, paste0(prefix, "_03_QQ_Arvind.png")),
      width = 800, height = 600)
  e_sorted <- sort(e_pos)
  m <- length(e_sorted)
  theo_q <- qarvind((seq_len(m) - 0.5) / m, theta)
  df3 <- data.frame(Theoretical = theo_q, Sample = e_sorted)
  p3 <- ggplot2::ggplot(df3, ggplot2::aes(x = .data$Theoretical,
                                           y = .data$Sample)) +
    ggplot2::geom_point(colour = "steelblue") +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "red",
                         linetype = "dashed") +
    ggplot2::labs(title = paste(fit$model_type, "- Q-Q Plot (Arvind)"),
                  x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_paper
  print(p3)
  grDevices::dev.off()

  # --- Plot 4: ACF of Residuals ---
  grDevices::png(file.path(output_dir, paste0(prefix, "_04_ACF.png")),
      width = 800, height = 600)
  acf_vals <- stats::acf(resid, plot = FALSE, lag.max = min(20, n - 1))
  df4 <- data.frame(Lag = acf_vals$lag[-1], ACF = acf_vals$acf[-1])
  p4 <- ggplot2::ggplot(df4, ggplot2::aes(x = .data$Lag, y = .data$ACF)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue", width = 0.5) +
    ggplot2::geom_hline(yintercept = c(-1.96, 1.96) / sqrt(n),
                        linetype = "dashed", colour = "red") +
    ggplot2::labs(title = paste(fit$model_type, "- ACF of Residuals"),
                  x = "Lag", y = "ACF") +
    theme_paper
  print(p4)
  grDevices::dev.off()

  # --- Plot 7: Actual vs Predicted Scatter ---
  grDevices::png(file.path(output_dir, paste0(prefix, "_07_Scatter.png")),
      width = 800, height = 600)
  df7 <- data.frame(Actual = Y, Predicted = Y_hat)
  p7 <- ggplot2::ggplot(df7, ggplot2::aes(x = .data$Actual,
                                           y = .data$Predicted)) +
    ggplot2::geom_point(colour = "steelblue") +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "red",
                         linetype = "dashed") +
    ggplot2::labs(title = paste(fit$model_type, "- Actual vs Predicted"),
                  x = "Actual", y = "Predicted") +
    theme_paper
  print(p7)
  grDevices::dev.off()

  # --- Plot 8: Residuals vs Fitted ---
  grDevices::png(file.path(output_dir, paste0(prefix, "_08_Resid_vs_Fit.png")),
      width = 800, height = 600)
  df8 <- data.frame(Fitted = Y_hat, Residuals = resid)
  p8 <- ggplot2::ggplot(df8, ggplot2::aes(x = .data$Fitted,
                                           y = .data$Residuals)) +
    ggplot2::geom_point(colour = "steelblue") +
    ggplot2::geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
    ggplot2::labs(title = paste(fit$model_type, "- Residuals vs Fitted"),
                  x = "Fitted Values", y = "Residuals") +
    theme_paper
  print(p8)
  grDevices::dev.off()

  # --- Plot 10: Arvind PDF Family ---
  grDevices::png(file.path(output_dir, paste0(prefix, "_10_PDF_Family.png")),
      width = 800, height = 600)
  thetas <- c(0.5, 1, 2, 5, 10)
  x_seq <- seq(0.01, 5, length.out = 300)
  df10 <- do.call(rbind, lapply(thetas, function(th) {
    data.frame(x = x_seq, y = darvind(x_seq, th),
               theta = paste0("theta = ", th))
  }))
  p10 <- ggplot2::ggplot(df10, ggplot2::aes(x = .data$x, y = .data$y,
                                              colour = .data$theta)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::labs(title = "Arvind Distribution PDF Family",
                  x = "x", y = "f(x)", colour = "") +
    theme_paper
  print(p10)
  grDevices::dev.off()

  # --- Plot 11: Arvind CDF Family ---
  grDevices::png(file.path(output_dir, paste0(prefix, "_11_CDF_Family.png")),
      width = 800, height = 600)
  df11 <- do.call(rbind, lapply(thetas, function(th) {
    data.frame(x = x_seq, y = parvind(x_seq, th),
               theta = paste0("theta = ", th))
  }))
  p11 <- ggplot2::ggplot(df11, ggplot2::aes(x = .data$x, y = .data$y,
                                              colour = .data$theta)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::labs(title = "Arvind Distribution CDF Family",
                  x = "x", y = "F(x)", colour = "") +
    theme_paper
  print(p11)
  grDevices::dev.off()

  # --- Plot 13: Residual Time Series ---
  grDevices::png(file.path(output_dir, paste0(prefix, "_13_Resid_TS.png")),
      width = 800, height = 600)
  df13 <- data.frame(Index = seq_len(n), Residuals = resid)
  p13 <- ggplot2::ggplot(df13, ggplot2::aes(x = .data$Index,
                                              y = .data$Residuals)) +
    ggplot2::geom_line(colour = "steelblue") +
    ggplot2::geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
    ggplot2::labs(title = paste(fit$model_type, "- Residual Time Series"),
                  x = "Time Index", y = "Residuals") +
    theme_paper
  print(p13)
  grDevices::dev.off()

  # --- Plot 14: ECDF vs Arvind CDF ---
  grDevices::png(file.path(output_dir, paste0(prefix, "_14_ECDF.png")),
      width = 800, height = 600)
  ecdf_vals <- stats::ecdf(e_pos)
  x_grid14 <- seq(min(e_pos), max(e_pos), length.out = 200)
  df14 <- data.frame(x = x_grid14,
                     ECDF = ecdf_vals(x_grid14),
                     Arvind = parvind(x_grid14, theta))
  p14 <- ggplot2::ggplot(df14, ggplot2::aes(x = .data$x)) +
    ggplot2::geom_step(ggplot2::aes(y = .data$ECDF, colour = "ECDF")) +
    ggplot2::geom_line(ggplot2::aes(y = .data$Arvind, colour = "Arvind CDF")) +
    ggplot2::labs(title = paste(fit$model_type, "- ECDF vs Arvind CDF"),
                  x = "x", y = "Probability", colour = "") +
    theme_paper
  print(p14)
  grDevices::dev.off()

  # --- Plot 19: P-P Plot ---
  grDevices::png(file.path(output_dir, paste0(prefix, "_19_PP_Plot.png")),
      width = 800, height = 600)
  e_sorted19 <- sort(e_pos)
  m19 <- length(e_sorted19)
  empirical_p <- (seq_len(m19) - 0.5) / m19
  theoretical_p <- parvind(e_sorted19, theta)
  df19 <- data.frame(Theoretical = theoretical_p, Empirical = empirical_p)
  p19 <- ggplot2::ggplot(df19, ggplot2::aes(x = .data$Theoretical,
                                              y = .data$Empirical)) +
    ggplot2::geom_point(colour = "steelblue") +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "red",
                         linetype = "dashed") +
    ggplot2::labs(title = paste(fit$model_type, "- P-P Plot"),
                  x = "Theoretical CDF", y = "Empirical CDF") +
    theme_paper
  print(p19)
  grDevices::dev.off()

  # --- Plot 21: Scale-Location ---
  grDevices::png(file.path(output_dir, paste0(prefix, "_21_Scale_Location.png")),
      width = 800, height = 600)
  df21 <- data.frame(Fitted = Y_hat,
                     SqrtAbsResid = sqrt(abs(resid)))
  p21 <- ggplot2::ggplot(df21, ggplot2::aes(x = .data$Fitted,
                                              y = .data$SqrtAbsResid)) +
    ggplot2::geom_point(colour = "steelblue") +
    ggplot2::geom_smooth(method = "loess", colour = "red", se = FALSE) +
    ggplot2::labs(title = paste(fit$model_type, "- Scale-Location"),
                  x = "Fitted Values", y = "sqrt(|Residuals|)") +
    theme_paper
  print(p21)
  grDevices::dev.off()

  message(sprintf("  %s: Diagnostic plots saved to %s",
              fit$model_type, output_dir))
  invisible(fit)
}
