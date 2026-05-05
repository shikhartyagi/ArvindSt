#' Summary and Comparison of Multiple Arvind Models
#'
#' Accepts multiple \code{ArvindFit} objects, computes diagnostics for each,
#' produces a unified comparison table, and prints the best model by RMSE,
#' R-squared, and AIC.
#'
#' @param ... one or more objects of class \code{"ArvindFit"}.
#' @param comparison_plots logical; if \code{TRUE} (default), generate
#'   comparison plots.
#' @param output_dir character; directory to save comparison plots.
#'   Defaults to a temporary directory.
#'
#' @return A data frame of diagnostic metrics (one row per model) is returned
#'   invisibly.
#'
#' @examples
#' dat <- simulate_arvind_data(n = 50, seed = 1)
#' m1 <- fit_rw1(Y ~ X1 + X2 + X3, dat, seed = 42)
#' summary_arvind(m1)
#'
#' @seealso [diagnostics_arvind()], [plot_arvind()]
#'
#' @export
summary_arvind <- function(...,
                           comparison_plots = TRUE,
                           output_dir = tempdir()) {
  fits <- list(...)
  comp_df <- do.call(rbind,
    lapply(fits, diagnostics_arvind))

  message("\n=== Model Comparison ===")
  print(comp_df, row.names = FALSE, digits = 4)
  message(sprintf("\n  Best RMSE: %s",
    comp_df$Model[which.min(comp_df$RMSE)]))
  message(sprintf("  Best R2:   %s",
    comp_df$Model[which.max(comp_df$R2)]))
  message(sprintf("  Best AIC:  %s",
    comp_df$Model[which.min(comp_df$AIC)]))

  if (comparison_plots && length(fits) > 1) {
    if (!dir.exists(output_dir))
      dir.create(output_dir, recursive = TRUE)

    theme_paper <- ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold", hjust = 0.5, size = 14),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "bottom")

    # C1: All Models Overlay
    grDevices::png(file.path(output_dir, "C1_All_Models_Overlay.png"),
        width = 1000, height = 600)
    n_max <- max(sapply(fits, function(f) f$n))
    df_all <- do.call(rbind, lapply(fits, function(f) {
      data.frame(
        Index = seq_len(f$n),
        Fitted = f$fitted,
        Model = f$model_type,
        stringsAsFactors = FALSE)
    }))
    df_obs <- data.frame(Index = seq_len(fits[[1]]$n),
                         Y = fits[[1]]$Y)
    p_c1 <- ggplot2::ggplot() +
      ggplot2::geom_line(data = df_obs,
                         ggplot2::aes(x = .data$Index, y = .data$Y),
                         colour = "black", linewidth = 1) +
      ggplot2::geom_line(data = df_all,
                         ggplot2::aes(x = .data$Index, y = .data$Fitted,
                                      colour = .data$Model),
                         linetype = "dashed") +
      ggplot2::labs(title = "All Models - Observed vs Fitted",
                    x = "Time Index", y = "Value", colour = "") +
      theme_paper
    print(p_c1)
    grDevices::dev.off()

    # C3: Error Metrics Bars
    grDevices::png(file.path(output_dir, "C3_Error_Metrics.png"),
        width = 800, height = 600)
    df_metrics <- data.frame(
      Model = comp_df$Model,
      RMSE = comp_df$RMSE,
      MAE = comp_df$MAE,
      stringsAsFactors = FALSE)
    df_long <- reshape2::melt(df_metrics, id.vars = "Model",
                              variable.name = "Metric", value.name = "Value")
    p_c3 <- ggplot2::ggplot(df_long,
                             ggplot2::aes(x = .data$Model, y = .data$Value,
                                          fill = .data$Metric)) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") +
      ggplot2::labs(title = "Error Metrics Comparison",
                    x = "", y = "Value", fill = "") +
      theme_paper
    print(p_c3)
    grDevices::dev.off()

    # C5: Residual Boxplots
    grDevices::png(file.path(output_dir, "C5_Residual_Boxplots.png"),
        width = 800, height = 600)
    df_resid <- do.call(rbind, lapply(fits, function(f) {
      data.frame(Model = f$model_type, Residuals = f$residuals,
                 stringsAsFactors = FALSE)
    }))
    p_c5 <- ggplot2::ggplot(df_resid,
                             ggplot2::aes(x = .data$Model,
                                          y = .data$Residuals,
                                          fill = .data$Model)) +
      ggplot2::geom_boxplot(alpha = 0.7, show.legend = FALSE) +
      ggplot2::labs(title = "Residual Boxplots", x = "", y = "Residuals") +
      theme_paper
    print(p_c5)
    grDevices::dev.off()

    # C8: R2 / Adj R2
    grDevices::png(file.path(output_dir, "C8_R2_Comparison.png"),
        width = 800, height = 600)
    df_r2 <- data.frame(
      Model = comp_df$Model,
      R2 = comp_df$R2,
      AdjR2 = comp_df$AdjR2,
      stringsAsFactors = FALSE)
    df_r2_long <- reshape2::melt(df_r2, id.vars = "Model",
                                  variable.name = "Metric",
                                  value.name = "Value")
    p_c8 <- ggplot2::ggplot(df_r2_long,
                             ggplot2::aes(x = .data$Model, y = .data$Value,
                                          fill = .data$Metric)) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") +
      ggplot2::labs(title = "R-squared Comparison",
                    x = "", y = "Value", fill = "") +
      theme_paper
    print(p_c8)
    grDevices::dev.off()

    message(sprintf("  Comparison plots saved to %s", output_dir))
  }

  invisible(comp_df)
}
