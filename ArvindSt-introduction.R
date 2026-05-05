## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----distribution-------------------------------------------------------------
library(ArvindSt)

# PDF at several points
darvind(c(0.5, 1, 2), theta = 1)

# CDF
parvind(1, theta = 2)

# Quantiles
qarvind(c(0.25, 0.5, 0.75), theta = 1)

# Random generation
set.seed(42)
x <- rarvind(1000, theta = 2)
summary(x)

## ----pdf-plot-----------------------------------------------------------------
x_seq <- seq(0.01, 4, length.out = 300)
thetas <- c(0.5, 1, 2, 5)

plot(NULL, xlim = c(0, 4), ylim = c(0, 1.5),
     xlab = "x", ylab = "f(x)", main = "Arvind PDF Family")
cols <- c("red", "blue", "darkgreen", "purple")
for (i in seq_along(thetas)) {
  lines(x_seq, darvind(x_seq, thetas[i]), col = cols[i], lwd = 2)
}
legend("topright", paste("theta =", thetas),
       col = cols, lwd = 2, cex = 0.8)

## ----mle----------------------------------------------------------------------
set.seed(42)
x <- rarvind(500, theta = 2)
fit <- fit_arvind_mle(x)
cat("Estimated theta:", fit$theta, "\n")
cat("True theta: 2\n")

## ----model-fit----------------------------------------------------------------
# Generate simulated data
dat <- simulate_arvind_data(n = 60, seed = 1)

# Fit RW1 model
m1 <- fit_rw1(Y ~ X1 + X2 + X3, dat, seed = 42)
cat("Model:", m1$model_type, "\n")
cat("Theta:", m1$theta, "\n")
cat("R-squared:", 1 - sum(m1$residuals^2) / sum((m1$Y - mean(m1$Y))^2), "\n")

## ----diagnostics--------------------------------------------------------------
d1 <- diagnostics_arvind(m1)
d1[, c("Model", "RMSE", "R2", "AIC", "KS_pvalue")]

## ----cv-----------------------------------------------------------------------
cv1 <- cv_arvind(m1, k_folds = 3, rolling = FALSE, seed = 42)
cat("Mean CV RMSE:", cv1$mean_cv_rmse, "\n")

