# ArvindSt

<!-- badges: start -->
<!-- badges: end -->

**ArvindSt** is an R package implementing the Arvind distribution and five novel stochastic regression models with Arvind-distributed errors.

## Installation

Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("shikhartyagi/ArvindSt")
```

## Features

- **Distribution Functions**: `darvind()`, `parvind()`, `qarvind()`, `rarvind()`, `rarvind_centred()`
- **Five Regression Models**:
  - `fit_rw1()` — Random Walk on Coefficients
  - `fit_tvlm()` — Time-Varying Coefficient Linear Model
  - `fit_simex()` — Simulation-Extrapolation
  - `fit_mixed()` — Mixed-Effects Regression
  - `fit_hmm()` — Regime-Switching (HMM)
- **Diagnostics**: `diagnostics_arvind()` (21 metrics), `plot_arvind()` (25 plots)
- **Forecasting**: `forecast_arvind()` with Monte Carlo prediction intervals
- **Cross-Validation**: `cv_arvind()` (k-fold and rolling-window)
- **Model Comparison**: `summary_arvind()`

## Quick Start

```r
library(ArvindSt)

# Load example data
data(climate_consumption)

# Define formula
frm <- Consumption ~ Precip + TempMaxAvg + TempMinAvg + HumidMax + HumidAvg

# Fit all five models
m1 <- fit_rw1(frm, climate_consumption)
m2 <- fit_tvlm(frm, climate_consumption)
m3 <- fit_simex(frm, climate_consumption, me_vars = c("Precip", "TempMaxAvg"))
m4 <- fit_mixed(frm, climate_consumption, group_var = "Season")
m5 <- fit_hmm(frm, climate_consumption, nstates = 2)

# Compare all models
summary_arvind(m1, m2, m3, m4, m5)
```

## Authors

- **Shikhar Tyagi** (maintainer) — shikhar1093tyagi@gmail.com
- **Arvind Pandey** — arvindmzu@gmail.com

## License

MIT
