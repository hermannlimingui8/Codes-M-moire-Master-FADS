# -----------------------------------------------------------------------------
# setup.R — installation
# -----------------------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  quantmod, ggplot2, tseries, psych, tidyr, stringr,
  TSA, forecast, FinTS, lmtest, moments,
  rugarch, zoo, dplyr, reshape2, tibble, PerformanceAnalytics
)
