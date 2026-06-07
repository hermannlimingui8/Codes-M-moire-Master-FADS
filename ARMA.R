# -----------------------------------------------------------------------------
# ARMA.R : Identification + estimation + validation ARMA (Box-Jenkins)
# -----------------------------------------------------------------------------

input_path <- file.path(tempdir(), "data_preprocessed.RDS")
if (!file.exists(input_path)) stop("⚠️ data_preprocessed.RDS manquant. Exécutez Data.R d'abord.")
data_obj <- readRDS(input_path)
log_returns <- data_obj$log_returns

message("🚀 Modélisation ARMA sur ", nrow(log_returns), " observations.")

# 1. Identification visuelle
par(mfrow = c(1, 2))
acf(log_returns, lag.max = 10, main = "ACF")
pacf(log_returns, lag.max = 10, main = "PACF")
par(mfrow = c(1, 1))

# eacf()
message("🔍 eacf() (AR jusqu'à 5, MA jusqu'à 5)...")
eacf(log_returns, ar.max = 5, ma.max = 5)

# 2. Sélection automatique
auto_models <- list(
  AIC  = auto.arima(log_returns, ic = "aic",  max.p = 5, max.q = 5, stepwise = FALSE),
  AICc = auto.arima(log_returns, ic = "aicc", max.p = 5, max.q = 5, stepwise = FALSE),
  BIC  = auto.arima(log_returns, ic = "bic",  max.p = 5, max.q = 5, stepwise = FALSE)
)

# comparaison des ordres
orders <- sapply(auto_models, function(m) paste(m$arma[1], m$arma[2], sep = ","))
cat("Ordres proposés :\n")
print(data.frame(Critère = names(orders), ARMA = orders), row.names = FALSE)

# Sélection finale BIC
best_auto <- auto_models[["BIC"]]
p <- best_auto$arma[1]; q <- best_auto$arma[2]
include_mean <- TRUE  # auto.arima() inclut la moyenne par défaut si significative

cat("\n✅ Modèle retenu (BIC) : ARMA(", p, ",", q, "), include.mean =", include_mean, "\n")

# 3. Estimation
fit <- Arima(log_returns, order = c(p, 0, q), include.mean = include_mean)
summary(fit)

# 4. Validation
message("\n🧪 Tests de validation...")
checkresiduals(fit) # graphiques + tests

resid <- residuals(fit)

# Tests sur les résidus
# Test de Jarque-Bera
cat("Jarque-Bera (H0: normalité) : p-value =", round(jarque.bera.test(resid)$p.value, 4), "\n")

# Test de Durbin-Watson
cat("Durbin-Watson (H0: pas d'autocorrélation)  : DW =", round(dwtest(lm(resid ~ 1))$statistic, 4), "\n")

# Test de Ljung-Box
for (lags in 1:5) {
  lbtest <- Box.test(resid, lag = lags, type = "Ljung-Box")
  cat("Ljung-Box(", lags, ") (H0: indépendance) : p-value =", round(lbtest$p.value, 4), "\n")
}

# Test de Box-Pierce
for (lags in 1:5) {
  lbtest <- Box.test(resid, lag = lags, type = "Box-Pierce")
  cat("Box-Pierce(", lags, ") (H0: indépendance) : p-value =", round(lbtest$p.value, 4), "\n")
}

# Test ARCH
for (lags in 1:10) {
  cat("ARCH(", lags, ") (H0: homoscédasticité) : p-value =", round(ArchTest(resid) $p.value, 4), "\n")
}

# 5. Sauvegarde — ajout explicite de include.mean
output_path <- file.path(tempdir(), "arma_results.RDS")
saveRDS(list(
  fit = fit,
  p = p, q = q,
  include.mean = include_mean,
  residuals = resid,
  log_returns = log_returns,
  auto_models = auto_models,
  ticker = data_obj$ticker
), file = output_path)

cat("💾 Résultats ARMA sauvegardés dans ", basename(output_path), "\n")