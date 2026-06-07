# -----------------------------------------------------------------------------
# Data.R : Téléchargement + EDA
# -----------------------------------------------------------------------------

# 1. Paramètres
ticker <- "SPY"  # SPY"
start_date <- "2000-01-01"
end_date <- "2024-12-31"

message("📥 Téléchargement des données Yahoo Finance pour ", ticker, "...")

data <- tryCatch(
  getSymbols(ticker, src = "yahoo", from = start_date, to = end_date, auto.assign = FALSE),
  error = function(e) {
    stop("❌ Échec du téléchargement : ", e$message, call. = FALSE)
  }
)

print(data)

prices <- na.omit(Ad(data))
log_returns <- na.omit(diff(log(prices)))
colnames(log_returns) <- "r"

message("✅ Données récupérées : ", nrow(log_returns), " rendements quotidiens.")

# 2. Tests de stationnarité
stationarity_tests <- function(x, name = "Série") {
  message("\n🔍 Tests de stationnarité pour ", name, ":")
  cat("ADF (H0: non-stationnarité): p-value =", round(adf.test(x)$p.value, 4), "\n")
  cat("KPSS (H0: stationnarité): p-value =", round(kpss.test(x)$p.value, 4), "\n")
  cat("PP  (H0: non-stationnarité): p-value =", round(pp.test(x)$p.value, 4), "\n")
}

stationarity_tests(prices, "prix")
stationarity_tests(log_returns, "log-rendements")

# 3. DEA

par(mar = c(5, 4, 4, 2) + 0.1)

p1 <- hist(log_returns, prob = TRUE,
           main = paste("Distribution des log-rendements (", ticker, ")"),
           xlab = "Rendement", col = "lightblue", border = "white")
lines(density(log_returns), col = "red", lwd = 2)

boxplot(log_returns, main = "Boxplot des log-rendements", 
        ylab = "Rendement", col = "skyblue", horizontal = TRUE)

acf(log_returns, lag.max = 10, main = "ACF des log-rendements")
pacf(log_returns, lag.max = 10, main = "PACF des log-rendements")

desc_stats <- psych::describe(log_returns)
print(desc_stats)

# Test d'asymétrie d'Agostino(H0:Skewness=0)
skew_test <- agostino.test(log_returns)
print(skew_test)

# Test d'asymétrie d'Aplatissement d'Anacombe-Glynn(H0:Kustosis=3)
kurt_test <- anscombe.test(log_returns)
print(kurt_test)

# Normalité
jb_test <- jarque.bera.test(log_returns)

cat("\nTest de skewness (Agostino, H0: 0)   : p =", round(skew_test$p.value, 4), "\n")
cat("Test de kurtosis (Anscombe, H0: 3)  : p =", round(kurt_test$p.value, 4), "\n")
cat("Test Jarque-Bera (H0: normalité)    : p =", round(jb_test$p.value, 4), "\n")

# 4. Sauvegarde
output_path <- file.path(tempdir(), "data_preprocessed.RDS")
saveRDS(list(
  log_returns = log_returns,
  ticker = ticker,
  start_date = start_date,
  end_date = end_date,
  desc_stats = desc_stats,
  jb_pvalue = jb_test$p.value,
  skew_test = skew_test,
  kurt_test = kurt_test
), file = output_path)

message("💾 Données sauvegardées dans ", output_path)