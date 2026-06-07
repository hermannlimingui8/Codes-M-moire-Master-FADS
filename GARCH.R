# -----------------------------------------------------------------------------
# GARCH.R : Modélisation de la volatilité conditionnelle (ARMA+GARCH)
# -----------------------------------------------------------------------------

arma_path <- file.path(tempdir(), "arma_results.RDS")
if (!file.exists(arma_path)) stop("⚠️ arma_results.RDS manquant. Exécutez ARMA.R d'abord.")
arma_obj <- readRDS(arma_path)

log_returns <- arma_obj$log_returns
p = 0 # <- arma_obj$p; 
q <- arma_obj$q; include_mean <- arma_obj$include.mean

# 1. Grille de modèles(réduits pour gain de temps-on peut étendre après)
grid <- expand.grid(
  vm = c("sGARCH", "gjrGARCH", "eGARCH"),      # modèles de variance
  dist = c("norm", "std", "sstd"),           # distributions     
  stringsAsFactors = FALSE
)

message("🎯 Estimation de ", nrow(grid), " modèles GARCH...")

# 2. Fonction robuste d'ajustement (retourne spec + fit)
fit_garch_model <- function(vm, dist, data, p, q, include.mean = TRUE) {
  spec <- tryCatch({
    ugarchspec(
      variance.model = list(model = vm, garchOrder = c(1, 1)),
      mean.model = list(armaOrder = c(p, q), include.mean = include.mean),
      distribution.model = dist
    )
  }, error = function(e) {
    warning("⚠️ Spécification échouée (", vm, " + ", dist, "): ", conditionMessage(e))
    return(NULL)
  })
  if (is.null(spec)) return(list(spec = NULL, fit = NULL))
  
  fit <- tryCatch({
    ugarchfit(spec, data = data, solver = "hybrid", fit.control = list(scale = 1))
  }, error = function(e) {
    warning("❌ Ajustement échoué (", vm, " + ", dist, "): ", conditionMessage(e))
    return(list(spec = spec, fit = NULL))
  })
  
  if (is.null(fit@fit$convergence) || fit@fit$convergence != 0) {
    warning("⚠️ Convergence douteuse (", vm, " + ", dist, "), code =", fit@fit$convergence)
  }
  
  return(list(spec = spec, fit = fit))
}

# 3. Boucle d'ajustement
models <- list()
pb <- txtProgressBar(min = 1, max = nrow(grid), style = 3)
for (i in seq_len(nrow(grid))) {
  vm <- as.character(grid$vm[i])
  dist <- as.character(grid$dist[i])
  key <- paste(vm, dist, sep = "_")
  
  res <- fit_garch_model(vm, dist, log_returns, p, q, include.mean = include_mean)
  spec <- res$spec; fit <- res$fit
  
  if (!is.null(fit) && !is.null(spec)) {
    ic <- infocriteria(fit)
    models[[key]] <- list(
      spec = spec,
      fit = fit,
      AIC = ic["Akaike", 1],
      BIC = ic["Bayes", 1]
    )
  }
  setTxtProgressBar(pb, i)
}
close(pb)

if (length(models) == 0) stop("Aucun modèle n'a convergé. Vérifiez les données ou élargissez les options.")

# 4. Comparaison par BIC
comp_df <- do.call(rbind, lapply(names(models), function(k) {
  data.frame(
    Modèle = k,
    AIC = models[[k]]$AIC,
    BIC = models[[k]]$BIC,
    stringsAsFactors = FALSE
  )
}))
comp_df <- comp_df[order(comp_df$BIC), ]
print(comp_df)

best_key <- comp_df$Modèle[1]
best_fit <- models[[best_key]]$fit
best_spec <- models[[best_key]]$spec  # ✅ Correction ici

cat("\n🏆 Meilleur modèle (BIC) :", best_key, "\n")
print(best_fit)

# 5. Diagnostics & Graphiques (Marges sécurisées)
old_par <- par(no.readonly = TRUE)  # Sauvegarde l'état graphique actuel
par(mfrow = c(2, 1), mar = c(3, 3, 2, 1) + 0.1)  # Marges réduites pour éviter l'erreur

plot(best_fit@model$modeldata$data, type = "l", col = "darkblue",
     main = paste("Rendements log —", arma_obj$ticker), ylab = "r_t", xlab = "")
plot(sigma(best_fit), type = "l", col = "red",
     main = paste("Volatilité conditionnelle —", best_key), ylab = expression(sigma[t]))
abline(h = mean(sigma(best_fit)), lty = 2, col = "gray")

par(mfrow = c(1, 1))  # Réinitialise la grille à 1 graphique

# 6. Sauvegarde complète & sécurisée
output_path <- file.path(tempdir(), "garch_results.RDS")

# Fallback robuste pour les dates si data_obj n'existe pas dans cet environnement
start_date <- if (exists("data_obj")) data_obj$start_date else index(log_returns)[1]
end_date   <- if (exists("data_obj")) data_obj$end_date   else index(log_returns)[length(log_returns)]

saveRDS(list(
  best_spec = best_spec,
  best_fit = best_fit,
  best_key = best_key,
  comparison_table = comp_df,
  log_returns = log_returns,
  arma_orders = c(p = p, q = q),
  include_mean = include_mean,
  model_info = list(
    ticker = arma_obj$ticker,
    start_date = start_date,
    end_date = end_date,
    date_run = Sys.time(),
    n_obs = length(log_returns)
  )
), file = output_path)

cat("💾 Résultats GARCH sauvegardés dans ", basename(output_path), "\n")