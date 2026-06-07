# ------------------------------------------------------------------------------
# Value-at-Risk.R : Estimation + Backtesting (eGARCH-sstd)
# ------------------------------------------------------------------------------
library(rugarch)
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Chargement des données
garch_path <- file.path(tempdir(), "garch_results.RDS")
if (!file.exists(garch_path)) stop("📦 garch_results.RDS manquant. Exécutez GARCH.R d'abord.")
garch_obj <- readRDS(garch_path)
log_returns <- garch_obj$log_returns

# 🔧 Extraction robuste (compatible xts/zoo/ts)
dates_all <- as.Date(index(log_returns))
r <- as.numeric(coredata(log_returns))

alpha_level <- 0.05
window_size <- 2000

message("📉 Estimation rolling GARCH (eGARCH-sstd)...")
spec_garch <- ugarchspec(
  variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "sstd"
)

roll_garch <- tryCatch({
  ugarchroll(spec_garch, data = log_returns, n.ahead = 1,
             n.start = window_size, refit.every = 25, refit.window = "moving",
             solver = "hybrid", calculate.VaR = TRUE, VaR.alpha = alpha_level)
}, error = function(e) stop("❌ Échec ugarchroll: ", conditionMessage(e)))

if ((convergence(roll_garch) > 0)) {
  message("🔁 Reprise des fenêtres non convergées...")
  roll_garch <- resume(roll_garch, solver = "gosolnp", solver.control = list(n.restarts = 10))
}

# Extraction VaR
df_garch <- as.data.frame(roll_garch, which = "VaR")
v_col <- grep(paste0("VaR\\.", alpha_level), names(df_garch), value = TRUE)
if (length(v_col) == 0) v_col <- grep("alpha\\(", names(df_garch), value = TRUE)[1]
if (length(v_col) == 0) stop("Impossible d'extraire la colonne VaR GARCH.")

VaR_garch <- as.numeric(df_garch[[v_col]])
realized <- as.numeric(df_garch$realized)
dates_garch <- as.Date(rownames(df_garch))

keep_idx <- is.finite(VaR_garch) & is.finite(realized)
df_garch_clean <- data.frame(date = dates_garch[keep_idx], 
                             VaR_GARCH = VaR_garch[keep_idx], 
                             realized = realized[keep_idx])

# 2. Méthodes alternatives
message("🔄 Calcul des VaR alternatives (rolling)...")
compute_rolling_var <- function(r, window, alpha) {
  n <- length(r)
  VaR_hist <- VaR_norm <- VaR_rm <- rep(NA, n)
  sigma_rm_sq <- var(r[1:window], na.rm = TRUE)
  lambda <- 0.94
  
  for (i in (window + 1):n) {
    win <- r[(i - window):(i - 1)]
    VaR_hist[i] <- quantile(win, probs = alpha, type = 8)
    mu <- mean(win); sig <- sd(win)
    VaR_norm[i] <- mu + sig * qnorm(alpha)
    
    if (i > window + 1) {
      sigma_rm_sq <- lambda * sigma_rm_sq + (1 - lambda) * r[i - 1]^2
    }
    VaR_rm[i] <- sqrt(sigma_rm_sq) * qnorm(alpha)
  }
  list(VaR_hist = VaR_hist, VaR_norm = VaR_norm, VaR_rm = VaR_rm)
}

roll_vars <- compute_rolling_var(r, window_size, alpha_level)

# 🔧 Alignement sécurisé des dates
idx_match <- match(df_garch_clean$date, dates_all)

df_comp <- data.frame(
  date = df_garch_clean$date,
  realized = df_garch_clean$realized,
  VaR_GARCH = df_garch_clean$VaR_GARCH,
  VaR_Hist = roll_vars$VaR_hist[idx_match],
  VaR_Norm = roll_vars$VaR_norm[idx_match],
  VaR_RM = roll_vars$VaR_rm[idx_match]
) %>% filter(complete.cases(.))

message("📊 Comparaison finale sur ", nrow(df_comp), " jours.")

# 1. Extraction des résultats rolling
df_g <- as.data.frame(roll_garch, which = "VaR")
df_g$date <- as.Date(rownames(df_g))
rownames(df_g) <- NULL

# 🔍 Diagnostic : affiche les colonnes réelles
cat("📊 Colonnes disponibles :", paste(names(df_g), collapse = " | "), "\n")

# 🔧 Détection intelligente de la colonne VaR
v_col <- names(df_g)[grep("VaR|alpha", names(df_g), ignore.case = TRUE)]
if (length(v_col) == 0) {
  # Fallback : prend la colonne qui n'est ni realized, ni Mu, ni Sigma, ni date
  v_col <- setdiff(names(df_g), c("realized", "Mu", "Sigma", "date"))[1]
}
if (is.na(v_col) || length(v_col) == 0) {
  stop("❌ Colonne VaR introuvable. Vérifiez que calculate.VaR = TRUE dans ugarchroll.")
}
v_col <- v_col[1]
cat("✅ Colonne VaR détectée :", v_col, "\n")

# 2. Préparer les VaR alternatives alignées sur log_returns
df_alt <- tibble::tibble(
  date     = as.Date(index(log_returns)),
  VaR_Hist = roll_vars$VaR_hist,
  VaR_Norm = roll_vars$VaR_norm,
  VaR_RM   = roll_vars$VaR_rm
)

# 3. Jointure relationnelle sécurisée
df_comp <- df_g %>%
  select(date, realized, garch_var = all_of(v_col)) %>%
  dplyr::left_join(df_alt, by = "date") %>%
  dplyr::rename(VaR_GARCH = garch_var) %>%
  dplyr::mutate(
    dplyr::across(c(realized, starts_with("VaR_")), ~ as.numeric(as.character(.x)))
  ) %>%
  dplyr::filter(dplyr::if_all(c(realized, starts_with("VaR_")), is.finite))

cat("📏 Lignes avant filtrage  :", nrow(df_g), "\n")
cat("✅ Lignes finales valides  :", nrow(df_comp), "\n")

if (nrow(df_comp) < 50) {
  stop("❌ df_comp vide ou trop petit. Vérifiez window_size et la qualité des données.")
}

# ==============================================================================
# 🔁 BACKTEST ROBUSTE - VERSION CORRIGÉE
# ==============================================================================

backtest_one <- function(name, VaR_vec, realized, alpha) {
  # Filtrage des valeurs finies
  keep <- is.finite(VaR_vec) & is.finite(realized)
  VaR_vec <- VaR_vec[keep]
  realized <- realized[keep]
  
  n <- length(realized)
  
  # Debug
  cat(sprintf("🔍 %s: n = %d, violations = %d\n", 
              name, n, sum(realized < VaR_vec, na.rm = TRUE)))
  
  if (n < 50) {
    warning(sprintf("⚠️  %s: trop peu d'observations (n=%d)", name, n))
    return(data.frame(
      Méthode = name, Violations = 0, Taux_viol = NA,
      p_UC = NA, p_CC = NA, stringsAsFactors = FALSE
    ))
  }
  
  # Indicateur de violation
  I <- as.integer(realized < VaR_vec)
  x <- sum(I)
  
  if (x == 0) {
    warning(sprintf("⚠️  %s: aucune violation", name))
    return(data.frame(
      Méthode = name, Violations = 0, Taux_viol = 0,
      p_UC = NA, p_CC = NA, stringsAsFactors = FALSE
    ))
  }
  
  p_hat <- x / n
  p_hat <- min(max(p_hat, 1e-9), 1 - 1e-9)
  
  # Fonction safe_xlog
  safe_xlog <- function(x_val, p_val) {
    if (x_val == 0) return(0.0)
    return(x_val * log(max(p_val, 1e-15)))
  }
  
  # Test de Kupiec
  L0_uc <- safe_xlog(n - x, 1 - alpha) + safe_xlog(x, alpha)
  L1_uc <- safe_xlog(n - x, 1 - p_hat) + safe_xlog(x, p_hat)
  LR_uc <- max(-2 * (L0_uc - L1_uc), 0.0)
  p_uc <- 1 - pchisq(LR_uc, df = 1)
  
  # Test de Christoffersen
  if (length(I) < 2) {
    p_cc <- NA
  } else {
    n00 <- sum(I[-length(I)] == 0 & I[-1] == 0)
    n01 <- sum(I[-length(I)] == 0 & I[-1] == 1)
    n10 <- sum(I[-length(I)] == 1 & I[-1] == 0)
    n11 <- sum(I[-length(I)] == 1 & I[-1] == 1)
    
    if ((n01 + n11) < 2) {
      LR_ind <- 0.0
    } else {
      pi1 <- if ((n00 + n01) > 0) n01 / (n00 + n01) else 1e-9
      pi2 <- if ((n10 + n11) > 0) n11 / (n10 + n11) else 1e-9
      pi <- (n01 + n11) / n
      
      pi1 <- min(max(pi1, 1e-9), 1 - 1e-9)
      pi2 <- min(max(pi2, 1e-9), 1 - 1e-9)
      pi <- min(max(pi, 1e-9), 1 - 1e-9)
      
      L0_ind <- safe_xlog(n00 + n10, 1 - pi) + safe_xlog(n01 + n11, pi)
      L1_ind <- (safe_xlog(n00, 1 - pi1) + safe_xlog(n01, pi1) +
                   safe_xlog(n10, 1 - pi2) + safe_xlog(n11, pi2))
      LR_ind <- max(-2 * (L0_ind - L1_ind), 0.0)
    }
    
    LR_cc <- LR_uc + LR_ind
    p_cc <- 1 - pchisq(LR_cc, df = 2)
  }
  
  return(data.frame(
    Méthode = name,
    Violations = x,
    Taux_viol = round(x / n, 4),
    p_UC = round(p_uc, 4),
    p_CC = round(p_cc, 4),
    stringsAsFactors = FALSE
  ))
}

# Appel des tests
cat("\n📊 Lancement du backtesting...\n")
back_df <- bind_rows(
  backtest_one("GARCH", df_comp$VaR_GARCH, df_comp$realized, alpha_level),
  backtest_one("Historique", df_comp$VaR_Hist, df_comp$realized, alpha_level),
  backtest_one("Normale", df_comp$VaR_Norm, df_comp$realized, alpha_level),
  backtest_one("RiskMetrics", df_comp$VaR_RM, df_comp$realized, alpha_level)
)

print(back_df, digits = 4)
# ==============================================================================
# ALIGNEMENT SÉCURISÉ + CALCUL DES PERTES (Version Blindée)
# ==============================================================================

# 1. Alignement strict des dates (élimine les NA silencieux de match())
message("🔗 Alignement des dates...")
dates_ref <- as.Date(rownames(df_garch))
dates_all_strict <- as.Date(dates_all)
df_garch_clean$date <- as.Date(df_garch_clean$date)

idx_match <- match(df_garch_clean$date, dates_all_strict)
cat("📊 Dates alignées :", sum(!is.na(idx_match)), "/", nrow(df_garch_clean), "\n")

if (all(is.na(idx_match))) stop("❌ Aucun match de dates. Vérifiez que log_returns et df_garch couvrent la même période.")

# 2. Construction du dataframe comparatif
df_comp <- tibble(
  date = df_garch_clean$date,
  realized = df_garch_clean$realized,
  VaR_GARCH = df_garch_clean$VaR_GARCH,
  VaR_Hist = roll_vars$VaR_hist[idx_match],
  VaR_Norm = roll_vars$VaR_norm[idx_match],
  VaR_RM = roll_vars$VaR_rm[idx_match]
)

cat("📈 VaR_Hist sur la plage alignée :", 
    sum(is.finite(roll_vars$VaR_hist[idx_match])), "valeurs valides\n")

# 3. Conversion explicite (tue le warning de coercion)
df_comp <- df_comp %>%
  mutate(
    realized = as.numeric(as.character(realized)),
    across(starts_with("VaR_"), ~ as.numeric(as.character(.x)))
  )

# 4. Diagnostic étape par étape
cat("📏 Lignes avant filtrage :", nrow(df_comp), "\n")
cat("📉 NAs par colonne :", paste(names(colSums(is.na(df_comp))), 
                                  colSums(is.na(df_comp)), collapse = ", "), "\n")

# Filtrage strict (NA, NaN, Inf)
df_clean <- df_comp %>%
  filter(if_all(c(realized, starts_with("VaR_")), is.finite))

cat("✅ Lignes valides pour calcul :", nrow(df_clean), "\n")
if (nrow(df_clean) == 0) {
  stop("❌ STOP: 0 ligne valide. Vérifiez que window_size <= length(log_returns) et que ugarchroll a convergé.")
}

# 5. Calcul des pertes (Quantile Loss)
alpha <- alpha_level
df_clean <- df_clean %>%
  mutate(
    Loss_GARCH = (realized - VaR_GARCH) * (alpha - as.integer(realized < VaR_GARCH)),
    Loss_Hist  = (realized - VaR_Hist)  * (alpha - as.integer(realized < VaR_Hist)),
    Loss_Norm  = (realized - VaR_Norm)  * (alpha - as.integer(realized < VaR_Norm)),
    Loss_RM    = (realized - VaR_RM)    * (alpha - as.integer(realized < VaR_RM))
  )

# 6. Résumé (syntaxe dplyr 1.1.0+ compatible)
loss_summary <- df_clean %>%
  select(starts_with("Loss_")) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Méthode", values_to = "Moyenne_perte") %>%
  mutate(Méthode = str_remove(Méthode, "Loss_")) %>%
  arrange(Moyenne_perte)

print(loss_summary)

# 6. Sauvegarde
final_save <- file.path(getwd(), paste0("VaR_Results_", gsub("[^A-Za-z0-9]", "_", garch_obj$model_info$ticker), "_", Sys.Date(), ".RDS"))
saveRDS(list(
  data = df_comp,
  backtest = back_df,
  losses = loss_summary,
  garch_model = garch_obj$best_fit,
  metadata = list(
    alpha = alpha_level,
    window_size = window_size,
    model_used = "eGARCH-sstd",
    date_run = Sys.time(),
    ticker = garch_obj$model_info$ticker
  )
), file = final_save)

cat("✅ Résultats sauvegardés dans :", final_save, "\n")