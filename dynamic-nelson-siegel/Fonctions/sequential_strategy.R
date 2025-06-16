###############################################################################
# FONCTION QUI PROCEDE A LA STRATEGIE SEQUENTIELLE POUR LE TEST ADF
###############################################################################

ADF_sequential <- function(xts_data) {
  library(urca)
  library(tseries)
  library(xts)
  
  # Initialisation des résultats pour la stratégie séquentielle
  ADF_sequential_stats <- data.frame(
    Variable = character(),
    ADF_Trend_Stat = numeric(),
    ADF_Trend_Critical_5pct = numeric(),
    Trend_TStat = numeric(),
    ADF_Drift_Stat = numeric(),
    ADF_Drift_Critical_5pct = numeric(),
    Drift_TStat = numeric(),
    ADF_None_Stat = numeric(),
    ADF_None_Critical_5pct = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Initialisation des résultats pour les p-values
  pvalues <- data.frame(
    Variable = character(),
    ADF_Stat = numeric(),
    Pvalue = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Boucle sur chaque colonne
  for (i in 1:ncol(xts_data)) {
    col_name <- colnames(xts_data)[i]
    ts_series <- as.numeric(xts_data[, i])
    ts_series <- ts_series[!is.na(ts_series)] 
    
    # Test avec tendance 
    ur_df_trend <- ur.df(ts_series, type = "trend", selectlags = "AIC")
    trend_stat <- ur_df_trend@teststat[1]
    trend_crit5 <- ur_df_trend@cval[1,2]
    trend_tstat <- ur_df_trend@testreg$coefficients[3,3]
    
    # Test avec drift 
    ur_df_drift <- ur.df(ts_series, type = "drift", selectlags = "AIC")
    drift_stat <- ur_df_drift@teststat[1]
    drift_crit5 <- ur_df_drift@cval[1,2]
    drift_tstat <- ur_df_drift@testreg$coefficients[1,3]
    
    # Test sans constante ni drift
    ur_df_none <- ur.df(ts_series, type = "none", selectlags = "AIC")
    none_stat <- ur_df_none@teststat[1]
    none_crit5 <- ur_df_none@cval[1,2]
    
    # Récupération des résultats pour la stratégie séquentielle
    ADF_sequential_stats <- rbind(ADF_sequential_stats, data.frame(
      Variable = col_name,
      ADF_Trend_Stat = trend_stat,
      ADF_Trend_Critical_5pct = trend_crit5,
      Trend_TStat = trend_tstat,
      ADF_Drift_Stat = drift_stat,
      ADF_Drift_Critical_5pct = drift_crit5,
      Drift_TStat = drift_tstat,
      ADF_None_Stat = none_stat,
      ADF_None_Critical_5pct = none_crit5
    ))
    
    adf_test <- adfTest(ts_series, type = "c") 
    adf_stat <- adf_test@test$statistic
    adf_pval <- adf_test@test$p.value
    
    pvalues <- rbind(pvalues, data.frame(
      Variable = col_name,
      ADF_Stat = adf_stat,
      Pvalue = adf_pval
    ))
  }
  
  list(
    ADF_sequential_stats = ADF_sequential_stats,
    Pvalues = pvalues
  )
}

###############################################################################
# FONCTION QUI PROCEDE A LA STRATEGIE SEQUENTIELLE POUR LE TEST KPSS
###############################################################################
KPSS_sequential <- function(xts_data) {
  library(urca)
  library(tseries)
  library(xts)
  
  # Initialisation des résultats pour KPSS
  KPSS_sequential_stats <- data.frame(
    Variable = character(),
    KPSS_Trend_Stat = numeric(),
    KPSS_Trend_Critical_5pct = numeric(),
    KPSS_Drift_Stat = numeric(),
    KPSS_Drift_Critical_5pct = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Initialisation des résultats pour les p-values
  pvalues <- data.frame(
    Variable = character(),
    KPSS_Stat = numeric(),
    Pvalue = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Boucle sur chaque colonne
  for (i in 1:ncol(xts_data)) {
    col_name <- colnames(xts_data)[i]
    ts_series <- as.numeric(xts_data[, i])
    ts_series <- ts_series[!is.na(ts_series)]
    
    # Test KPSS avec tendance
    kpss_trend <- ur.kpss(ts_series, type = "tau", lags = "short")
    trend_stat <- kpss_trend@teststat
    trend_crit5 <- kpss_trend@cval[2]  
    
    # Test KPSS avec drift
    kpss_level <- ur.kpss(ts_series, type = "mu", lags = "short")
    level_stat <- kpss_level@teststat
    level_crit5 <- kpss_level@cval[2]  # Valeur critique à 5%
    
    # Récupération des résultats
    KPSS_sequential_stats <- rbind(KPSS_sequential_stats, data.frame(
      Variable = col_name,
      KPSS_Trend_Stat = trend_stat,
      KPSS_Trend_Critical_5pct = trend_crit5,
      KPSS_Level_Stat = level_stat,
      KPSS_Level_Critical_5pct = level_crit5
    ))
    
    # Test KPSS avec p-value 
    kpss_test <- kpss.test(ts_series, null = "Level")
    kpss_stat <- kpss_test$statistic
    kpss_pval <- kpss_test$p.value
    
    pvalues <- rbind(pvalues, data.frame(
      Variable = col_name,
      KPSS_Stat = kpss_stat,
      Pvalue = kpss_pval
    ))
  }
  
  list(
    KPSS_sequential_stats = KPSS_sequential_stats,
    Pvalues = pvalues
  )
}