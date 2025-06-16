###############################################################################
#MAIN ESTIMATION DYNAMIC NELSON SIEGEL
###############################################################################

###############################################################################
# 1. Import des bibliothèques, nettoyage et définition du répertoire
###############################################################################

# Import des bibliothèques
{
  if(!require(rstudioapi)) install.packages("rstudioapi"); library(rstudioapi)
  if(!require(readxl)) install.packages('readxl'); library(readxl)
  if(!require(plotly)) install.packages('plotly'); library(plotly)
  if(!require(dplyr)) install.packages('dplyr'); library(dplyr)
  if(!require(lubridate)) install.packages('lubridate'); library(lubridate)
  if(!require(RColorBrewer)) install.packages('RColorBrewer'); library(RColorBrewer)
  if(!require(xts)) install.packages('xts'); library(xts)
  if(!require(zoo)) install.packages("zoo"); library(zoo)
  if(!require(forecast)) install.packages("forecast"); library(forecast)
  if(!require(lmtest)) install.packages("lmtest"); library(lmtest)
  if(!require(tseries)) install.packages("tseries"); library(tseries)
  if(!require(FinTS)) install.packages("FinTS"); library(FinTS)
  if(!require(urca)) install.packages("urca"); library(urca)
  if(!require(vars)) install.packages("vars"); library(vars)
  if(!require(fUnitRoots)) install.packages("fUnitRoots"); library(fUnitRoots)
  if(!require(Metrics)) install.packages("Metrics"); library(Metrics)
  if(!require(mvtnorm)) install.packages("mvtnorm"); library(mvtnorm)
  if(!require(tidyr)) install.packages("tidyr"); library(tidyr)
  if(!require(ggplot2)) install.packages("ggplot2"); library(ggplot2)
  if(!require(writexl)) install.packages("writexl"); library(writexl)
}


# Nettoyage et répertoire
set.seed(100)
rm(list=ls())
graphics.off()
setwd(dir=dirname(rstudioapi::getSourceEditorContext()$path))

# Récupération des fonctions
{
  source("Fonctions/utils.R")
  source("Fonctions/classic_NS.R")
  source("Fonctions/DNS_Diebold_Li.R")
  source("Fonctions/sequential_strategy.R")
  source("Fonctions/TVP_VAR.R")
  source("Fonctions/curve_shocks.R")
}

###############################################################################
# 2. Récupération et affichage de la structure par terme des taux d'intérêt
###############################################################################
{
  daily_data <- as.data.frame(read_excel("Data/Structure taux.xlsx"))
  daily_data[, 2:ncol(daily_data)] <- daily_data[, 2:ncol(daily_data)] * 100
  dates <- as.Date(daily_data$Dates)
  daily_data[, -1] <- lapply(daily_data[, -1], na.approx)
  maturities <- c(1/12, 3/12, 1, 2, 5, 10, 15, 20)
  names(maturities) <- colnames(daily_data)[-1]
  surface <- plot_surface(data = daily_data,
                          maturities = maturities,
                          title = "Structure par terme des taux d'intérêt - 2010-2025") 

  # Conversion en xts pour format série temporelle
  daily_data_xts <- xts(daily_data[,-1], order.by = daily_data$Dates)
  
  # Passage en mensuel
  data_xts <- apply.monthly(daily_data_xts, last)
  data <- data.frame(Date = index(data_xts), coredata(data_xts))
}

x11()
plot.xts(data_xts[, c("EUR3M", "SWAP5Y", "SWAP10Y")],
         main = "Évolution de taux caractéristiques sur la période",
         multi.panel = FALSE,
         col = c("blue", "red", "black"),
         lwd = 2,
         ylab = "Taux (%)",
         legend.loc = "bottomleft")
surface

###############################################################################
# 3. Estimation des facteurs du modèle de Nelson Siegel standard
###############################################################################
{
  #On crée une matrice pour stocker les facteurs
  ns_factors_matrix <- matrix(NA, nrow = nrow(data_xts), ncol = 3)
  colnames(ns_factors_matrix) <- c("Level", "Slope", "Curvature")

  # et pour pour stocker les R²
  r_squared_vec <- numeric(nrow(data_xts))

  # Estimation des facteurs à chaque date
  for(i in 1:nrow(data_xts)) {
    yields_i <- as.numeric(data_xts[i, ])
    ns_est <- estimate_NS_single_date(yields_i, maturities)
    ns_factors_matrix[i, ] <- ns_est$beta
    r_squared_vec[i] <- ns_est$r_squared
  }
  
  # Conversion en objet xts
  ns_factors <- xts(ns_factors_matrix, order.by = index(data_xts))
  
  cat("Estimation NS classique terminée \n")
  cat("R² moyen:", round(mean(r_squared_vec, na.rm = TRUE), 4), "\n")

  # Graphique des 3 facteurs NS
  x11()
  plot.xts(ns_factors,
         main = "Facteurs Nelson-Siegel estimés",
         multi.panel = TRUE,
         yaxis.same = FALSE,
         col = c("darkblue", "darkred", "darkgreen"))


  # Graphiques de fitting 
  x11()
  plot_NS_fit_avg(data_xts, ns_factors, maturities,
                col_obs = "red", col_fit = "blue",
                legend_pos = "topleft")
  
  dates_random <- sort(sample(1:nrow(data_xts), 4))
  x11()
  plot_NS_fit_dates(data_xts, ns_factors, maturities, dates_random,
                  col_obs = "red", col_fit = "blue",
                  legend_pos = "topleft")
}

###############################################################################
# 4. DNS selon Diebold & Li (2006) -> AR(1) et VAR
###############################################################################

# AR ou ARMA sur les facteurs
{
  optimal_orders_level <- estimate_optimal_arma(ns_factors$Level, max_p = 5, max_q = 5,
                                    information_criterion = "bic")
  optimal_orders_slope <- estimate_optimal_arma(ns_factors$Slope, max_p = 5, max_q = 5,
                                          information_criterion = "bic")
  optimal_orders_curvature <- estimate_optimal_arma(ns_factors$Curvature, max_p = 5, max_q = 5,
                                          information_criterion = "bic")
  test_res_AR_level <- test_residuals(optimal_orders_level$best_model$residuals, 
                                      'Level')
  test_res_AR_slope <- test_residuals(optimal_orders_slope$best_model$residuals, 
                                      'Slope')
  test_res_AR_curvature <- test_residuals(optimal_orders_curvature$best_model$residuals, 
                                          'Curvature')
}

# VAR
{
  VAR_model <- estimate_optimal_var(ns_factors_matrix, max_p = 1,
                                    information_criterion = "bic")
  summary(VAR_model$best_model)
}

###############################################################################
# 5. COINTEGRATION
###############################################################################
{
  # Tests de stationnarité 
  results_ADF_factors <- ADF_sequential(ns_factors)
  results_KPSS_factors <- KPSS_sequential(ns_factors)

  # Tests de cointégration
  factors_ts <- ts(ns_factors[, c("Level", "Slope", "Curvature")], start = c(2010, 1), frequency = 12)
  VARselect(factors_ts, lag.max = 10, type = "const")$selection
  trace_johansen_test <- ca.jo(factors_ts,
                             type = "trace",    
                             ecdet = "const",   
                             K = 2)  
  eigen_johansen_test <- ca.jo(factors_ts,
                             type = "eigen",    
                             ecdet = "const",   
                             K = 2)  
  cat("=========TEST DE LA TRACE=========")
  summary(trace_johansen_test)
  cat("=========TEST DE LA VALEUR PROPRE=========")
  summary(eigen_johansen_test)
}

###############################################################################
# 6. PERFORMANCES PREDICTIVES DES MODELES DYNAMIQUES SIMPLES
###############################################################################
{
  # On considère un modèle estimé sur 2010 - 2022
  train <- window(factors_ts, end = c(2021,12))
  test <- window(factors_ts, start = c(2022, 1))

  # Estimation sur la période d'entraînement
  ar_level_train <- Arima(train[,"Level"], order = c(1,0,0))
  ar_slope_train <- Arima(train[,"Slope"], order = c(1,0,0))
  ar_curvature_train <- Arima(train[,"Curvature"], order = c(1,0,0))
  var_train <- VAR(train, p = 1, type = "const")
  johansen_train <- ca.jo(train, K = 2, type = "trace", ecdet = "const")
  vecm_train <- cajorls(johansen_train, r = 1)
  vecm_train_var <- vec2var(johansen_train, r = 1)

  # Prévisions
  ar_pred <- ts.union(
    forecast(ar_level_train, h = length(test))$mean,
    forecast(ar_slope_train, h = length(test))$mean,
    forecast(ar_curvature_train, h = length(test))$mean
  )

  var_pred <- predict(var_train, n.ahead = length(test))
  var_pred_df <- ts(sapply(var_pred$fcst, function(x) x[, 1]), start = c(2022, 1), frequency = 12)
  vecm_pred <- predict(vecm_train_var, n.ahead = length(test))
  vecm_pred_df <- ts(sapply(vecm_pred$fcst, function(x) x[, 1]), start = c(2022, 1), frequency = 12)

  rmse_ar <- colMeans((test - ar_pred)^2, na.rm = TRUE)^0.5
  rmse_var <- colMeans((test - var_pred_df)^2, na.rm = TRUE)^0.5
  rmse_vecm <- colMeans((test - vecm_pred_df)^2, na.rm = TRUE)^0.5

  results_prev <- rolling_rmse_dm_evaluation(ns_factors, 2)
  print(results_prev$dm)
}

###############################################################################
# 7. TVP-VAR DMS / DMA
###############################################################################
# Etape 0 : Tests de rupture et obtention des dates de rupture
{
  break_idx <- breakpoints(Level ~ Slope + Curvature, data = as.data.frame(ns_factors))
  break_indices <- breakpoints(break_idx)$breakpoints
  break_dates <- index(ns_factors)[break_indices]
  print(break_dates)
}

# Etape 1 : estimation 
{
  k <- ncol(ns_factors_matrix)
  d <- k^2
  lambda_grid <- c(1, 0.99, 0.985, 0.98, 0.97, 0.975, 0.95)
  X0 <- rep(0, d)
  P0 <- diag(10, d)
  H0 <- diag(0.1, k)
  res_dms_dma <- run_dms_dma(Y = ns_factors_matrix,
                           lambda_grid = lambda_grid,
                           kappa = 0.98,
                           alpha = 0.99,
                           X0 = X0,
                           P0 = P0,
                           H0 = H0)
}

# Etape 2 : comparaisons
{
  pi_df <- as.data.frame(res_dms_dma$pi)
  colnames(pi_df) <- paste0("lambda_", res_dms_dma$lambdas)
  pi_df$time <- seq(as.Date("2010-01-01"), by = "month", length.out = nrow(res_dms_dma$pi))
  
  pi_long <- pivot_longer(pi_df, -time, names_to = "lambda", values_to = "pi")

  ggplot(pi_long, aes(x = time, y = pi, color = lambda)) +
    geom_line() +
    theme_minimal() +
    labs(title = "Évolution des probabilités DMS/DMA", y = "Poids", x = "Temps")

  matplot(cbind(res_dms_dma$pred_dms[,1], res_dms_dma$pred_dma[,2]),
        type = "l", col = c("black", "blue", "red"), lty = 1,
        ylab = "Level", xlab = "Temps", main = "Prévision : DMS vs DMA vs Réel")
  legend("topright", legend = c("Réel", "DMS", "DMA"), col = c("black", "blue", "red"), lty = 1)

  matplot(cbind(ns_factors_matrix[,1], res_dms_dma$pred_dma[,1]),
        type = "l", col = c("black", "red"), lty = 1,
        ylab = "Level", xlab = "Temps", main = "Prévision : DMS vs DMA vs Réel")

  # Performances prédictives
  rmse_dma_dms <- evaluate_dms_dma_performance(res_dms_dma, ns_factors_matrix)
  rmse_fixed <- results_prev$rmse
  cat("====RMSE SSM==== \n")
  print(rmse_dma_dms)
  cat("====RMSE fixed====\n")
  print(rmse_fixed)
}

# Etape 3 : forecast
{
  dma_baseline <- forecast_path(res_dms_dma$X_list_all, res_dms_dma$pi, 
                           ns_factors_matrix[nrow(ns_factors_matrix),], h = 20, type = "dma")
  baseline_scenario <- reconstruct_yield_curve(dma_baseline, maturities, lambda = 0.0609)
  colnames(baseline_scenario) <- names(maturities)
  future_dates <- seq(from = data$Date[nrow(data)] + months(1), 
                     by = "month", 
                     length.out = nrow(baseline_scenario))
  yield_proj_for_plot <- data.frame(
    Dates = future_dates, 
    baseline_scenario             
  )
  
  plot_surface(yield_proj_for_plot, maturities, 
               title = "Projection de la structure par terme avec TVP-VAR-DMA")
}

# On affiche les coefficients d'un modèle TVP pour observer si oui ou non ils varient réellement
x11();plot_tvp_coefficients(res_dms_dma$X_list_all[[2]])

###############################################################################
# 8. STRESS DE TAUX
###############################################################################
# Paramètres pour les modèles de stress
lambda_grid_stress <- c(1, 0.995, 0.99)
kappa_stress <- 0.98
alpha_stress <- 0.99

###############################################################################
# 8.1 GÉNÉRATION DE TOUS LES SCÉNARIOS DE STRESS
###############################################################################
{
  all_stress_results <- generate_all_stress_scenarios(
    data_xts = data_xts,
    ns_factors_matrix = ns_factors_matrix,
    maturities = maturities,
    X0 = X0,
    P0 = P0,
    H0 = H0
  )
  all_scenarios <- all_stress_results$all_scenarios
  all_ns_factors_proj <- all_stress_results$all_ns_factors_proj
  
  cat("Scénarios générés:\n")
  print(names(all_scenarios))
}

###############################################################################
# 8.2 SCÉNARIO BASELINE
###############################################################################
{
  # Prévision baseline avec TVP-VAR-DMA 
  dma_baseline <- forecast_path(
    res_dms_dma$X_list_all, 
    res_dms_dma$pi, 
    ns_factors_matrix[nrow(ns_factors_matrix), ], 
    h = 20, 
    type = "dma"
  )
  
  # Reconstruction de la courbe des taux baseline
  baseline_scenario <- reconstruct_yield_curve(
    dma_baseline, 
    maturities, 
    lambda = 0.0609
  )
  colnames(baseline_scenario) <- names(maturities)
  
  # Dates futures pour les projections
  future_dates <- seq(
    from = data$Date[nrow(data)] + months(1), 
    by = "month", 
    length.out = nrow(baseline_scenario)
  )
  
  # Préparation pour l'affichage
  baseline_df <- data.frame(
    Dates = future_dates, 
    baseline_scenario
  )
  
  cat("Scénario baseline calculé pour", nrow(baseline_scenario), "périodes\n")
}

###############################################################################
# 8.3 AFFICHAGE DES SURFACES 3D POUR TOUS LES SCÉNARIOS
###############################################################################
{
  cat("Affichage du scénario baseline...\n")
  surface_baseline <- plot_surface(
    baseline_df, 
    maturities, 
    title = "Scénario Baseline - Projection TVP-VAR-DMA"
  )
  print(surface_baseline)
  
  # Surfaces pour tous les scénarios de stress
  cat("Affichage des scénarios de stress...\n")
  plot_all_surfaces(all_scenarios, maturities)
}

###############################################################################
# 8.4 ANALYSE DES FACTEURS NELSON-SIEGEL SOUS STRESS
###############################################################################
{
  key_scenarios <- c("up_parallel_tvpvar", "down_parallel_tvpvar", 
                     "up_short_tvpvar", "up_long_tvpvar")
  
  # Graphique des facteurs NS
  par(mfrow = c(2, 2))
  
  for (scenario in key_scenarios) {
    if (scenario %in% names(all_ns_factors_proj)) {
      factors_data <- all_ns_factors_proj[[scenario]]
      
      title_formatted <- gsub("_", " ", scenario)
      title_formatted <- gsub("tvpvar", "TVP-VAR", title_formatted)
      title_formatted <- tools::toTitleCase(title_formatted)
      
      matplot(index(factors_data), factors_data, 
              type = "l", lty = 1, lwd = 2,
              col = c("blue", "red", "green"),
              main = paste("Facteurs NS -", title_formatted),
              xlab = "Date", ylab = "Valeur du facteur")
      legend("topright", legend = c("Level", "Slope", "Curvature"),
             col = c("blue", "red", "green"), lty = 1, lwd = 2, cex = 0.8)
    }
  }
  
  par(mfrow = c(1, 1))
}