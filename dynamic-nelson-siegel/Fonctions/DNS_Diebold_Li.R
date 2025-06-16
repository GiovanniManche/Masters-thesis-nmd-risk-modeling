################################################################################
#FONCTION POUR TROUVER L'ORDRE OPTIMAL D'UN ARMA(P,Q) A L'AIDE DES CRITERES 
# D'INFORMATION ET L'ESTIMER
################################################################################
estimate_optimal_arma <- function(variable, max_p = 5, max_q = 5,
                              information_criterion = "aic"){
  
  if (!information_criterion %in% c("aic", "bic")){
    stop("Le critère d'information doit être AIC ou BIC")
  }
  
  results <- data.frame()
  best_information_criterion <- Inf 
  best_p <- 0
  best_q <- 0
  best_model <- NULL
  cat("Critère :", toupper(information_criterion), "\n")
  
  # Boucle sur l'ordre AR
  for(p in 0:max_p) {
    # Boucle sur l'ordre MA
    for(q in 0:max_q) {
      
      if(p == 0 & q == 0) next
      
      tryCatch({
        # Estimation du modèle ARMA(p,q)
        model <- Arima(variable, order = c(p, 0, q))
        
         
        ic_value <- switch(information_criterion,
                           "aic" = AIC(model),
                           "bic" = BIC(model))
        
        results <- rbind(results, data.frame(
          p = p, q = q, 
          aic = AIC(model),
          bic = BIC(model),
          loglik = model$loglik,
          sigma2 = model$sigma2
        ))
        
        # Mise à jour du meilleur modèle
        if(ic_value < best_information_criterion) {
          best_information_criterion <- ic_value
          best_p <- p
          best_q <- q
          best_model <- model
        }
        
      }, error = function(e) {
        # En cas d'erreur d'estimation, on passe au suivant
        print("ERREUR")
        NULL
      })
    }
  }
  
  results <- results[order(results[[information_criterion]]), ]
  cat("\nTop 5 modèles selon", toupper(information_criterion), " :\n")
  print(head(results[, c("p", "q", information_criterion)], 5))
  
  # Retourner les résultats
  return(list(
    best_p = best_p,
    best_q = best_q,
    best_model = best_model,
    best_ic = best_information_criterion,
    all_results = results,
    criterion = information_criterion
  ))
}

################################################################################
# FONCTION GENERALE DE TESTS SUR LES RESIDUS
################################################################################
test_residuals <- function(residuals, series_name = "Series", 
                           lags_ljb = 10, lags_arch = 10,
                           alpha = 0.05) {
  
  if(!is.numeric(residuals)) {
    stop("Les résidus doivent être un vecteur numérique")
  }
  
  if(any(is.na(residuals))) {
    cat("ATTENTION : Des valeurs manquantes ont été détectées. Suppression automatique\n")
    residuals <- residuals[!is.na(residuals)]
  }
  
  n_obs <- length(residuals)
  if(n_obs < 10) {
    stop("Pas assez d'observations (minimum 10)")
  }
  
  max_lags_ljb <- min(lags_ljb, floor(n_obs/4))
  max_lags_arch <- min(lags_arch, floor(n_obs/10))
  
  test_results <- list()
  test_results$series_name <- series_name
  test_results$n_obs <- n_obs
  test_results$alpha <- alpha
  
 
  cat("TESTS DE DIAGNOSTIC DES RESIDUS\n")
 
  cat("Série :", series_name, "\n")
  cat("Observations :", n_obs, "\n")
  
  
  ##############################################################################
  # 1. TEST DE LJUNG-BOX
  ##############################################################################
  
  ljung_test <- Box.test(residuals, lag = max_lags_ljb, type = "Ljung-Box")
  test_results$ljung_box <- list(
    statistic = ljung_test$statistic,
    p_value = ljung_test$p.value,
    lags = max_lags_ljb,
    reject_h0 = ljung_test$p.value < alpha
  )
  
  cat("1. TEST DE LJUNG-BOX \n")
  cat("   Statistique :", round(ljung_test$statistic, 4), "\n")
  cat("   p-value :", round(ljung_test$p.value, 4), "\n")
  cat("   Conclusion:", ifelse(ljung_test$p.value < alpha, 
                               "Rejet de H0 : Autocorrélation détectée !", 
                               "Non rejet de H0 : Pas d'autocorrélation"), "\n\n")
  
  ##############################################################################
  # 2. TEST DE JARQUE-BERA 
  ##############################################################################
  
  jb_test <- jarque.bera.test(residuals)
  test_results$jarque_bera <- list(
    statistic = jb_test$statistic,
    p_value = jb_test$p.value,
    reject_h0 = jb_test$p.value < alpha
  )
  
  cat("2. TEST DE JARQUE-BERA \n")
  cat("   Statistique :", round(jb_test$statistic, 4), "\n")
  cat("   p-value :", round(jb_test$p.value, 4), "\n")
  cat("   Conclusion:", ifelse(jb_test$p.value < alpha, 
                               "Rejet de H0 : Résidus non-normaux", 
                               "Non rejet de H0 : Résidus normaux"), "\n\n")
  
  ##############################################################################
  # 3. TEST D'ABSENCE D'EFFET ARCH
  ##############################################################################
  
  arch_test <- ArchTest(residuals, lags = max_lags_arch)
  test_results$arch_test <- list(
    statistic = as.numeric(arch_test$statistic),
    p_value = arch_test$p.value,
    lags = max_lags_arch,
    method = "ARCH LM Test",
    reject_h0 = arch_test$p.value < alpha
  )
  
  cat("3. TEST D'ABSENCE D'EFFET ARCH \n")
  cat("   Statistique :", round(as.numeric(arch_test$statistic), 4), "\n")
  cat("   p-value :", round(arch_test$p.value, 4), "\n")
  cat("   Lags testés :", max_lags_arch, "\n")
  cat("   Conclusion :", ifelse(arch_test$p.value < alpha, 
                                "Rejet de H0 : présence d'un effet ARCH", 
                                "Non rejet de H0 : absence d'effet ARCH"), "\n\n")
  
  ##############################################################################
  # 4. RESUME FINAL
  ##############################################################################
  
 
  cat("RÉSUMÉ DES TESTS\n")
 
  
  passed_tests <- 0
  total_tests <- 3
  
  if(!test_results$ljung_box$reject_h0) passed_tests <- passed_tests + 1
  cat("Test de Ljung-Box :", ifelse(!test_results$ljung_box$reject_h0, "PASS", "FAIL"), "\n")
  
  if(!test_results$jarque_bera$reject_h0) passed_tests <- passed_tests + 1
  cat("Test de Jarque-Bera :", ifelse(!test_results$jarque_bera$reject_h0, "PASS", "FAIL"), "\n")
  
  if(!test_results$arch_test$reject_h0) passed_tests <- passed_tests + 1
  cat("Test d'absence d'effet ARCH :", ifelse(!test_results$arch_test$reject_h0, "PASS", "FAIL"), "\n")
  
  
  cat("Score global :", passed_tests, "/", total_tests, 
      "(", round(100*passed_tests/total_tests, 1), "%)\n")
  
  # Retourner les résultats
  test_results$passed_tests <- passed_tests
  test_results$total_tests <- total_tests
  test_results$success_rate <- passed_tests / total_tests
  
  return(invisible(test_results))
}


################################################################################
# FONCTION POUR TROUVER L'ORDRE OPTIMAL D'UN VAR A L'AIDE DES CRITERES 
# D'INFORMATION ET L'ESTIMER
################################################################################

estimate_optimal_var <- function(data_ts, max_p = 10, 
                                 information_criterion = "aic") {
  
  if (!information_criterion %in% c("aic", "bic")) {
    stop("Le critère d'information doit être 'aic' ou 'bic'")
  }
  
  if (!is.data.frame(data_ts) && !is.ts(data_ts) && !is.matrix(data_ts)) {
    stop("Les données doivent être une série temporelle multivariée (ts ou data.frame)")
  }
  
  results <- data.frame()
  best_ic_value <- Inf
  best_p <- NA
  best_model <- NULL
  
  cat("Critère :", toupper(information_criterion), "\n")
  
  for (p in 1:max_p) {
    tryCatch({
      var_model <- VAR(y = data_ts, p = p, type = "const")
      
      ic_value <- switch(information_criterion,
                         "aic" = AIC(var_model),
                         "bic" = BIC(var_model))
      
      results <- rbind(results, data.frame(
        p = p,
        aic = AIC(var_model),
        bic = BIC(var_model)
      ))
      
      if (ic_value < best_ic_value) {
        best_ic_value <- ic_value
        best_p <- p
        best_model <- var_model
      }
      
    }, error = function(e) {
      message(paste("Erreur à l'estimation pour p =", p, ":", e$message))
    })
  }
  
  results <- results[order(results[[information_criterion]]), ]

  cat("\nTop 5 modèles selon", toupper(information_criterion), ":\n")
  print(head(results[, c("p", information_criterion)], 5))
  
  # Retour
  return(list(
    best_p = best_p,
    best_model = best_model,
    best_ic = best_ic_value,
    all_results = results,
    criterion = information_criterion
  ))
}


###############################################################################
# FONCTION QUI EFFECTUE DES PREVISIONS ROULANTES ET CALCULE LA RMSE ASSOCIEE
###############################################################################
rolling_rmse_dm_evaluation <- function(data, K_jo, window_size = 120, horizons = c(1, 3, 6, 12)) {
  library(forecast)
  library(vars)
  library(urca)
  
  n <- nrow(data)
  start_idx <- window_size + max(horizons)
  rmse_summary <- data.frame()
  dm_summary <- list()
  
  for (h in horizons) {
    cat("== Horizon:", h, "mois ==\n")
    
    n_forecasts <- n - start_idx + 1
    forecasts_ar <- matrix(NA, n_forecasts, 3)
    forecasts_var <- matrix(NA, n_forecasts, 3)
    forecasts_vecm <- matrix(NA, n_forecasts, 3)
    actuals <- matrix(NA, n_forecasts, 3)
    
    for (i in 1:n_forecasts) {
      end_train <- start_idx + i - 1
      forecast_idx <- end_train + h
      if (forecast_idx > n) break
      
      train_data <- data[1:end_train, ]
      actual_value <- data[forecast_idx, ]
      
      tryCatch({
        # AR
        ar_level <- Arima(train_data[,"Level"], order = c(1,0,0))
        ar_slope <- Arima(train_data[,"Slope"], order = c(1,0,0))
        ar_curv  <- Arima(train_data[,"Curvature"], order = c(1,0,0))
        forecasts_ar[i, ] <- c(
          forecast(ar_level, h = h)$mean[h],
          forecast(ar_slope, h = h)$mean[h],
          forecast(ar_curv,  h = h)$mean[h]
        )
        
        # VAR
        var_model <- VAR(train_data, p = 1, type = "const")
        var_pred <- predict(var_model, n.ahead = h)
        forecasts_var[i, ] <- c(
          var_pred$fcst$Level[h, 1],
          var_pred$fcst$Slope[h, 1],
          var_pred$fcst$Curvature[h, 1]
        )
        
        # VECM
        johansen <- ca.jo(train_data, K = K_jo, type = "trace", ecdet = "const")
        vecm_model <- vec2var(johansen, r = 1)
        vecm_pred <- predict(vecm_model, n.ahead = h)
        forecasts_vecm[i, ] <- c(
          vecm_pred$fcst$Level[h, 1],
          vecm_pred$fcst$Slope[h, 1],
          vecm_pred$fcst$Curvature[h, 1]
        )
        
        actuals[i, ] <- as.numeric(actual_value)
      }, error = function(e) {})
    }
    
    valid <- complete.cases(forecasts_ar, forecasts_var, forecasts_vecm, actuals)
    fa <- forecasts_ar[valid, ]; fv <- forecasts_var[valid, ]; fe <- forecasts_vecm[valid, ]
    act <- actuals[valid, ]
    
    # RMSE
    rmse_ar <- sqrt(colMeans((fa - act)^2, na.rm = TRUE))
    rmse_var <- sqrt(colMeans((fv - act)^2, na.rm = TRUE))
    rmse_vecm <- sqrt(colMeans((fe - act)^2, na.rm = TRUE))
    
    rmse_summary <- rbind(rmse_summary, data.frame(
      Horizon = h,
      AR_Level = rmse_ar[1], AR_Slope = rmse_ar[2], AR_Curvature = rmse_ar[3],
      VAR_Level = rmse_var[1], VAR_Slope = rmse_var[2], VAR_Curvature = rmse_var[3],
      VECM_Level = rmse_vecm[1], VECM_Slope = rmse_vecm[2], VECM_Curvature = rmse_vecm[3]
    ))
    
    # Diebold-Mariano
    dm_test <- function(e1, e2) {
      sapply(1:3, function(j) {
        if (nrow(e1) < 12) return(NA)
        res <- tryCatch(
          dm.test(e1[, j], e2[, j], alternative = "two.sided")$p.value,
          error = function(e) NA
        )
        return(res)
      })
    }
    
    dm_summary[[paste0("h", h)]] <- list(
      AR_vs_VAR = dm_test((act - fa)^2, (act - fv)^2),
      AR_vs_VECM = dm_test((act - fa)^2, (act - fe)^2),
      VAR_vs_VECM = dm_test((act - fv)^2, (act - fe)^2)
    )
  }
  
  # Moyenne des colonnes RMSE
  rmse_cols <- names(rmse_summary)[-1]  # exclut 'Horizon'
  avg_row <- c(Horizon = NA, colMeans(rmse_summary[, rmse_cols], na.rm = TRUE))
  rmse_summary <- rbind(rmse_summary, avg_row)
  rownames(rmse_summary)[nrow(rmse_summary)] <- "Average"
  
  cat("\n=== RMSE Résumé ===\n")
  print(round(rmse_summary, 4))
  
  return(list(rmse = rmse_summary, dm = dm_summary))
}