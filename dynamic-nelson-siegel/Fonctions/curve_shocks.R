###############################################################################
# FONCTION EFFECTUANT LES PREVISIONS POUR UN SCENARIO DE STRESS
###############################################################################

run_stress_projection <- function(data_xts, ns_factors_matrix, maturities, 
                                  direction = c(1, -1), shock_type = "parallel", model = "tvpvar", 
                                  lambda_grid = c(1, 0.995, 0.99), lambda_ns = 0.0609,
                                  kappa = 0.98, alpha = 0.99, X0, P0, H0, h = 20,
                                  return_factors = FALSE) {
  
  # 1. On définit le type de choc
  shock_profile <- direction * switch(shock_type,
                                      "parallel" = rep(2.00, length(maturities)),
                                      "short" = approx(x = c(0.25, 2, 10, 20), y = c(2.5, 1.0, 0.0, 0.0),
                                                       xout = maturities, rule = 2)$y,
                                      "long" = approx(x = c(0.25, 1, 10, 20), y = c(0.0, 0.0, 1.5, 2.5),
                                                      xout = maturities, rule = 2)$y,
                                      stop("Type de choc non reconnu.")
  )
  
  # 2. On choque les derniers taux
  stressed_rates <- data_xts
  stressed_rates[NROW(stressed_rates), ] <- stressed_rates[NROW(stressed_rates), ] + shock_profile
  
  # 3. On refitte les facteurs NS sur les taux choqués
  ns_est <- estimate_NS_single_date(as.numeric(stressed_rates[NROW(stressed_rates), ]), maturities, lambda = lambda_ns)
  stressed_factors <- ns_factors_matrix
  stressed_factors[nrow(stressed_factors), ] <- ns_est$beta
  
  # 4. Forecast selon le modèle choisi
  if (model == "tvpvar") {
    res <- run_dms_dma(Y = stressed_factors,
                       lambda_grid = lambda_grid,
                       kappa = kappa,
                       alpha = alpha,
                       X0 = X0,
                       P0 = P0,
                       H0 = H0)
    forecast_factors <- forecast_path(res$X_list_all, res$pi,
                                      stressed_factors[nrow(stressed_factors), ], h = h, type = "dma")
    
  } else if (model == "vecm") {
    library(urca); library(vars)
    vecm <- ca.jo(stressed_factors, type = "trace", ecdet = "const", K = 2)
    vecm_var <- vec2var(vecm, r = 1)
    forecast_obj <- predict(vecm_var, n.ahead = h)
    forecast_factors <- sapply(forecast_obj$fcst, function(x) x[, "fcst"])
    colnames(forecast_factors) <- c("Level", "Slope", "Curvature")
  } else {
    stop("Modèle non reconnu. Choisir 'tvpvar' ou 'vecm'.")
  }
  
  # 5. On reconstitue la structure par terme à partir des facteurs projetés
  proj_yields <- reconstruct_yield_curve(forecast_factors, maturities, lambda = lambda_ns)
  colnames(proj_yields) <- names(maturities)
  
  future_dates <- seq(from = index(data_xts)[NROW(data_xts)] + months(1), by = "month", length.out = h)
  yield_proj <- data.frame(Dates = future_dates, proj_yields)
  
  if (return_factors) {
    proj_factors_xts <- xts(forecast_factors, order.by = future_dates)
    return(list(
      yields = yield_proj,
      factors = proj_factors_xts
    ))
  } else {
    return(yield_proj)
  }
}

###############################################################################
# FONCTION QUI GENERE TOUTES LES PROJECTIONS POUR TOUS LES SCENARIOS
###############################################################################
generate_all_stress_scenarios <- function(data_xts, ns_factors_matrix, maturities, X0, P0, H0) {
  scenarios <- list()
  ns_factors_proj <- list()
  
  shock_types <- c("parallel", "short", "long")
  directions <- c("up", "down")
  models <- c("tvpvar", "vecm")
  
  for (shock in shock_types) {
    for (dir in directions) {
      for (model in models) {
        key <- paste0(dir, "_", shock, "_", model)
        
        result <- run_stress_projection(
          data_xts, ns_factors_matrix, maturities,
          shock_type = shock, direction = ifelse(dir == "up", 1, -1),
          model = model,
          X0 = X0, P0 = P0, H0 = H0,
          return_factors = TRUE  
        )
        
        scenarios[[key]] <- result$yields
        ns_factors_proj[[key]] <- result$factors
      }
    }
  }
  
  return(list(
    all_scenarios = scenarios,
    all_ns_factors_proj = ns_factors_proj
  ))
}

###############################################################################
# FONCTION QUI CALCULE LES TAUX SUR DES MATURITIES PLUS FRAGMENTEES
###############################################################################
reconstruct_NS_curve_dense <- function(beta_matrix, lambda = 0.0609, max_maturity = 20, step = 1/12) {
  beta_matrix <- as.matrix(beta_matrix) 
  tau <- seq(step, max_maturity, by = step)
  
  # Matrice des loadings
  loadings <- cbind(
    Level = rep(1, length(tau)),
    Slope = (1 - exp(-lambda * tau)) / (lambda * tau),
    Curvature = ((1 - exp(-lambda * tau)) / (lambda * tau)) - exp(-lambda * tau)
  )
  # Pour chaque ligne, on calcule la courbe NS complète
  curve_matrix <- beta_matrix %*% t(loadings)  # dim = (dates) × (tau)
  
  colnames(curve_matrix) <- paste0("Y_", round(tau * 12), "M")  # ex : Y_12M, Y_24M...
  return(curve_matrix)
}
