###############################################################################
# FONCTION QUI PROCEDE AU FILTRE DE KALMAN RECURSIF AVEC UNE MATRICE DE COVARIANCE
# COMME FONCTION D'UN FORGETTING FACTOR LAMBDA
###############################################################################
kalman_filter_forgetting_factors <- function(Y, lambda = 0.99, kappa = 1.0,
                                      H0 = NULL, X0 = NULL, P0 = NULL) {
  T <- nrow(Y)
  k <- ncol(Y)
  d <- k^2
  
  # Initialisation
  if (is.null(H0)) H_t <- diag(0.1, k) else H_t <- H0
  if (is.null(X0)) X_t <- rep(0, d) else X_t <- X0
  if (is.null(P0)) P_t <- diag(10, d) else P_t <- P0
  
  # Stockage
  X_list <- vector("list", T)
  P_list <- vector("list", T)
  F_list <- vector("list", T)
  S_list <- vector("list", T)
  H_list <- vector("list", T)
  
  # Filtre récursif
  for (t in 2:T) {
    Y_t <- Y[t, ]
    Y_tm1 <- Y[t - 1, , drop = FALSE]
    
    O_t <- kronecker(Y_tm1, diag(k))
    M_t <- diag(d)
    C_t <- rep(0, d)
    D_t <- rep(0, k)
    
    # 1. Prévision t|t-1 avec forgetting factor
    X_pred <- M_t %*% X_t + C_t
    P_pred <- (1 / lambda) * (M_t %*% P_t %*% t(M_t))  
    F_t <- O_t %*% X_pred + D_t
    S_t <- O_t %*% P_pred %*% t(O_t) + H_t
    e_t <- Y_t - F_t
    
    # 2. Mise à jour
    # Gain de Kalman
    K_t <- P_pred %*% t(O_t) %*% solve(S_t)
    X_t <- X_pred + K_t %*% e_t
    P_t <- P_pred - K_t %*% O_t %*% P_pred
    
    # Mise à jour de H_t via kappa
    H_t <- kappa * H_t + (1 - kappa) * (e_t %*% t(e_t))
    
    X_list[[t]] <- X_t
    P_list[[t]] <- P_t
    F_list[[t]] <- F_t
    S_list[[t]] <- S_t
    H_list[[t]] <- H_t
  }
  
  return(list(
    X = X_list,
    P = P_list,
    F = F_list,
    S = S_list,
    H = H_list
  ))
}

###############################################################################
# FONCTION QUI CALCULE LA PERFORMANCE PREDICTIVE DU DMA/DMS SUR EXPANDING WINDOW
###############################################################################
run_dms_dma <- function(Y, 
                        lambda_grid = c(1.00, 0.99, 0.97, 0.95),
                        kappa = 0.98,
                        alpha = 0.99,
                        X0 = NULL,
                        P0 = NULL,
                        H0 = NULL) {
  

  
  T <- nrow(Y)
  k <- ncol(Y)
  d <- k^2
  J <- length(lambda_grid)
  
  if (is.null(X0)) X0 <- rep(0, d)
  if (is.null(P0)) P0 <- diag(10, d)
  if (is.null(H0)) H0 <- diag(0.1, k)

  # 1. Filtre de Kalman pour chaque modèle
  kalman_outputs <- vector("list", J)
  for (j in 1:J) {
    kalman_outputs[[j]] <- kalman_filter_forgetting_factors(
      Y = Y,
      lambda = lambda_grid[j],
      kappa = kappa,
      X0 = X0,
      P0 = P0,
      H0 = H0
    )
  }

  # Initialisation des objets pour DMA/DMS
  pi_t <- rep(1 / J, J)
  pi_store <- matrix(NA, T, J)
  pred_store <- array(NA, dim = c(T, k, J))
  pred_dma <- matrix(NA, T, k)
  pred_dms <- matrix(NA, T, k)
  
  # 2. Calcul des probabilités
  # Mécanisme récursif
  for (t in 2:T) {
    y_t <- Y[t, ]
    pred_dens <- rep(NA, J)
    # On récupère les résultats des prévisions de chaque modèle
    for (j in 1:J) {
      res <- kalman_outputs[[j]]
      F_tj <- res$F[[t]]
      S_tj <- res$S[[t]]
      
      pred_store[t, , j] <- F_tj
      pred_dens[j] <- dmvnorm(y_t, mean = as.vector(F_tj), sigma = S_tj, log = FALSE)
    }

    # On met à jour les probabilités (= les poids)
    pi_pred <- pi_t^alpha
    pi_pred <- pi_pred / sum(pi_pred)
    numerateur <- pi_pred * pred_dens
    pi_t <- numerateur / sum(numerateur)
    
    pi_store[t, ] <- pi_t
    j_star <- which.max(pi_t)
    pred_dms[t, ] <- pred_store[t, , j_star]
    pred_dma[t, ] <- as.vector(pred_store[t, , ] %*% pi_t)
    
  }

  # Récupération des matrices de coefficients pour la prévision
  X_list_all <- lapply(kalman_outputs, function(res) res$X)
  
  return(list(
    pi = pi_store,
    pred_dms = pred_dms,
    pred_dma = pred_dma,
    pred_store = pred_store,
    lambdas = lambda_grid,
    X_list_all = X_list_all
  ))
}

###############################################################################
# FONCTION QUI PROCEDE A LA STRATEGIE SEQUENTIELLE POUR LE TEST ADF
###############################################################################
evaluate_dms_dma_performance <- function(res_dms_dma, ns_factors_matrix, 
                                         window_size = 120, horizons = c(1, 3, 6, 12)) {
  
  n <- nrow(ns_factors_matrix)
  k <- ncol(ns_factors_matrix)  # 3 facteurs
  start_idx <- window_size + max(horizons)
  rmse_summary <- data.frame()
  
  for (h in horizons) {
    cat("== Horizon:", h, "mois ==\n")
    
    n_forecasts <- n - start_idx - h + 1
    if (n_forecasts <= 0) {
      warning(paste("Pas assez de données pour horizon", h))
      next
    }
    
    forecasts_dms <- matrix(NA, n_forecasts, k)
    forecasts_dma <- matrix(NA, n_forecasts, k)
    actuals <- matrix(NA, n_forecasts, k)
    
    for (i in 1:n_forecasts) {
      forecast_origin <- start_idx + i - 1  # Point de prévision
      target_date <- forecast_origin + h     # Date cible
      
      if (target_date <= n) {
        forecasts_dms[i, ] <- res_dms_dma$pred_dms[forecast_origin, ]
        forecasts_dma[i, ] <- res_dms_dma$pred_dma[forecast_origin, ]
        actuals[i, ] <- ns_factors_matrix[target_date, ]
      }
    }
    
    # Nettoyage des données manquantes
    valid_rows <- complete.cases(forecasts_dms, forecasts_dma, actuals)
    f_dms <- forecasts_dms[valid_rows, , drop = FALSE]
    f_dma <- forecasts_dma[valid_rows, , drop = FALSE]
    act <- actuals[valid_rows, , drop = FALSE]
    
    if (nrow(f_dms) == 0) {
      warning(paste("Aucune prévision valide pour horizon", h))
      next
    }
    
    # Calcul RMSE par facteur
    rmse_dms <- sqrt(colMeans((f_dms - act)^2, na.rm = TRUE))
    rmse_dma <- sqrt(colMeans((f_dma - act)^2, na.rm = TRUE))

    # Stockage des résultats
    rmse_summary <- rbind(rmse_summary, data.frame(
      Horizon = h,
      DMS_Level = rmse_dms[1],
      DMS_Slope = rmse_dms[2], 
      DMS_Curvature = rmse_dms[3],
      DMA_Level = rmse_dma[1],
      DMA_Slope = rmse_dma[2],
      DMA_Curvature = rmse_dma[3],
      N_forecasts = nrow(f_dms)
    ))
  }
  
  rmse_cols <- c("DMS_Level", "DMS_Slope", "DMS_Curvature",
                 "DMA_Level", "DMA_Slope", "DMA_Curvature")
  avg_row <- c(Horizon = NA,
               colMeans(rmse_summary[, rmse_cols], na.rm = TRUE),
               N_forecasts = NA)
  
  rmse_summary <- rbind(rmse_summary, avg_row)
  rownames(rmse_summary)[nrow(rmse_summary)] <- "Average"
  
  return(rmse_summary)
}




###############################################################################
# FONCTION DE PRÉVISION RÉCURSIVE PONCTUELLE (HORIZON h)
###############################################################################
forecast_recursive <- function(X_list_list, pi_matrix, y_t, h, type = c("dma", "dms")) {
  type <- match.arg(type)
  
  J <- length(X_list_list)
  k <- length(y_t)
  T <- nrow(pi_matrix)
  
  A_t_list <- lapply(1:J, function(j) {
    matrix(X_list_list[[j]][[T]], nrow = k, byrow = FALSE)
  })
  
  y_forecast <- y_t
  
  for (step in 1:h) {
    pred_j <- sapply(1:J, function(j) {
      A_t_list[[j]] %*% y_forecast
    })
    
    y_forecast <- switch(type,
                         dma = as.vector(pred_j %*% pi_matrix[T, ]),
                         dms = as.vector(pred_j[, which.max(pi_matrix[T, ])])
    )
  }
  
  return(y_forecast)
}

###############################################################################
# FONCTION DE PRÉVISION RÉCURSIVE COMPLÈTE (TRAJECTOIRE SUR h PÉRIODES)
###############################################################################
forecast_path <- function(X_list_list, pi_matrix, y_t, h, type = c("dma", "dms")) {
  type <- match.arg(type)
  
  J <- length(X_list_list)
  k <- length(y_t)
  T <- nrow(pi_matrix)
  
  A_t_list <- lapply(1:J, function(j) {
    matrix(X_list_list[[j]][[T]], nrow = k, byrow = FALSE)
  })
  
  y_forecast <- y_t
  path <- matrix(NA, nrow = h, ncol = k)
  
  for (step in 1:h) {
    pred_j <- sapply(1:J, function(j) {
      A_t_list[[j]] %*% y_forecast
    })
    
    y_forecast <- switch(type,
                         dma = as.vector(pred_j %*% pi_matrix[T, ]),
                         dms = as.vector(pred_j[, which.max(pi_matrix[T, ])])
    )
    
    path[step, ] <- y_forecast
  }
  
  if (!is.null(colnames(pi_matrix))) {
    colnames(path) <- paste0("Factor_", 1:k)
  } else {
    colnames(path) <- c("Level", "Slope", "Curvature")[1:k]
  }
  
  rownames(path) <- paste0("T+", 1:h)
  
  return(path)
}

###############################################################################
# FONCTION DE RECONSTRUCTION DE LA COURBE DES TAUX À PARTIR DES FACTEURS NS
###############################################################################
reconstruct_yield_curve <- function(ns_path, maturities, lambda = 0.0609) {
  if (!is.matrix(ns_path) && !is.data.frame(ns_path)) {
    stop("ns_path doit être une matrice ou un data.frame")
  }
  
  if (ncol(ns_path) != 3) {
    stop("ns_path doit avoir exactement 3 colonnes (Level, Slope, Curvature)")
  }
  
  h <- nrow(ns_path)
  n <- length(maturities)
  
  compute_loadings <- function(tau, lambda) {
    l1 <- 1
    l2 <- (1 - exp(-lambda * tau)) / (lambda * tau)
    l3 <- l2 - exp(-lambda * tau)
    return(c(l1, l2, l3))
  }
  
  L <- t(sapply(maturities, function(tau) {
    compute_loadings(tau, lambda)
  }))
  
  if (nrow(L) != n || ncol(L) != 3) {
    stop("Erreur dans la construction de la matrice des loadings")
  }
  
  yield_curves <- matrix(NA, nrow = h, ncol = n)
  
  for (t in 1:h) {
    betas <- as.numeric(ns_path[t, ])
    yield_curves[t, ] <- as.vector(L %*% betas)
  }
  
  if (!is.null(names(maturities))) {
    colnames(yield_curves) <- names(maturities)
  } else {
    colnames(yield_curves) <- paste0("Tau_", maturities, "Y")
  }
  
  rownames(yield_curves) <- paste0("T+", 1:h)
  
  return(yield_curves)
}

###############################################################################
# FONCTION DE VISUALISATION DES COEFFICIENTS TVP-VAR
###############################################################################
plot_tvp_coefficients <- function(X_list, var_names = c("Level", "Slope", "Curvature"),
                                  title_suffix = "", lambda_value = NULL) {
  if (length(X_list) == 0) {
    stop("X_list ne peut pas être vide")
  }
  
  k <- length(var_names)
  T <- length(X_list)
  d <- k^2
  
  if (!is.null(X_list[[T]]) && length(X_list[[T]]) != d) {
    warning("Dimension incohérente détectée dans X_list")
  }
  
  coef_mat <- matrix(NA, nrow = T, ncol = d)
  
  for (t in 1:T) {
    if (!is.null(X_list[[t]]) && length(X_list[[t]]) == d) {
      A_t <- matrix(X_list[[t]], nrow = k, byrow = FALSE)
      coef_mat[t, ] <- as.vector(A_t)
    }
  }
  
  df <- data.frame(time = 1:T, coef_mat)
  
  coef_names <- paste0(
    rep(var_names, each = k),
    "_on_", 
    rep(var_names, times = k)
  )
  colnames(df)[-1] <- coef_names
  
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    stop("Le package 'reshape2' est requis pour cette fonction")
  }
  
  df_long <- reshape2::melt(df, id.vars = "time", 
                            variable.name = "coefficient", 
                            value.name = "value")
  
  main_title <- "Évolution des coefficients dynamiques (TVP-VAR)"
  if (!is.null(lambda_value)) {
    main_title <- paste0(main_title, ", λ = ", lambda_value)
  }
  if (title_suffix != "") {
    main_title <- paste0(main_title, " - ", title_suffix)
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Le package 'ggplot2' est requis pour cette fonction")
  }
  
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = time, y = value, color = coefficient)) +
    ggplot2::geom_line(size = 0.8, alpha = 0.8) +
    ggplot2::labs(
      title = main_title,
      x = "Temps (périodes)",
      y = "Valeur du coefficient",
      color = "Coefficient"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 11)
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(ncol = 3))
  
  return(p)
}