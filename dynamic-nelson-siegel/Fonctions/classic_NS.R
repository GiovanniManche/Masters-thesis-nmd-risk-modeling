################################################################################
# FONCTION D'ESTIMATION DU MODELE NELSON-SIEGEL CLASSIQUE A DATE DONNEE
################################################################################

estimate_NS_single_date <- function(yields, maturities, lambda = 0.0609){
  # Construction des loadings Nelson-Siegel
  loadings <- cbind(
    Level = rep(1, length(maturities)),  # Level : toutes maturités = 1
    Slope = (1 - exp(-lambda * maturities)) / (lambda * maturities),  # Slope
    Curvature = ((1 - exp(-lambda * maturities)) / (lambda * maturities)) - exp(-lambda * maturities)  # Curvature
  )
  
  # Régression OLS : yields = loadings * betas + erreur
  reg <- lm(yields ~ loadings - 1)  # -1 pour supprimer l'intercept
  
  return(list(
    beta = coef(reg),
    fitted = fitted(reg),
    residuals = residuals(reg),
    r_squared = summary(reg)$r.squared,
    loadings = loadings
  ))
}

################################################################################
# FONCTION PERMETTANT DE COMPARER GRAPHIQUEMENT LA STRUCTURE MOYENNE REELLE ET 
# ESTIMEE
################################################################################
plot_NS_fit_avg <- function(data_xts, ns_factors, maturities,
                            lambda = 0.0609,
                            col_obs = "red", col_fit = "blue",
                            legend_pos = "bottomright",
                            main_title = "Courbe des taux : Moyenne observée vs Estimée (Nelson-Siegel)",
                            xlab = "Maturité (années)",
                            ylab = "Taux moyen (%)") {
  
  observed_mean <- apply(data_xts, 2, mean, na.rm = TRUE)
  mean_factors <- apply(ns_factors, 2, mean, na.rm = TRUE)
  
  loadings <- cbind(
    Level = rep(1, length(maturities)),
    Slope = (1 - exp(-lambda * maturities)) / (lambda * maturities),
    Curvature = ((1 - exp(-lambda * maturities)) / (lambda * maturities)) - exp(-lambda * maturities)
  )
  
  estimated_mean <- loadings %*% mean_factors
  df <- data.frame(
    Maturite = maturities,
    Observee = as.numeric(observed_mean),
    Estimee_NS = as.numeric(estimated_mean)
  )
  
  x11()
  plot(df$Maturite, df$Observee, type = "p",
       main = main_title, xlab = xlab, ylab = ylab,
       ylim = range(c(df$Observee, df$Estimee_NS)))
  lines(df$Maturite, df$Estimee_NS, col = col_fit, lwd = 2)
  points(df$Maturite, df$Estimee_NS, col = col_fit)
  
  legend(legend_pos,
         legend = c("Moyenne observée", "Moyenne estimée NS"),
         col = c(col_obs, col_fit),
         lty = c(0, 1), lwd = c(0, 2))
}

plot_NS_fit_dates <- function(data_xts, ns_factors, maturities, dates,
                              lambda = 0.0609,
                              col_obs = "red", col_fit = "blue",
                              legend_pos = "bottomright",
                              main_prefix = "Ajustement NS - ",
                              xlab = "Maturité (années)",
                              ylab = "Taux (%)") {
  
  # Fonction interne pour reconstruire la courbe NS
  reconstruct_curve <- function(date_index) {
    factors <- as.numeric(ns_factors[date_index, ])
    loadings <- cbind(
      Level = rep(1, length(maturities)),
      Slope = (1 - exp(-lambda * maturities)) / (lambda * maturities),
      Curvature = ((1 - exp(-lambda * maturities)) / (lambda * maturities)) - exp(-lambda * maturities)
    )
    estimated_curve <- loadings %*% factors
    list(
      date = index(data_xts)[date_index],  # Correction : utiliser index()
      observed = as.numeric(data_xts[date_index, ]),
      estimated = as.numeric(estimated_curve)
    )
  }
  
  results <- lapply(dates, reconstruct_curve)
  
  x11()
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
  
  for (res in results) {
    # Calcul RMSE pour cette date
    rmse_date <- sqrt(mean((res$observed - res$estimated)^2))
    
    plot(maturities, res$observed,  # Correction : utiliser maturities
         type = "p", pch = 16, col = col_obs, cex = 1.2,
         main = paste0(main_prefix, format(res$date, "%Y-%m-%d")),
         xlab = xlab, ylab = ylab,
         ylim = range(c(res$observed, res$estimated)))
    
    lines(maturities, res$estimated, col = col_fit, lwd = 2)
    points(maturities, res$estimated, pch = 17, col = col_fit, cex = 1.2)
    
    legend(legend_pos,
           legend = c("Observé", "Estimé NS"),
           col = c(col_obs, col_fit, "white"),
           pch = c(16, 17, NA),
           lty = c(0, 1, 0), lwd = c(0, 2, 0),
           cex = 0.8)
  }
  
  par(mfrow = c(1, 1))
}