###############################################################################
# ENSEMBLE DE FONCTIONS D'AFFICHAGE, GRAPHIQUES, ETC
###############################################################################

###############################################################################
# FONCTION POUR AFFICHER UNE SURFACE A PARTIR DE DONNEES IMPORTEES
###############################################################################
plot_surface <- function(data, maturities, 
                                title = "Surface des Taux d'Intérêt",
                                colorscale = "custom",
                                log_scale = TRUE,
                                font_size = 12) {
  
  # Passage au format matriciel
  rate_matrix <- as.matrix(data[, names(maturities)])
  
  # Echelle de couleur
  if (colorscale == "custom") {
    colors <- list(
      c(0, "blue"),
      c(0.25, "cyan"), 
      c(0.5, "yellow"),
      c(0.75, "orange"),
      c(1, "red")
    )
  } else {
    colors <- colorscale  #"Viridis", "Plasma", "Cividis", etc.
  }
  
  # Type d'échelle pour l'axe X (maturités)
  x_axis_type <- if(log_scale) "log" else "linear"
  
  # Surface 3D
  p1 <- plot_ly(
    x = ~maturities,                    # Axe X : Maturités
    y = ~data$Dates,                    # Axe Y : Dates
    z = ~rate_matrix,                   # Axe Z : Taux
    type = "surface",
    colorscale = colors
  ) %>%
    layout(
      title = title,
      scene = list(
        xaxis = list(title = "Maturité (années)", type = x_axis_type),
        yaxis = list(title = "Date"),
        zaxis = list(title = "Taux (%)")
      ),
      font = list(size = font_size)
    )
  
  return(p1)
}

###############################################################################
# FONCTION POUR AFFICHER TOUTES LES SURFACES DE LA LISTE PASSEE EN ARGUMENT
###############################################################################
plot_all_surfaces <- function(scenario_list, maturities) {
  for (name in names(scenario_list)) {
    title <- gsub("_", " ", toupper(name))
    p <- plot_surface(scenario_list[[name]], maturities, title = title)
    print(p)
  }
}

###############################################################################
# FONCTION POUR EXPORTER LES SCENARIOS SUR EXCEL (1 fichier/scénario)
###############################################################################
export_scenarios_to_excel <- function(scenario_list, filepath_prefix = "stress_projection") {
  if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl")
  library(writexl)
  
  for (name in names(scenario_list)) {
    df <- data.frame(Date = index(scenario_list[[name]]), coredata(scenario_list[[name]]))
    writexl::write_xlsx(df, path = paste0(filepath_prefix, "_", name, ".xlsx"))
  }
}

###############################################################################
# FONCTION POUR EXPORTER LES SCENARIOS SUR EXCEL (1 fichier, une feuille/scénario)
###############################################################################
export_stress_scenarios_excel_multisheet <- function(scenario_list, filename = "Chocs transitoires.xlsx") {
  if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl")
  library(writexl)
  
  # Construire une liste de data.frames nommés pour les sheets
  sheets <- list()
  
  for (name in names(scenario_list)) {
    df <- data.frame(Date = index(scenario_list[[name]]), coredata(scenario_list[[name]]))
    
    # Nettoyage du nom pour Excel
    sheet_name <- gsub("_", " ", name)            # underscore -> espace
    sheet_name <- gsub("tvpvar", "TVP", sheet_name)
    sheet_name <- gsub("vecm", "VECM", sheet_name)
    sheet_name <- tools::toTitleCase(sheet_name)  # majuscules propres
    
    sheets[[sheet_name]] <- df
  }
  
  writexl::write_xlsx(sheets, path = filename)
}