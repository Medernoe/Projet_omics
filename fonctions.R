#==================================================================================================
#Author : Noé Méderlet 
#contact : noe.mederlet@univ-rouen.fr
#github : https://github.com/Medernoe/Projet_omics
#organism : Master Bims M2, université de rouen 
#project : Création d'une application interactive dédiée à l'analyse de données transcriptomiques,
#développée dans le cadre d'un projet universitaire du Master 2 de Bioinformatique de l'Université de Rouen.
#==================================================================================================

# Classe la significativité des gènes selon la pvalue et logFC
significativity <- function(data, log2FC_cutoff, P_cutoff){ 
  
  data$Significance <- "Not significant"
  
  # Cas où pval = 1 => significativité ne depend que de log2FC
  if (P_cutoff == 1) {
    upregulated <- data$log2FC > log2FC_cutoff
    downregulated <- data$log2FC < -log2FC_cutoff
    
    data$Significance[upregulated] <- "Upregulated"
    data$Significance[downregulated] <- "Downregulated"
  } 
  else {
    # Cas classique avec seuil pval > 0
    upregulated <- data$pval < P_cutoff & data$log2FC >= log2FC_cutoff
    downregulated <- data$pval < P_cutoff & data$log2FC <= -log2FC_cutoff
    
    data$Significance[upregulated] <- "Upregulated"
    data$Significance[downregulated] <- "Downregulated"
  }
  
  # Définition des niveaux du facteur selon les catégories présentes
  present_levels <- unique(data$Significance)
  if ("Not significant" %in% present_levels) {
    data$Significance <- factor(data$Significance, levels=c("Not significant","Upregulated","Downregulated"))
  } else {
    data$Significance <- factor(data$Significance, levels=c("Upregulated","Downregulated"))
  }
  
  return(data)
}




# Plot volcano 
plot_volcano <- function(data, 
                         log2FC_cutoff,
                         P_cutoff, 
                         seuil_v, 
                         seuil_h, 
                         title = "None", 
                         highlight_row = NULL) {
  
  # Base du plot 
  p <- ggplot(data, aes(x = log2FC, y = -log10(pval), 
                        color = Significance, 
                        text = paste("Gene:", GeneName))) +  
    geom_point(alpha = 0.4, size = 0.8) +  
    labs(title = title,
         x = "Log2 Fold Change",
         y = "-log10(adj P-value)",
         color = "Significance") +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank()  # Réduit le nombre d'éléments à rendre
    )
  
  # Définir les couleurs selon le nombre de niveaux
  n_levels <- length(levels(data$Significance))
  if (n_levels == 3) {
    p <- p + scale_color_manual(values = c("gray", "red", "blue"))
  } else if (n_levels == 2) {
    p <- p + scale_color_manual(values = c("red", "blue"))
  }
  
  # Ligne horizontale (seuil p-value)
  if (seuil_h && P_cutoff > 0 && P_cutoff < 1) {
    ythr <- -log10(P_cutoff)
    p <- p + geom_hline(yintercept = ythr, linetype = "dashed", color = "gray40", linewidth = 0.5)
  }
  
  # Lignes verticales (seuil log2FC)
  if (seuil_v && log2FC_cutoff > 0 && log2FC_cutoff < max(data$log2FC, na.rm = TRUE)) {
    xthr <- c(-log2FC_cutoff, log2FC_cutoff)
    p <- p + geom_vline(xintercept = xthr, linetype = "dashed", color = "gray40", linewidth = 0.5)
  }
  
  # Vérification du highlight
  if (!is.null(highlight_row) && 
      length(highlight_row) > 0 && 
      !is.na(highlight_row) &&
      highlight_row > 0 && 
      highlight_row <= nrow(data)) {  
    
      # Ajouter un point plus grand pour le gène sélectionné
      p <- p + geom_point(data = data[highlight_row, , drop = FALSE],
                          aes(x = log2FC, y = -log10(pval)),
                          color = "purple", 
                          size = 2, 
                          shape = 16) +
        # Ajouter une étiquette avec le nom du gène
        geom_text(data = data[highlight_row, , drop = FALSE],
                  aes(x = log2FC, y = -log10(pval), label = GeneName),
                  vjust = -1.5,
                  hjust = 0.5,
                  size = 4,
                  fontface = "bold",
                  color = "black")
    }
  
  return(p)
}