#==================================================================================================
#Author : Noé Méderlet 
#contact : noe.mederlet@univ-rouen.fr
#github : https://github.com/Medernoe/Projet_omics
#organism : Master Bims M2, université de rouen 
#project : Création d'une application interactive dédiée à l’analyse de données transcriptomiques,
#développée dans le cadre d’un projet universitaire du Master 2 de Bioinformatique de l’Université de Rouen.
#==================================================================================================#fonction

#fonction qui selon la pvalue et logFC classe la significativité des gènes 
significativity <- function(data, log2FC_cutoff, P_cutoff){ 
  
  # Initialisation
  data$Significance <- "Not significant"
  
  # Cas où p_value = 0 => tous significatifs
  if (P_cutoff == 0) {
    data$Significance[data$log2FC > 0] <- "Upregulated"
    data$Significance[data$log2FC < 0] <- "Downregulated"
  } 
  else {
    # Cas classique avec seuil p_value > 0
    data$Significance[data$p_value < P_cutoff & data$log2FC >= log2FC_cutoff] <- "Upregulated"
    data$Significance[data$p_value < P_cutoff & data$log2FC <= -log2FC_cutoff] <- "Downregulated"
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


#fonction qui selon la pvalue et logFC classe la significativité des gènes 
plot_volcano <- function(data, log2FC_cutoff, P_cutoff,
                         title = "Volcano Plot",
                         seuil_v = FALSE, seuil_h = FALSE) {
  
  # Base du plot
  p <- ggplot(data, aes(x = log2FC, y = -log10(p_value), color = Significance)) +
    geom_point(alpha = 0.5, size = 1) +
    labs(title = title,
         x = "Log2 Fold Change",
         y = "-log10(adj P-value)",
         color = "Significance") +
    theme_minimal()
  
  # Définir les couleurs selon le nombre de niveaux
  n_levels <- length(levels(data$Significance))
  if (n_levels == 3) {
    p <- p + scale_color_manual(values = c("gray", "red", "blue"))
  } else if (n_levels == 2) {
    p <- p + scale_color_manual(values = c("red", "blue"))
  }
  
  # Ligne horizontale (seuil p-value)
  if (seuil_h & P_cutoff>0) {
    ythr <- -log10(P_cutoff)
    p <- p + geom_hline(yintercept = ythr, linetype = "dashed")
  }
  
  # Lignes verticales (seuil log2FC)
  if (seuil_v & log2FC_cutoff>0) {
    xthr <- c(-log2FC_cutoff, log2FC_cutoff)
    p <- p + geom_vline(xintercept = xthr, linetype = "dashed")
  }
  
  return(p)
}
