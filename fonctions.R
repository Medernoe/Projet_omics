#==================================================================================================
#Author : Noé Méderlet 
#contact : noe.mederlet@gmail.com
#github : 
#organism : Master Bims M2, université de rouen 
#project : 


#==================================================================================================


#fonction
significativity <- function(data, log2FC_cutoff, P_cutoff){ 
  data$Significance <- "Not significant"
  data$Significance[data$p_value < P_cutoff & data$log2FC >  log2FC_cutoff] <- "Upregulated"
  data$Significance[data$p_value < P_cutoff & data$log2FC < -log2FC_cutoff] <- "Downregulated"
  data$Significance <- factor(data$Significance, levels=c("Not significant","Upregulated","Downregulated"))
  
  return(data)
}


plot_volcano <- function(data, log2FC_cutoff, P_cutoff,
                         title = "Volcano Plot",
                         seuil_v = FALSE, seuil_h = FALSE) {
  
  # Base du plot
  p <- ggplot(data, aes(x = log2FC, y = -log10(p_value), color = Significance)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c("gray", "red", "blue")) +
    labs(title = title,
         x = "Log2 Fold Change",
         y = "-log10(adj P-value)",
         color = "Significance") +
    theme_minimal()
  
  # Ajout conditionnel de la ligne horizontale
  if (seuil_h) {
    ythr <- -log10(P_cutoff)
    p <- p + geom_hline(yintercept = ythr, linetype = "dashed")
  }
  
  # Ajout conditionnel de la ligne verticale
  if (seuil_v) {
    xthr <- c(-log2FC_cutoff, log2FC_cutoff)
    p <- p + geom_vline(xintercept = xthr, linetype = "dashed")
  }
  
  return(p)
}

