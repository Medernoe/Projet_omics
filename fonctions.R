#==================================================================================================
#Author : Noé Méderlet 
#contact : noe.mederlet@univ-rouen.fr
#github : https://github.com/Medernoe/Projet_omics
#organism : Master Bims M2, université de rouen 
#project : Création d'une application interactive dédiée à l'analyse de données transcriptomiques,
#développée dans le cadre d'un projet universitaire du Master 2 de Bioinformatique de l'Université de Rouen.
#==================================================================================================

# ========================== DEG ========================== 

# Classe la significativité des gènes selon la pvalue et logFC
# entree : tableau avec colonne spe, un seuil pvalue et logFC 
# sortie : retourne un data modifier avec une colonne de significativité selon les seuils 
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
# entree : un tableau et des seuils 
# un vplot 
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



# ========================== Enrichissement ========================== 


# ---- ORA -----

# ORA : GO 
# Entrée : 
#   - gene_list : un vecteur de caractères contenant les IDs des gènes significatifs.
#   - label : titre de l'analyse
#   - org_db : base de données d'annotations de l'organisme (ex: org.Hs.eg.db).
#   - ontology : ontologies à tester ("BP", "CC", "MF"). # // spécifier l'ontologie depuis l'input utilisateur pour accélérer le processus
#   - p_adj : méthode d'ajustement de la p-value (défaut "BH").
#   - p_cutoff : seuil de significativité pour la p-value (défaut 0.05).
#   - key_type : format des IDs de gènes en entrée (ex: "SYMBOL" ou "ENTREZID" = valeur).
# Sortie : Une liste contenant les objets pour chaque ontologie testée.
run_ORA_go <- function(gene_list, label,
                       org_db = org.Hs.eg.db,
                       ontology = c("BP", "CC", "MF"),
                       p_adj = "BH", p_cutoff = 0.05,
                       key_type = "SYMBOL") {
  
  result <- list()
  for (GO_term in ontology) {
    message(paste("Calcul de l'enrichissement ORA GO pour :", GO_term))
    
    ego <- clusterProfiler::enrichGO(
      gene          = gene_list, 
      OrgDb         = org_db, 
      keyType       = key_type, 
      ont           = GO_term,  
      pAdjustMethod = p_adj, 
      pvalueCutoff  = p_cutoff,
      qvalueCutoff  = 1, # On force à 1 pour ne filtrer QUE sur la p-value // peut etre modifier pour filtrer sur les 2 
      readable      = TRUE
    )
    
    # Stocker dans la liste avec le nom de l'ontologie
    result[[GO_term]] <- ego
  }
  
  return(result)
}


# ORA : KEGG 
# Entrée : 
#   - gene_list : un vecteur contenant les IDs des gènes (ENTREZID pour KEGG).
#   - label : identifiant ou nom de l'analyse.
#   - organism : code organisme pour KEGG (ex: "hsa" pour humain, "mmu" pour souris). # // lier à l'input utilisateur "organisme"
#   - p_adj : méthode d'ajustement de la p-value.
#   - p_cutoff : seuil de significativité pour la p-value. 
# Sortie : Un objet `enrichResult` contenant l'enrichissement KEGG.
run_ORA_KEGG <- function(gene_list, label,
                         organism = "hsa",
                         p_adj = "BH", p_cutoff = 0.05) {
  
  message("Calcul de l'enrichissement ORA pour : KEGG")
  
  ekegg <- clusterProfiler::enrichKEGG(
    gene          = gene_list, 
    organism      = organism, 
    pAdjustMethod = p_adj, 
    pvalueCutoff  = p_cutoff,
    qvalueCutoff  = 1 # On force à 1 pour ne filtrer QUE sur la p-value
  )
  
  return(ekegg)
}


# ORA : Pathways (Reactome)
# Entrée : 
#   - gene_list : un vecteur contenant les IDs des gènes (ENTREZID pour Reactome).
#   - label : identifiant ou nom de l'analyse.
#   - organism : nom de l'organisme (ex: "human", "mouse"). # // lier à l'input utilisateur "organisme"
#   - p_adj : méthode d'ajustement de la p-value.
#   - p_cutoff : seuil de significativité pour la p-value. # // je veux seuiller sur pvalue !!
# Sortie : Un objet `enrichResult` contenant l'enrichissement Reactome.
run_ORA_pathway <- function(gene_list, label,
                            organism = "human",
                            p_adj = "BH", p_cutoff = 0.05) {
  
  message("Calcul de l'enrichissement ORA pour : Reactome Pathways")
  
  epath <- ReactomePA::enrichPathway(
    gene          = gene_list, 
    organism      = organism, 
    pAdjustMethod = p_adj, 
    pvalueCutoff  = p_cutoff,
    qvalueCutoff  = 1, # On force à 1 pour ne filtrer QUE sur la p-value
    readable      = TRUE # Permet de repasser en SYMBOL dans les résultats si possible
  )
  
  return(epath)
}



# ----- GSEA -----

# GSEA : GO 
# Entrée : 
#   - ranked_gene_list : un vecteur numérique nommé, trié par ordre décroissant (ex: Log2FC), dont les noms sont les IDs des gènes.
#   - label : identifiant ou nom de l'analyse.
#   - org_db : base de données de l'organisme.
#   - ontology : vecteur des ontologies à tester ("BP", "CC", "MF"). // spécifier l'ontologie depuis l'input utilisateur pour accélérer le processus
#   - p_adj : méthode d'ajustement.
#   - p_cutoff : seuil de significativité de la p-value ajustée.
#   - key_type : type d'ID des gènes (ex: "SYMBOL").
# Sortie : Une liste contenant les objets `gseaResult` pour chaque ontologie.
run_gsea_go <- function(ranked_gene_list, label,
                        org_db = org.Hs.eg.db,
                        ontology = c("BP", "CC", "MF"),
                        p_adj = "BH", p_cutoff = 0.05,
                        key_type = "SYMBOL") {
  
  result <- list()
  for (GO_term in ontology) {
    message(paste("Calcul GSEA GO pour :", GO_term))
    
    gse <- clusterProfiler::gseGO(
      geneList      = ranked_gene_list,
      OrgDb         = org_db,
      keyType       = key_type,
      ont           = GO_term,
      pAdjustMethod = p_adj,
      pvalueCutoff  = p_cutoff,
      verbose       = FALSE
    )
    result[[GO_term]] <- gse
  }
  return(result)
}


# GSEA : KEGG 
# Entrée : 
#   - ranked_gene_list : un vecteur numérique nommé et trié (noms = ENTREZID obligatoirement pour KEGG).
#   - label : identifiant ou nom de l'analyse.
#   - organism : code organisme ("hsa" pour humain, "mmu" pour souris). // lier à l'input utilisateur "organisme"
#   - p_adj : méthode d'ajustement.
#   - p_cutoff : seuil de significativité de la p-value. 
# Sortie : Un objet `gseaResult` contenant les résultats GSEA KEGG.
run_gsea_kegg <- function(ranked_gene_list, label,
                          organism = "hsa", 
                          p_adj = "BH", p_cutoff = 0.05) {
  
  message("Calcul GSEA KEGG")
  gse_kegg <- clusterProfiler::gseKEGG(
    geneList      = ranked_gene_list,
    organism      = organism,
    pAdjustMethod = p_adj,
    pvalueCutoff  = p_cutoff,
    verbose       = FALSE
  )
  return(gse_kegg)
}


# GSEA : Reactome 
# Entrée : 
#   - ranked_gene_list : un vecteur numérique nommé et trié (noms = ENTREZID obligatoirement pour Reactome).
#   - label : identifiant ou nom de l'analyse.
#   - organism : nom de l'organisme ("human", "mouse"). // lier à l'input utilisateur "organisme"
#   - p_adj : méthode d'ajustement.
#   - p_cutoff : seuil de significativité de la p-value.
# Sortie : Un objet `gseaResult` contenant les résultats GSEA Reactome.
run_gsea_reactome <- function(ranked_gene_list, label,
                              organism = "human", 
                              p_adj = "BH", p_cutoff = 0.05) {
  
  message("Calcul GSEA Reactome")
  gse_reac <- ReactomePA::gsePathway(
    geneList      = ranked_gene_list,
    organism      = organism,
    pAdjustMethod = p_adj,
    pvalueCutoff  = p_cutoff,
    verbose       = FALSE
  )
  return(gse_reac)
}


# ========================= Visualisation ========================== 


#  ----- ORA / GSEA -----

# Plot : Dotplot
# Entrée : 
#   - enrich_obj : Objet de résultat (ORA ou GSEA)
#   - label : Titre du graphique
#   - top_n : Nombre de catégories à afficher // Lier à un numeric/slider input au front
# Sortie : Objet ggplot
generate_dotplot <- function(enrich_obj, label = 'Enrichissement', top_n = 10) {
  p <- enrichplot::dotplot(enrich_obj, showCategory = top_n) +
    ggplot2::ggtitle(paste0("Dotplot – ", label)) +
    ggplot2::theme_minimal()
  return(p)
}


# Plot : Cnetplot (Réseau Gènes-Concepts)
# Entrée : 
#   - enrich_obj : Objet de résultat (ORA ou GSEA)
#   - label : Titre du graphique 
#   - top_n : Nombre de concepts à afficher // Lier à un numeric input (souvent plus bas, ex: 5)
# Sortie : Objet ggplot
generate_cnetplot <- function(enrich_obj, label = 'Enrichissement', top_n = 5) {
  p <- enrichplot::cnetplot(enrich_obj, 
                            showCategory = top_n, 
                            circular = FALSE, 
                            colorEdge = TRUE) +
    ggplot2::ggtitle(paste0("Réseau Gènes-Concepts – ", label)) +
    ggplot2::theme_minimal()
  return(p)
}


# Plot : Emapplot (Carte de similarité)
# Entrée : 
#   - enrich_obj : Objet de résultat (ORA ou GSEA)
#   - label : Titre du graphique
#   - top_n : Nombre de catégories à afficher // Lier à un numeric/slider input
# Sortie : Objet ggplot
generate_emapplot <- function(enrich_obj, label = 'Enrichissement', top_n = 10) {
  # L'Emapplot nécessite d'abord le calcul de la similarité des termes
  sim_obj <- enrichplot::pairwise_termsim(enrich_obj)
  
  p <- enrichplot::emapplot(sim_obj, showCategory = top_n) +
    ggplot2::ggtitle(paste0("Carte de similarité (Emap) – ", label)) +
    ggplot2::theme_minimal()
  return(p)
}


# Plot : Upsetplot (Intersections des gènes)
# Entrée : 
#   - enrich_obj : Objet de résultat (ORA ou GSEA)
#   - label : Titre du graphique
# Sortie : Objet ggplot / upset
generate_upsetplot <- function(enrich_obj, label = 'Enrichissement') {
  p <- enrichplot::upsetplot(enrich_obj) +
    ggplot2::ggtitle(paste0("Intersections des gènes – ", label)) +
    ggplot2::theme_minimal()
  return(p)
}


# Plot : Heatplot (Heatmap Gènes-Termes)
# Entrée : 
#   - enrich_obj : Objet de résultat (ORA ou GSEA)
#   - label : Titre du graphique
#   - top_n : Nombre de catégories à afficher // Lier à un numeric/slider input
# Sortie : Objet ggplot
generate_heatplot <- function(enrich_obj, label = 'Enrichissement', top_n = 10) {
  p <- enrichplot::heatplot(enrich_obj, showCategory = top_n) +
    ggplot2::ggtitle(paste0("Heatmap Gènes-Termes – ", label)) +
    ggplot2::theme_minimal()
  return(p)
}


# ----- ORA specific ----- 
# Plot : Barplot
# Entrée : 
#   - enrich_obj : Objet de résultat (Uniquement ORA) // cacher/ ne pas afficher ce plot au front si l'utilisateur choisit GSEA
#   - label : Titre du graphique
#   - top_n : Nombre de catégories à afficher // Lier à un numeric/slider input
# Sortie : Objet ggplot
generate_barplot <- function(enrich_obj, label = 'Enrichissement', top_n = 10) {
  p <- barplot(enrich_obj, showCategory = top_n) +
    ggplot2::ggtitle(paste0("Barplot – ", label)) +
    ggplot2::theme_minimal()
  return(p)
}


# Plot : Goplot
# Entrée : 
#   - enrich_obj : Objet de résultat (Uniquement ORA) // cacher/ ne pas afficher ce plot au front si l'utilisateur choisit GSEA
#   - label : Titre du graphique
#   - top_n : Nombre de catégories à afficher // Lier à un numeric/slider input
# Sortie : Objet ggplot
generate_goplot <- function(enrich_obj, label = 'Enrichissement', top_n = 10) {
  p <- goplot(enrich_obj, showCategory = top_n) +
    ggplot2::ggtitle(paste0("Barplot – ", label)) +
    ggplot2::theme_minimal()
  return(p)
}



# ----- GSEA specific ----- 
# Plot : Ridgeplot
# Entrée : 
#   - gse_obj : Objet de résultat (Uniquement GSEA) // Afficher uniquement si l'utilisateur choisit GSEA
#   - label : Titre du graphique
#   - top_n : Nombre de catégories à afficher
# Sortie : Objet ggplot
generate_ridgeplot <- function(gse_obj, label = 'GSEA', top_n = 10) {
  p <- enrichplot::ridgeplot(gse_obj, showCategory = top_n) +
    ggplot2::ggtitle(paste0("Ridgeplot (Distribution log2FC) – ", label)) +
    ggplot2::theme_minimal()
  return(p)
}


# Plot : GSEA Plot (Type 2 - Multiples pathways)
# Entrée : 
#   - gse_obj : Objet de résultat (Uniquement GSEA) // Afficher uniquement si l'utilisateur choisit GSEA
#   - gene_set_ids : Vecteur d'indices des pathways à afficher (ex: 1:3 pour les 3 premiers) // Lier à une sélection multiple au front
#   - show_pvalue : Afficher ou non la table des p-values // Lier à une checkbox au front
# Sortie : Objet ggplot complexe
generate_gseaplot2 <- function(gse_obj, gene_set_ids = 1:3, show_pvalue = TRUE) {
  p <- enrichplot::gseaplot2(gse_obj, 
                             geneSetID = gene_set_ids, 
                             pvalue_table = show_pvalue)
  return(p)
}


# Plot : GSEA Rank
# Entrée : 
#   - gse_obj : Objet de résultat (Uniquement GSEA) // Afficher uniquement si l'utilisateur choisit GSEA
#   - gene_set_id : Index du pathway spécifique à afficher (un seul) // Lier à un input numérique ou un menu déroulant
# Sortie : Objet ggplot
generate_gsearank <- function(gse_obj, gene_set_id = 1) {
  # Récupère le nom dynamique du pathway pour le titre
  pathway_title <- gse_obj[gene_set_id, "Description"]
  
  p <- enrichplot::gsearank(gse_obj, geneSetID = gene_set_id, title = pathway_title)
  return(p)
}




# ----- Custom MAnhattan plot ----- 

# Plot : Manhattan Plot personnalisé
# Entrée : 
#   - enrich_list : Une liste nommée d'objets d'enrichissement (ex: list("GO-BP" = res$BP, "KEGG" = res_kegg)) // Le front doit assembler les résultats sélectionnés dans une liste nommée avant d'appeler la fonction
#   - label : Titre du graphique
#   - p_cutoff : Seuil pour tracer la ligne de significativité # // Lier à l'input utilisateur "p_cutoff"
#   - cap_y : Valeur maximale pour tronquer l'axe Y (ex: 15) pour éviter que les p-values extrêmes n'écrasent le plot // Optionnel, lier à un numeric input "Y max" ou laisser NULL
# Sortie : Objet ggplot
generate_manhattan_plot <- function(enrich_list, label = "Enrichissement Global", p_cutoff = 0.05, cap_y = NULL) {
  
  # 1. Extraction et combinaison des données de la liste
  df_list <- lapply(names(enrich_list), function(source_name) {
    obj <- enrich_list[[source_name]]
    if (is.null(obj)) return(NULL)
    
    # Convertir l'objet endata.frame classique
    df <- as.data.frame(obj)
    if (nrow(df) == 0) return(NULL)
    
    df$Source <- source_name
    df$logP <- -log10(df$p.adjust)
    
    # Gérer la taille des points (Count pour ORA, setSize pour GSEA)
    if ("Count" %in% colnames(df)) {
      df$Size <- as.numeric(df$Count)
    } else if ("setSize" %in% colnames(df)) {
      df$Size <- as.numeric(df$setSize)
    } else {
      df$Size <- 1 # Sécurité
    }
    
    return(df)
  })
  
  # Fusionner la liste en un seul data.frame
  plot_data <- dplyr::bind_rows(df_list)
  
  if (nrow(plot_data) == 0) {
    warning("Aucun résultat significatif à afficher pour le Manhattan Plot.")
    return(ggplot2::ggplot() + ggplot2::ggtitle(paste0("Aucun résultat - ", label)) + ggplot2::theme_void())
  }
  
  # 2. Préparation des axes (Création de l'index X artificiel)
  plot_data <- plot_data %>%
    dplyr::arrange(Source, dplyr::desc(logP)) %>%
    dplyr::mutate(Index = dplyr::row_number())
  
  # Plafonner (cap) les valeurs Y extrêmes si demandé par l'utilisateur
  if (!is.null(cap_y)) {
    plot_data$logP <- ifelse(plot_data$logP > cap_y, cap_y, plot_data$logP)
  }
  
  # Calculer le centre de chaque groupe pour placer les étiquettes sur l'axe X
  axis_data <- plot_data %>%
    dplyr::group_by(Source) %>%
    dplyr::summarize(Center = mean(Index), N = dplyr::n(), .groups = 'drop') %>%
    dplyr::mutate(Label = paste0(Source, "\n(", N, ")"))
  
  # Calculer la position de la ligne de seuil
  sig_line <- -log10(p_cutoff)
  
  # 3. Création du graphique avec ggplot2
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Index, y = logP, color = Source, size = Size)) +
    ggplot2::geom_point(alpha = 0.8) +
    
    # Ligne de seuil de significativité
    ggplot2::geom_hline(yintercept = sig_line, linetype = "dashed", color = "grey50") +
    
    # Personnalisation des axes
    ggplot2::scale_x_continuous(breaks = axis_data$Center, labels = axis_data$Label) +
    ggplot2::scale_size_continuous(range = c(2, 6), guide = "none") + # Empêche la légende de taille de surcharger le plot
    
    # Labels et thèmes
    ggplot2::labs(
      title = paste0("Manhattan Plot – ", label),
      x = "",
      y = expression("-log"[10]*"(p.adjust)"),
      color = "Base de données"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 1, face = "bold"),
      panel.grid.major.x = ggplot2::element_blank(), # Retire la grille verticale pour un effet GWAS
      panel.grid.minor.x = ggplot2::element_blank(),
      legend.position = "right"
    )
  
  # Ajouter une note visuelle si les valeurs ont été tronquées (cappées)
  if (!is.null(cap_y)) {
    p <- p + ggplot2::annotate("text", x = max(plot_data$Index), y = cap_y + 0.2,
                               label = "Valeurs plafonnées", hjust = 1, size = 3.5, color = "grey30")
  }
  
  return(p)
}