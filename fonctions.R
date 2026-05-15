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

# ---- Convertion ID ----

# Pour ORA
convert_symbols_to_entrez <- function(symbols, org_db) {
  converted <- tryCatch(
    clusterProfiler::bitr(symbols,
                          fromType = "SYMBOL",
                          toType   = "ENTREZID",
                          OrgDb    = org_db,
                          drop     = TRUE),
    error = function(e) NULL
  )
  if (is.null(converted) || nrow(converted) == 0) return(NULL)
  converted$ENTREZID
}

# Pour GSEA - Conserve le ranking
convert_ranked_to_entrez <- function(ranked_gene_list, org_db) {
  symbols <- names(ranked_gene_list)
  mapping <- tryCatch(
    clusterProfiler::bitr(symbols,
                          fromType = "SYMBOL",
                          toType   = "ENTREZID",
                          OrgDb    = org_db,
                          drop     = TRUE),
    error = function(e) NULL
  )
  if (is.null(mapping) || nrow(mapping) == 0) return(NULL)
  
  # On crée un nouveau vecteur avec les ENTREZID comme noms
  ranked_entrez <- ranked_gene_list[mapping$SYMBOL]
  names(ranked_entrez) <- mapping$ENTREZID
  
  # Re-tri et déduplication (au cas où plusieurs SYMBOL → même ENTREZID)
  ranked_entrez <- ranked_entrez[!duplicated(names(ranked_entrez))]
  sort(ranked_entrez, decreasing = TRUE)
}


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
                       p_adj = "BH",
                       p_cutoff = 0.05,
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
                            showCategory = top_n) +
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
generate_upsetplot <- function(enrich_obj, label = 'Enrichissement', top_n = 10) {
  # n contrôle le nombre de catégories (termes GO) affichées dans l'upsetplot.
  # Sans ce paramètre, enrichplot utilise une valeur par défaut indépendante du slider UI.
  p <- enrichplot::upsetplot(enrich_obj, n = top_n) +
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

# ----- Custom MAnhattan plot ----- 

# Plot : Manhattan Plot personnalisé
# Entrée : 
#   - enrich_list : Une liste nommée d'objets d'enrichissement (ex: list("GO-BP" = res$BP, "KEGG" = res_kegg)) // Le front doit assembler les résultats sélectionnés dans une liste nommée avant d'appeler la fonction
#   - label : Titre du graphique
#   - p_cutoff : Seuil pour tracer la ligne de significativité # // Lier à l'input utilisateur "p_cutoff"
#   - cap_y : Valeur maximale pour tronquer l'axe Y (ex: 15) pour éviter que les p-values extrêmes n'écrasent le plot // Optionnel, lier à un numeric input "Y max" ou laisser NULL
# Sortie : Objet ggplot
generate_manhattan_plot <- function(go_results,
                                    label        = "GO Enrichissement",
                                    p_cutoff     = 0.05,
                                    top_n_labels = 8L,
                                    cap_y        = NULL) {
  
  # 1) Extraction : convertir chaque objet GO en data.frame uniforme
  .extract_go <- function(obj, ontology_name) {
    if (is.null(obj)) return(NULL)
    df <- tryCatch(as.data.frame(obj), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0L) return(NULL)
    
    # p.adjust en priorité, sinon pvalue brute
    p_col <- intersect(c("p.adjust", "pvalue"), colnames(df))[1]
    if (is.na(p_col)) return(NULL)
    
    # Taille du gene set : Count (ORA) ou setSize (GSEA)
    size_col <- if ("Count" %in% colnames(df)) {
      df$Count
    } else if ("setSize" %in% colnames(df)) {
      df$setSize
    } else {
      rep(1L, nrow(df))
    }
    
    data.frame(
      term_id    = df$ID,
      term_name  = df$Description,
      ontology   = ontology_name,
      p_adj      = as.numeric(df[[p_col]]),
      gene_count = as.integer(size_col),
      stringsAsFactors = FALSE
    )
  }
  
  # 2) Fusion des ontologies disponibles
  df_raw <- do.call(rbind, lapply(names(go_results), function(ont) {
    .extract_go(go_results[[ont]], ont)
  }))
  
  if (is.null(df_raw) || nrow(df_raw) == 0L) return(NULL)
  
  # 3) Nettoyage
  df_raw <- df_raw[!is.na(df_raw$p_adj) & is.finite(df_raw$p_adj) & df_raw$p_adj > 0, ]
  if (nrow(df_raw) == 0L) return(NULL)
  
  # 4) Force l'ordre BP -> CC -> MF (toujours dans cet ordre, même si certaines manquent)
  ontology_order <- c("BP", "CC", "MF")
  df_raw$ontology <- factor(df_raw$ontology, 
                            levels = intersect(ontology_order, unique(df_raw$ontology)))
  
  # 5) Calcul du -log10(p) avec cap éventuel
  df_raw$neg_log_p <- -log10(df_raw$p_adj)
  if (!is.null(cap_y)) {
    df_raw$neg_log_p <- pmin(df_raw$neg_log_p, cap_y)
  }
  
  # 6) Position x : on attribue une position au sein de chaque ontologie
  df_raw <- df_raw[order(df_raw$ontology, df_raw$p_adj), ]
  df_raw$x_pos <- seq_len(nrow(df_raw))
  
  # 7) Top termes à annoter
  top_terms <- df_raw[df_raw$p_adj < p_cutoff, ]
  top_terms <- top_terms[order(top_terms$p_adj), ]
  top_terms <- head(top_terms, top_n_labels)
  
  # 8) Centre de chaque ontologie pour positionner les labels d'axe x
  cat_centers <- aggregate(x_pos ~ ontology, data = df_raw, FUN = mean)
  
  # 9) Palette fixe pour les 3 ontologies (cohérence visuelle)
  ontology_colors <- c("BP" = "#1A5FA8",   # bleu
                       "CC" = "#3B6D11",   # vert
                       "MF" = "#7A4200")   # orange
  
  # 10) Labels lisibles pour la légende
  ontology_labels <- c("BP" = "Processus Biologique",
                       "CC" = "Composant Cellulaire",
                       "MF" = "Fonction Moléculaire")
  
  # 11) Construction du ggplot
  p <- ggplot(df_raw, aes(x = x_pos, y = neg_log_p, 
                          color = ontology,
                          size = gene_count,
                          text = paste0("Terme : ", term_name,
                                        "<br>ID : ", term_id,
                                        "<br>Ontologie : ", ontology,
                                        "<br>p.adjust : ", signif(p_adj, 3),
                                        "<br>Gènes : ", gene_count))) +
    geom_point(alpha = 0.7) +
    geom_hline(yintercept = -log10(p_cutoff), 
               linetype = "dashed", color = "red", linewidth = 0.5) +
    scale_color_manual(values = ontology_colors, 
                       labels = ontology_labels,
                       name = "Ontologie GO") +
    scale_size_continuous(range = c(1.5, 6), name = "Nb gènes") +
    scale_x_continuous(breaks = cat_centers$x_pos, 
                       labels = cat_centers$ontology) +
    labs(
      title = paste0("Manhattan plot — ", label),
      x = "Ontologie GO",
      y = "-log10(p.adjust)"
    ) +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.x        = element_text(size = 11, face = "bold")
    )
  
  return(p)
}