#===========================Presentation========================================
# Authors : Noé Méderlet, Mehdi Tachekort, Mathieu Cartier, Valentin Fourdreane
# contact : noe.mederlet@univ-rouen.fr, mehdi.tachekort@univ-rouen.fr, 
#           mathieu.cartier@univ-rouen.fr, valentin.fourdraine@univ-rouen.fr
# github : https://github.com/Medernoe/Projet_omics
# organism : Master Bims M2, université de rouen 
# project : Création d'une application interactive dédiée à l'analyse de données 
#           transcriptomiques, développée dans le cadre d'un projet universitaire 
#           du Master 2 de Bioinformatique de l'Université de Rouen.
#==============================================================================#

#=========Chargement du global.R===============================================
source("global.R")

#=========================SERVEUR==============================================

function(input, output, session) {
  
  ####============================INITIALISATION================================
  
  # Cache le loader hamster une fois l'app prête
  waiter::waiter_hide()
  
  
  ####============================CHARGEMENT DES DONNÉES========================
  
  # Colonnes minimales requises dans le CSV uploadé
  required_columns <- c("GeneName", "log2FC", "pval")
  
  #####=======================Lecture du fichier CSV============================
  raw_data <- reactive({
    req(input$deg_file)
    
    ext <- tools::file_ext(input$deg_file$name)
    if (tolower(ext) != "csv") {
      shinyalert(
        title = "Format non valide",
        text  = "Veuillez importer un fichier .csv uniquement.",
        type  = "error"
      )
      return(NULL)
    }
    
    df <- tryCatch(
      read.csv(input$deg_file$datapath, sep = ";"),
      error = function(e) {
        shinyalert(
          title = "Erreur de lecture",
          text  = "Impossible de lire le fichier. Vérifiez qu'il s'agit d'un CSV valide.",
          type  = "error"
        )
        return(NULL)
      }
    )
    
    if (is.null(df)) return(NULL)
    
    if (!all(required_columns %in% colnames(df))) {
      shinyalert(
        title = "Colonnes manquantes",
        text  = paste0(
          "Le fichier doit contenir les colonnes suivantes : ",
          paste(required_columns, collapse = ", ")
        ),
        type = "error"
      )
      return(NULL)
    }
    
    df
  })
  
  #####=======================Debounce des sliders==============================
  fc_threshold_debounced     <- debounce(reactive(input$fc_threshold), 300)
  pvalue_threshold_debounced <- debounce(reactive(input$pvalue_threshold), 300)
  
  #####=======================Classification des gènes==========================
  processed_data <- reactive({
    if (is.null(raw_data())) return(NULL)
    
    significativity(
      data          = raw_data(),
      log2FC_cutoff = fc_threshold_debounced(),
      P_cutoff      = pvalue_threshold_debounced()
    )
  })
  
  #####=========================MAPPING ESPÈCE -> BASES=========================
  
  # OrgDb pour GO (clusterProfiler::enrichGO/gseGO)
  species_to_orgdb <- list(
    "Homo sapiens"            = org.Hs.eg.db::org.Hs.eg.db,
    "Mus musculus"            = org.Mm.eg.db::org.Mm.eg.db,
    "Drosophila melanogaster" = org.Dm.eg.db::org.Dm.eg.db
  )
  
  # Codes KEGG (clusterProfiler::enrichKEGG/gseKEGG)
  species_to_kegg <- c(
    "Homo sapiens"            = "hsa",
    "Mus musculus"            = "mmu",
    "Drosophila melanogaster" = "dme"
  )
  
  # Codes Reactome (ReactomePA::enrichPathway/gsePathway)
  species_to_reactome <- c(
    "Homo sapiens"            = "human",
    "Mus musculus"            = "mouse",
    "Drosophila melanogaster" = "fly"
  )
  
  #####=======================Sélection des gènes par direction=================
  # Filtre les gènes selon la direction (up / down / both)
  get_directional_genes <- function(df, direction) {
    df_sig <- df[as.character(df$Significance) != "Not significant", ]
    
    if (direction == "up") {
      return(df_sig$GeneName[as.character(df_sig$Significance) == "Upregulated"])
    } else if (direction == "down") {
      return(df_sig$GeneName[as.character(df_sig$Significance) == "Downregulated"])
    }
    df_sig$GeneName
  }
  
  
  ####============================ONGLET DEG====================================
  
  #####=======================Volcano plot======================================
  
  create_volcano <- reactive({
    req(processed_data())
    
    selected_row <- input$deg_table_rows_selected
    
    w <- make_waiter(
      id = "volcano_plot",
      message = "Calcul Pathway GSEA en cours... (plusieurs minutes possibles)"
    )
    w$show()
    on.exit(w$hide())
    
    plot_volcano(
      data          = processed_data(),
      log2FC_cutoff = fc_threshold_debounced(),
      P_cutoff      = pvalue_threshold_debounced(),
      seuil_v       = input$show_vline,
      seuil_h       = input$show_hline,
      title         = input$deg_title,
      highlight_row = selected_row
    )
  })
  
  output$show_volcano_error <- reactive({
    is.null(input$deg_file) || is.null(raw_data())
  })
  outputOptions(output, "show_volcano_error", suspendWhenHidden = FALSE)
  
  output$volcano_error_img <- renderImage({
    list(
      src         = "www/erreur_format.jpg",
      contentType = "image/jpeg",
      width       = "70%",
      height      = "auto",
      alt         = "Format de fichier attendu"
    )
  }, deleteFile = FALSE)
  
  output$volcano_plot <- renderPlotly({
    req(input$deg_file)
    
    ggplotly(create_volcano(), tooltip = c("x", "y", "colour")) %>%
      layout(dragmode = "zoom", hovermode = "closest") %>%
      config(
        displayModeBar         = input$deg_toolbox,
        modeBarButtonsToAdd    = list("drawrect", "eraseshape"),
        modeBarButtonsToRemove = list("toImage"),
        displaylogo            = FALSE
      ) %>%
      plotly::toWebGL()
  })
  
  output$downloadVolcano <- downloadHandler(
    filename = function() paste0("Volcano_plot_", Sys.Date(), ".png"),
    content = function(file) {
      ggsave(file, plot = create_volcano(), width = 12, height = 8, dpi = 300)
    }
  )
  
  #####=======================Tableau des données===============================
  
  output$show_table_error <- reactive({
    is.null(input$deg_file) || is.null(raw_data())
  })
  outputOptions(output, "show_table_error", suspendWhenHidden = FALSE)
  
  output$table_error_text <- renderText({
    "Veuillez charger un fichier CSV au format attendu pour explorer les données"
  })
  
  output$deg_table <- renderDT({
    req(processed_data())
    
    datatable(
      processed_data(),
      selection = "single",
      options = list(
        pageLength    = 10,
        scrollX       = TRUE,
        deferRender   = TRUE,
        scroller      = TRUE
      )
    )
  })
  
  #####=======================Préparation du ranking GSEA======================
  # Construit le vecteur ranked utilisé par toutes les analyses GSEA
  gsea_ranked <- reactive({
    req(processed_data())
    
    df <- processed_data()
    df <- df[!is.na(df$log2FC) & !is.na(df$GeneName), ]
    df <- df[order(-abs(df$log2FC)), ]
    df <- df[!duplicated(df$GeneName), ]
    
    ranked <- df$log2FC
    names(ranked) <- df$GeneName
    sort(ranked, decreasing = TRUE)
  })
  
  
  ####============================ONGLET GO ORA=================================
  
  #####=======================Calcul GO ORA====================================
  
  go_ora_results <- eventReactive(input$run_go_ora, {
    
    req(processed_data())
    req(length(input$go_ora_ontology) > 0)
    
    # Loader sur le plot pendant le calcul
    w <- make_waiter(
      id = "go_ora_plot",
      message = paste0("Calcul GO ORA en cours sur ", 
                       length(input$go_ora_ontology), " ontologie(s)...")
    )
    w$show()
    on.exit(w$hide())  # Garantit que le loader se cache même en cas d'erreur
    
    df <- processed_data()
    significant_genes <- get_directional_genes(df, input$go_ora_direction)
    
    if (length(significant_genes) == 0) {
      shinyalert(
        title = "Aucun gène significatif",
        text  = paste0("Aucun gène ne correspond à la direction sélectionnée (",
                       input$go_ora_direction, "). Ajustez les seuils ou la direction."),
        type  = "warning"
      )
      return(NULL)
    }
    
    org_db <- species_to_orgdb[[input$species]]
    if (is.null(org_db)) {
      shinyalert(title = "Espèce non supportée",
                 text  = paste0("Aucune base GO pour : ", input$species),
                 type  = "error")
      return(NULL)
    }
    
    showNotification(
      paste0("Calcul GO ORA en cours sur ", length(input$go_ora_ontology), " ontologie(s)..."),
      duration = NULL, id = "go_ora_running", type = "message"
    )
    
    ego_list <- tryCatch(
      run_ORA_go(
        gene_list = significant_genes,
        label     = "GO ORA",
        org_db    = org_db,
        ontology  = input$go_ora_ontology,
        p_adj     = "BH",
        p_cutoff  = 0.05,
        key_type  = "SYMBOL"
      ),
      error = function(e) NULL
    )
    
    removeNotification("go_ora_running")
    
    if (is.null(ego_list)) {
      shinyalert(title = "Erreur de calcul", type = "error")
      return(NULL)
    }
    
    valid <- !sapply(ego_list, is.null)
    if (!any(valid)) {
      shinyalert(
        title = "Aucun gène mappé",
        text  = paste0("Vérifiez que l'espèce '", input$species, 
                       "' correspond bien à vos identifiants."),
        type  = "error"
      )
      return(NULL)
    }
    
    if (!all(valid)) {
      showNotification(
        paste0("Aucun résultat pour : ", paste(names(valid)[!valid], collapse = ", ")),
        duration = 6, type = "warning"
      )
    }
    
    showNotification("Calcul GO ORA terminé", duration = 4, type = "message")
    ego_list[valid]
  })
  
  #####=======================Update sélecteur ontologie GO ORA===============
  
  observeEvent(go_ora_results(), {
    req(go_ora_results())
    available <- names(go_ora_results())
    
    labels_map <- c("BP" = "Processus Biologique (BP)",
                    "CC" = "Composant Cellulaire (CC)",
                    "MF" = "Fonction Moléculaire (MF)")
    
    choices <- setNames(available, labels_map[available])
    
    updateSelectInput(session, "go_ora_displayed_ontology",
                      choices = choices, selected = available[1])
  })
  
  #####=======================Plot GO ORA======================================
  
  output$go_ora_plot <- renderPlotly({
    
    if (is.null(input$deg_file) || is.null(raw_data()) || input$run_go_ora == 0) {
      return(
        plotly_empty(type = "scatter", mode = "markers") %>%
          layout(
            xaxis = list(visible = FALSE),
            yaxis = list(visible = FALSE),
            annotations = list(list(
              text = "Sélectionnez la ou les ontologie(s) puis cliquez sur « Lancer GO ORA »",
              showarrow = FALSE,
              font = list(size = 16, color = "#666")
            ))
          )
      )
    }
    
    req(go_ora_results())
    req(input$go_ora_displayed_ontology)
    
    ego <- go_ora_results()[[input$go_ora_displayed_ontology]]
    if (is.null(ego)) return(NULL)
    
    label <- if (is.null(input$go_ora_title) || input$go_ora_title == "") {
      paste0("GO ORA - ", input$go_ora_displayed_ontology)
    } else input$go_ora_title
    
    top_n <- input$go_ora_top_n_terms
    
    selected_plot <- switch(input$go_ora_selected_plot,
                            "Dotplot"   = generate_dotplot(ego, label = label, top_n = top_n),
                            "Barplot"   = generate_barplot(ego, label = label, top_n = top_n),
                            "Cnetplot"  = generate_cnetplot(ego, label = label, top_n = top_n),
                            "Emapplot"  = generate_emapplot(ego, label = label, top_n = top_n),
                            "Goplot"    = generate_goplot(ego, label = label, top_n = top_n),
                            "Upsetplot" = generate_upsetplot(ego, label = label),
                            "Heatplot"  = generate_heatplot(ego, label = label, top_n = top_n),
                            "Manhattan" = generate_manhattan_plot(go_results = go_ora_results(),label = label,
                              p_cutoff     = 0.05,
                              top_n_labels = min(top_n, 15))
    )
    
    if (is.null(selected_plot)) return(NULL)
    
    ggplotly(selected_plot, tooltip = c("x", "y", "text")) %>%
      layout(dragmode = "zoom", hovermode = "closest") %>%
      config(displayModeBar = input$go_ora_toolbox,
             modeBarButtonsToRemove = list("toImage"),
             displaylogo = FALSE) %>%
      plotly::toWebGL()
  })
  
  output$go_ora_status <- renderText({
    if (input$run_go_ora == 0) return("")
    results <- go_ora_results()
    if (is.null(results)) return("Aucun résultat disponible")
    paste0("Calcul terminé (", length(results), " ontologie(s))")
  })
  
  output$downloadGoOra <- downloadHandler(
    filename = function() {
      paste0("GO_ORA_", input$go_ora_selected_plot, "_",
             input$go_ora_displayed_ontology, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      ego <- go_ora_results()[[input$go_ora_displayed_ontology]]
      label <- if (is.null(input$go_ora_title) || input$go_ora_title == "") {
        paste0("GO ORA - ", input$go_ora_displayed_ontology)
      } else input$go_ora_title
      top_n <- input$go_ora_top_n_terms
      
      selected_plot <- switch(input$go_ora_selected_plot,
                              "Dotplot"   = generate_dotplot(ego, label = label, top_n = top_n),
                              "Barplot"   = generate_barplot(ego, label = label, top_n = top_n),
                              "Cnetplot"  = generate_cnetplot(ego, label = label, top_n = top_n),
                              "Emapplot"  = generate_emapplot(ego, label = label, top_n = top_n),
                              "Goplot"    = generate_goplot(ego, label = label, top_n = top_n),
                              "Upsetplot" = generate_upsetplot(ego, label = label),
                              "Heatplot"  = generate_heatplot(ego, label = label, top_n = top_n)
      )
      ggsave(file, plot = selected_plot, width = 12, height = 8, dpi = 300)
    }
  )
  
  
  ####============================ONGLET GO GSEA================================
  
  #####=======================Calcul GO GSEA===================================
  
  go_gsea_results <- eventReactive(input$run_go_gsea, {
    
    req(processed_data())
    req(length(input$go_gsea_ontology) > 0)
    
    # Loader
    w <- make_waiter(
      id = "go_gsea_plot",
      message = "Calcul GO GSEA en cours... (plusieurs minutes possibles)"
    )
    w$show()
    on.exit(w$hide())
    
    ranked <- gsea_ranked()
    if (length(ranked) == 0) {
      shinyalert(title = "Aucun gène valide", type = "warning")
      return(NULL)
    }
    
    org_db <- species_to_orgdb[[input$species]]
    if (is.null(org_db)) {
      shinyalert(title = "Espèce non supportée", type = "error")
      return(NULL)
    }
    
    showNotification(
      paste0("Calcul GO GSEA en cours sur ", length(input$go_gsea_ontology), 
             " ontologie(s)... Peut prendre plusieurs minutes."),
      duration = NULL, id = "go_gsea_running", type = "message"
    )
    
    gse_list <- tryCatch(
      run_gsea_go(
        ranked_gene_list = ranked,
        label            = "GO GSEA",
        org_db           = org_db,
        ontology         = input$go_gsea_ontology,
        p_adj            = "BH",
        p_cutoff         = 0.05,
        key_type         = "SYMBOL"
      ),
      error = function(e) NULL
    )
    
    removeNotification("go_gsea_running")
    
    if (is.null(gse_list)) {
      shinyalert(title = "Erreur de calcul", type = "error")
      return(NULL)
    }
    
    valid <- sapply(gse_list, function(g) !is.null(g) && nrow(as.data.frame(g)) > 0)
    if (!any(valid)) {
      shinyalert(title = "Aucun terme enrichi", type = "warning")
      return(NULL)
    }
    
    if (!all(valid)) {
      showNotification(
        paste0("Aucun résultat pour : ", paste(names(valid)[!valid], collapse = ", ")),
        duration = 6, type = "warning"
      )
    }
    
    showNotification("Calcul GO GSEA terminé", duration = 4, type = "message")
    gse_list[valid]
  })
  
  #####=======================Update sélecteur ontologie GO GSEA==============
  
  observeEvent(go_gsea_results(), {
    req(go_gsea_results())
    available <- names(go_gsea_results())
    
    labels_map <- c("BP" = "Processus Biologique (BP)",
                    "CC" = "Composant Cellulaire (CC)",
                    "MF" = "Fonction Moléculaire (MF)")
    
    choices <- setNames(available, labels_map[available])
    
    updateSelectInput(session, "go_gsea_displayed_ontology",
                      choices = choices, selected = available[1])
  })
  
  #####=======================Filtrage direction GSEA (post-calcul)===========
  
  # Applique le filtre NES (up / down / both) sur l'objet gseaResult
  filter_gsea_by_direction <- function(gse, direction) {
    if (is.null(gse) || nrow(as.data.frame(gse)) == 0) return(gse)
    if (direction == "both") return(gse)
    
    df <- gse@result
    keep <- if (direction == "up") df$NES > 0 else df$NES < 0
    
    if (!any(keep)) return(NULL)
    
    gse@result <- df[keep, , drop = FALSE]
    gse
  }
  
  #####=======================Plot GO GSEA====================================
  
  output$go_gsea_plot <- renderPlotly({
    
    if (is.null(input$deg_file) || is.null(raw_data()) || input$run_go_gsea == 0) {
      return(
        plotly_empty(type = "scatter", mode = "markers") %>%
          layout(
            xaxis = list(visible = FALSE),
            yaxis = list(visible = FALSE),
            annotations = list(list(
              text = "Sélectionnez la ou les ontologie(s) puis cliquez sur « Lancer GO GSEA »",
              showarrow = FALSE,
              font = list(size = 16, color = "#666")
            ))
          )
      )
    }
    
    req(go_gsea_results())
    req(input$go_gsea_displayed_ontology)
    
    gse <- go_gsea_results()[[input$go_gsea_displayed_ontology]]
    if (is.null(gse)) return(NULL)
    
    # Filtrage par direction NES
    gse <- filter_gsea_by_direction(gse, input$go_gsea_direction)
    if (is.null(gse) || nrow(as.data.frame(gse)) == 0) {
      return(
        plotly_empty() %>%
          layout(annotations = list(list(
            text = "Aucun pathway dans cette direction",
            showarrow = FALSE, font = list(size = 16, color = "#666")
          )))
      )
    }
    
    label <- if (is.null(input$go_gsea_title) || input$go_gsea_title == "") {
      paste0("GO GSEA - ", input$go_gsea_displayed_ontology)
    } else input$go_gsea_title
    
    top_n <- input$go_gsea_top_n_terms
    
    selected_plot <- switch(input$go_gsea_selected_plot,
                            "Dotplot"   = generate_dotplot(gse, label = label, top_n = top_n),
                            "Cnetplot"  = generate_cnetplot(gse, label = label, top_n = top_n),
                            "Emapplot"  = generate_emapplot(gse, label = label, top_n = top_n),
                            "Upsetplot" = generate_upsetplot(gse, label = label),
                            "Heatplot"  = generate_heatplot(gse, label = label, top_n = top_n),
                            "Ridgeplot" = generate_ridgeplot(gse, label = label, top_n = top_n),
                            "GSEAplot2" = generate_gseaplot2(gse, gene_set_ids = 1:min(3, nrow(as.data.frame(gse)))),
                            "GSEArank"  = generate_gsearank(gse, gene_set_id = 1)
    )
    
    if (is.null(selected_plot)) return(NULL)
    
    ggplotly(selected_plot, tooltip = c("x", "y", "text")) %>%
      layout(dragmode = "zoom", hovermode = "closest") %>%
      config(displayModeBar = input$go_gsea_toolbox,
             modeBarButtonsToRemove = list("toImage"),
             displaylogo = FALSE) %>%
      plotly::toWebGL()
  })
  
  output$go_gsea_status <- renderText({
    if (input$run_go_gsea == 0) return("")
    results <- go_gsea_results()
    if (is.null(results)) return("Aucun résultat disponible")
    paste0("Calcul terminé (", length(results), " ontologie(s))")
  })
  
  output$downloadGoGsea <- downloadHandler(
    filename = function() {
      paste0("GO_GSEA_", input$go_gsea_selected_plot, "_",
             input$go_gsea_displayed_ontology, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      gse <- go_gsea_results()[[input$go_gsea_displayed_ontology]]
      gse <- filter_gsea_by_direction(gse, input$go_gsea_direction)
      
      label <- if (is.null(input$go_gsea_title) || input$go_gsea_title == "") {
        paste0("GO GSEA - ", input$go_gsea_displayed_ontology)
      } else input$go_gsea_title
      top_n <- input$go_gsea_top_n_terms
      
      selected_plot <- switch(input$go_gsea_selected_plot,
                              "Dotplot"   = generate_dotplot(gse, label = label, top_n = top_n),
                              "Cnetplot"  = generate_cnetplot(gse, label = label, top_n = top_n),
                              "Emapplot"  = generate_emapplot(gse, label = label, top_n = top_n),
                              "Upsetplot" = generate_upsetplot(gse, label = label),
                              "Heatplot"  = generate_heatplot(gse, label = label, top_n = top_n),
                              "Ridgeplot" = generate_ridgeplot(gse, label = label, top_n = top_n),
                              "GSEAplot2" = generate_gseaplot2(gse, gene_set_ids = 1:min(3, nrow(as.data.frame(gse)))),
                              "GSEArank"  = generate_gsearank(gse, gene_set_id = 1)
      )
      ggsave(file, plot = selected_plot, width = 12, height = 8, dpi = 300)
    }
  )
  
  
  ####============================ONGLET PATHWAY ORA============================
  
  #####=======================Calcul Pathway ORA===============================
  
  pathway_ora_results <- eventReactive(input$run_pathway_ora, {
    
    req(processed_data())
    req(length(input$pathway_ora_databases) > 0)
    
    # Loader
    w <- make_waiter(
      id = "pathway_ora_plot",
      message = "Calcul GO GSEA en cours... (plusieurs minutes possibles)"
    )
    w$show()
    on.exit(w$hide())
    
    df <- processed_data()
    significant_genes <- get_directional_genes(df, input$pathway_ora_direction)
    
    if (length(significant_genes) == 0) {
      shinyalert(title = "Aucun gène significatif", type = "warning")
      return(NULL)
    }
    
    org_db <- species_to_orgdb[[input$species]]
    if (is.null(org_db)) {
      shinyalert(title = "Espèce non supportée", type = "error")
      return(NULL)
    }
    
    # Conversion SYMBOL -> ENTREZID (requis pour KEGG/Reactome)
    entrez_ids <- convert_symbols_to_entrez(significant_genes, org_db)
    if (is.null(entrez_ids) || length(entrez_ids) == 0) {
      shinyalert(
        title = "Conversion impossible",
        text  = "Aucun de vos gènes n'a pu être converti en ENTREZID.",
        type  = "error"
      )
      return(NULL)
    }
    
    showNotification(
      paste0("Calcul Pathway ORA sur ", length(input$pathway_ora_databases), " base(s)..."),
      duration = NULL, id = "pathway_ora_running", type = "message"
    )
    
    results <- list()
    
    # KEGG
    if ("KEGG" %in% input$pathway_ora_databases) {
      kegg_code <- species_to_kegg[[input$species]]
      results$KEGG <- tryCatch(
        run_ORA_KEGG(
          gene_list = entrez_ids,
          label     = "KEGG ORA",
          organism  = kegg_code,
          p_adj     = "BH",
          p_cutoff  = 0.05
        ),
        error = function(e) {
          showNotification(paste0("Erreur KEGG : ", e$message), type = "warning")
          NULL
        }
      )
    }
    
    # Reactome
    if ("Reactome" %in% input$pathway_ora_databases) {
      reactome_code <- species_to_reactome[[input$species]]
      results$Reactome <- tryCatch(
        run_ORA_pathway(
          gene_list = entrez_ids,
          label     = "Reactome ORA",
          organism  = reactome_code,
          p_adj     = "BH",
          p_cutoff  = 0.05
        ),
        error = function(e) {
          showNotification(paste0("Erreur Reactome : ", e$message), type = "warning")
          NULL
        }
      )
    }
    
    removeNotification("pathway_ora_running")
    
    valid <- sapply(results, function(r) !is.null(r) && nrow(as.data.frame(r)) > 0)
    if (!any(valid)) {
      shinyalert(title = "Aucun pathway enrichi", type = "warning")
      return(NULL)
    }
    
    showNotification("Calcul Pathway ORA terminé", duration = 4, type = "message")
    results[valid]
  })
  
  #####=======================Update sélecteur DB Pathway ORA================
  
  observeEvent(pathway_ora_results(), {
    req(pathway_ora_results())
    available <- names(pathway_ora_results())
    
    updateSelectInput(session, "pathway_ora_displayed_db",
                      choices = available, selected = available[1])
  })
  
  #####=======================Plot Pathway ORA================================
  
  output$pathway_ora_plot <- renderPlotly({
    
    if (is.null(input$deg_file) || is.null(raw_data()) || input$run_pathway_ora == 0) {
      return(
        plotly_empty() %>%
          layout(
            xaxis = list(visible = FALSE),
            yaxis = list(visible = FALSE),
            annotations = list(list(
              text = "Sélectionnez la ou les base(s) puis cliquez sur « Lancer Pathway ORA »",
              showarrow = FALSE,
              font = list(size = 16, color = "#666")
            ))
          )
      )
    }
    
    req(pathway_ora_results())
    req(input$pathway_ora_displayed_db)
    
    res <- pathway_ora_results()[[input$pathway_ora_displayed_db]]
    if (is.null(res)) return(NULL)
    
    label <- if (is.null(input$pathway_ora_title) || input$pathway_ora_title == "") {
      paste0("Pathway ORA - ", input$pathway_ora_displayed_db)
    } else input$pathway_ora_title
    
    top_n <- input$pathway_ora_top_n_terms
    
    selected_plot <- switch(input$pathway_ora_selected_plot,
                            "Dotplot"   = generate_dotplot(res, label = label, top_n = top_n),
                            "Barplot"   = generate_barplot(res, label = label, top_n = top_n),
                            "Cnetplot"  = generate_cnetplot(res, label = label, top_n = top_n),
                            "Emapplot"  = generate_emapplot(res, label = label, top_n = top_n),
                            "Upsetplot" = generate_upsetplot(res, label = label),
                            "Heatplot"  = generate_heatplot(res, label = label, top_n = top_n)
    )
    
    if (is.null(selected_plot)) return(NULL)
    
    ggplotly(selected_plot, tooltip = c("x", "y", "text")) %>%
      layout(dragmode = "zoom", hovermode = "closest") %>%
      config(displayModeBar = input$pathway_ora_toolbox,
             modeBarButtonsToRemove = list("toImage"),
             displaylogo = FALSE) %>%
      plotly::toWebGL()
  })
  
  output$pathway_ora_status <- renderText({
    if (input$run_pathway_ora == 0) return("")
    results <- pathway_ora_results()
    if (is.null(results)) return("Aucun résultat disponible")
    paste0("Calcul terminé (", length(results), " base(s))")
  })
  
  output$downloadPathwayOra <- downloadHandler(
    filename = function() {
      paste0("Pathway_ORA_", input$pathway_ora_selected_plot, "_",
             input$pathway_ora_displayed_db, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      res <- pathway_ora_results()[[input$pathway_ora_displayed_db]]
      label <- if (is.null(input$pathway_ora_title) || input$pathway_ora_title == "") {
        paste0("Pathway ORA - ", input$pathway_ora_displayed_db)
      } else input$pathway_ora_title
      top_n <- input$pathway_ora_top_n_terms
      
      selected_plot <- switch(input$pathway_ora_selected_plot,
                              "Dotplot"   = generate_dotplot(res, label = label, top_n = top_n),
                              "Barplot"   = generate_barplot(res, label = label, top_n = top_n),
                              "Cnetplot"  = generate_cnetplot(res, label = label, top_n = top_n),
                              "Emapplot"  = generate_emapplot(res, label = label, top_n = top_n),
                              "Upsetplot" = generate_upsetplot(res, label = label),
                              "Heatplot"  = generate_heatplot(res, label = label, top_n = top_n)
      )
      ggsave(file, plot = selected_plot, width = 12, height = 8, dpi = 300)
    }
  )
  
  
  ####============================ONGLET PATHWAY GSEA===========================
  
  #####=======================Calcul Pathway GSEA==============================
  
  pathway_gsea_results <- eventReactive(input$run_pathway_gsea, {
    
    req(processed_data())
    req(length(input$pathway_gsea_databases) > 0)
    
    # Loader
    w <- make_waiter(
      id = "pathway_gsea_plot",
      message = "Calcul GO GSEA en cours... (plusieurs minutes possibles)"
    )
    w$show()
    on.exit(w$hide())
    
    ranked <- gsea_ranked()
    if (length(ranked) == 0) {
      shinyalert(title = "Aucun gène valide", type = "warning")
      return(NULL)
    }
    
    org_db <- species_to_orgdb[[input$species]]
    if (is.null(org_db)) {
      shinyalert(title = "Espèce non supportée", type = "error")
      return(NULL)
    }
    
    # Conversion du ranking SYMBOL -> ENTREZID
    ranked_entrez <- convert_ranked_to_entrez(ranked, org_db)
    if (is.null(ranked_entrez) || length(ranked_entrez) == 0) {
      shinyalert(
        title = "Conversion impossible",
        text  = "Aucun de vos gènes n'a pu être converti en ENTREZID.",
        type  = "error"
      )
      return(NULL)
    }
    
    showNotification(
      paste0("Calcul Pathway GSEA sur ", length(input$pathway_gsea_databases),
             " base(s)... Peut prendre plusieurs minutes."),
      duration = NULL, id = "pathway_gsea_running", type = "message"
    )
    
    results <- list()
    
    # KEGG
    if ("KEGG" %in% input$pathway_gsea_databases) {
      kegg_code <- species_to_kegg[[input$species]]
      results$KEGG <- tryCatch(
        run_gsea_kegg(
          ranked_gene_list = ranked_entrez,
          label            = "KEGG GSEA",
          organism         = kegg_code,
          p_adj            = "BH",
          p_cutoff         = 0.05
        ),
        error = function(e) {
          showNotification(paste0("Erreur KEGG : ", e$message), type = "warning")
          NULL
        }
      )
    }
    
    # Reactome
    if ("Reactome" %in% input$pathway_gsea_databases) {
      reactome_code <- species_to_reactome[[input$species]]
      results$Reactome <- tryCatch(
        run_gsea_reactome(
          ranked_gene_list = ranked_entrez,
          label            = "Reactome GSEA",
          organism         = reactome_code,
          p_adj            = "BH",
          p_cutoff         = 0.05
        ),
        error = function(e) {
          showNotification(paste0("Erreur Reactome : ", e$message), type = "warning")
          NULL
        }
      )
    }
    
    removeNotification("pathway_gsea_running")
    
    valid <- sapply(results, function(r) !is.null(r) && nrow(as.data.frame(r)) > 0)
    if (!any(valid)) {
      shinyalert(title = "Aucun pathway enrichi", type = "warning")
      return(NULL)
    }
    
    showNotification("Calcul Pathway GSEA terminé", duration = 4, type = "message")
    results[valid]
  })
  
  #####=======================Update sélecteur DB Pathway GSEA===============
  
  observeEvent(pathway_gsea_results(), {
    req(pathway_gsea_results())
    available <- names(pathway_gsea_results())
    
    updateSelectInput(session, "pathway_gsea_displayed_db",
                      choices = available, selected = available[1])
  })
  
  #####=======================Plot Pathway GSEA===============================
  
  output$pathway_gsea_plot <- renderPlotly({
    
    if (is.null(input$deg_file) || is.null(raw_data()) || input$run_pathway_gsea == 0) {
      return(
        plotly_empty() %>%
          layout(
            xaxis = list(visible = FALSE),
            yaxis = list(visible = FALSE),
            annotations = list(list(
              text = "Sélectionnez la ou les base(s) puis cliquez sur « Lancer Pathway GSEA »",
              showarrow = FALSE,
              font = list(size = 16, color = "#666")
            ))
          )
      )
    }
    
    req(pathway_gsea_results())
    req(input$pathway_gsea_displayed_db)
    
    gse <- pathway_gsea_results()[[input$pathway_gsea_displayed_db]]
    if (is.null(gse)) return(NULL)
    
    # Filtrage par direction NES
    gse <- filter_gsea_by_direction(gse, input$pathway_gsea_direction)
    if (is.null(gse) || nrow(as.data.frame(gse)) == 0) {
      return(
        plotly_empty() %>%
          layout(annotations = list(list(
            text = "Aucun pathway dans cette direction",
            showarrow = FALSE, font = list(size = 16, color = "#666")
          )))
      )
    }
    
    label <- if (is.null(input$pathway_gsea_title) || input$pathway_gsea_title == "") {
      paste0("Pathway GSEA - ", input$pathway_gsea_displayed_db)
    } else input$pathway_gsea_title
    
    top_n <- input$pathway_gsea_top_n_terms
    
    selected_plot <- switch(input$pathway_gsea_selected_plot,
                            "Dotplot"   = generate_dotplot(gse, label = label, top_n = top_n),
                            "Cnetplot"  = generate_cnetplot(gse, label = label, top_n = top_n),
                            "Emapplot"  = generate_emapplot(gse, label = label, top_n = top_n),
                            "Upsetplot" = generate_upsetplot(gse, label = label),
                            "Heatplot"  = generate_heatplot(gse, label = label, top_n = top_n),
                            "Ridgeplot" = generate_ridgeplot(gse, label = label, top_n = top_n),
                            "GSEAplot2" = generate_gseaplot2(gse, gene_set_ids = 1:min(3, nrow(as.data.frame(gse)))),
                            "GSEArank"  = generate_gsearank(gse, gene_set_id = 1)
    )
    
    if (is.null(selected_plot)) return(NULL)
    
    ggplotly(selected_plot, tooltip = c("x", "y", "text")) %>%
      layout(dragmode = "zoom", hovermode = "closest") %>%
      config(displayModeBar = input$pathway_gsea_toolbox,
             modeBarButtonsToRemove = list("toImage"),
             displaylogo = FALSE) %>%
      plotly::toWebGL()
  })
  
  output$pathway_gsea_status <- renderText({
    if (input$run_pathway_gsea == 0) return("")
    results <- pathway_gsea_results()
    if (is.null(results)) return("Aucun résultat disponible")
    paste0("Calcul terminé (", length(results), " base(s))")
  })
  
  output$downloadPathwayGsea <- downloadHandler(
    filename = function() {
      paste0("Pathway_GSEA_", input$pathway_gsea_selected_plot, "_",
             input$pathway_gsea_displayed_db, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      gse <- pathway_gsea_results()[[input$pathway_gsea_displayed_db]]
      gse <- filter_gsea_by_direction(gse, input$pathway_gsea_direction)
      
      label <- if (is.null(input$pathway_gsea_title) || input$pathway_gsea_title == "") {
        paste0("Pathway GSEA - ", input$pathway_gsea_displayed_db)
      } else input$pathway_gsea_title
      top_n <- input$pathway_gsea_top_n_terms
      
      selected_plot <- switch(input$pathway_gsea_selected_plot,
                              "Dotplot"   = generate_dotplot(gse, label = label, top_n = top_n),
                              "Cnetplot"  = generate_cnetplot(gse, label = label, top_n = top_n),
                              "Emapplot"  = generate_emapplot(gse, label = label, top_n = top_n),
                              "Upsetplot" = generate_upsetplot(gse, label = label),
                              "Heatplot"  = generate_heatplot(gse, label = label, top_n = top_n),
                              "Ridgeplot" = generate_ridgeplot(gse, label = label, top_n = top_n),
                              "GSEAplot2" = generate_gseaplot2(gse, gene_set_ids = 1:min(3, nrow(as.data.frame(gse)))),
                              "GSEArank"  = generate_gsearank(gse, gene_set_id = 1)
      )
      ggsave(file, plot = selected_plot, width = 12, height = 8, dpi = 300)
    }
  )
  
}