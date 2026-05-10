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
  
  ####============================CHARGEMENT DES DONNÉES========================
  
  # Colonnes minimales requises dans le CSV uploadé
  required_columns <- c("GeneName", "log2FC", "pval")
  
  #####=======================Lecture du fichier CSV============================
  raw_data <- reactive({
    req(input$deg_file)
    
    # Vérification de l'extension du fichier
    ext <- tools::file_ext(input$deg_file$name)
    if (tolower(ext) != "csv") {
      shinyalert(
        title = "Format non valide",
        text  = "Veuillez importer un fichier .csv uniquement.",
        type  = "error"
      )
      return(NULL)
    }
    
    # Lecture sécurisée du CSV (séparateur ;)
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
    
    # Vérification des colonnes obligatoires
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
  # Évite que le calcul se relance à chaque pixel quand l'utilisateur fait 
  # glisser un slider (latence de 300 ms)
  fc_threshold_debounced     <- debounce(reactive(input$fc_threshold), 300)
  pvalue_threshold_debounced <- debounce(reactive(input$pvalue_threshold), 300)
  
  #####=======================Classification des gènes==========================
  # Ajoute une colonne "Significance" : Upregulated / Downregulated / Not significant
  processed_data <- reactive({
    if (is.null(raw_data())) return(NULL)
    
    significativity(
      data          = raw_data(),
      log2FC_cutoff = fc_threshold_debounced(),
      P_cutoff      = pvalue_threshold_debounced()
    )
  })
  
  #####=========================MAPPING ESPÈCE -> ORG.DB========================
  # Table de correspondance entre le nom de l'espèce affiché côté UI
  # et l'objet OrgDb à passer à clusterProfiler
  species_to_orgdb <- list(
    "Homo sapiens"            = org.Hs.eg.db::org.Hs.eg.db,
    "Mus musculus"            = org.Mm.eg.db::org.Mm.eg.db,
    "Drosophila melanogaster" = org.Dm.eg.db::org.Dm.eg.db
  )
  
  
  ####============================ONGLET DEG====================================
  
  #####=======================Volcano plot======================================
  
  # Génère l'objet ggplot du volcano (utilisé par le rendu et le download)
  create_volcano <- reactive({
    req(processed_data())
    
    selected_row <- input$deg_table_rows_selected  # ligne cliquée dans le tableau
    
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
  
  # Affiche l'image d'erreur si pas de fichier chargé
  output$show_volcano_error <- reactive({
    is.null(input$deg_file) || is.null(raw_data())
  })
  outputOptions(output, "show_volcano_error", suspendWhenHidden = FALSE)
  
  # Image affichée tant qu'aucun fichier n'est chargé
  output$volcano_error_img <- renderImage({
    list(
      src         = "www/erreur_format.jpg",
      contentType = "image/jpeg",
      width       = "70%",
      height      = "auto",
      alt         = "Format de fichier attendu"
    )
  }, deleteFile = FALSE)
  
  # Rendu interactif du volcano plot (ggplot -> plotly)
  output$volcano_plot <- renderPlotly({
    req(input$deg_file)
    
    ggplotly(create_volcano(), tooltip = c("x", "y", "colour")) %>%
      layout(
        dragmode   = "zoom",
        hovermode  = "closest"
      ) %>%
      config(
        displayModeBar         = input$deg_toolbox,
        modeBarButtonsToAdd    = list("drawrect", "eraseshape"),
        modeBarButtonsToRemove = list("toImage"),
        displaylogo            = FALSE
      ) %>%
      plotly::toWebGL()  # WebGL pour gérer les milliers de points
  })
  
  # Téléchargement du volcano plot en PNG
  output$downloadVolcano <- downloadHandler(
    filename = function() {
      paste0("Volcano_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      ggsave(file, plot = create_volcano(), width = 12, height = 8, dpi = 300)
    }
  )
  
  #####=======================Tableau des données===============================
  
  # Affiche le message d'erreur si pas de fichier chargé
  output$show_table_error <- reactive({
    is.null(input$deg_file) || is.null(raw_data())
  })
  outputOptions(output, "show_table_error", suspendWhenHidden = FALSE)
  
  output$table_error_text <- renderText({
    "Veuillez charger un fichier CSV au format attendu pour explorer les données"
  })
  
  # Tableau interactif DataTables
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
  
  
  ####============================ONGLET ENRICHISSEMENT ORA=====================
  
  #####=======================Calcul ORA========================================
  
  # eventReactive : ne se déclenche QUE quand l'utilisateur clique sur le bouton.
  # Le résultat est ensuite mis en cache : changer l'ontologie affichée ou le
  # type de plot ne relance PAS le calcul.
  ora_results <- eventReactive(input$run_ora, {
    
    req(processed_data())
    req(length(input$ora_go_ontology) > 0)
    
    df <- processed_data()
    significant_genes <- df$GeneName[as.character(df$Significance) != "Not significant"]
    
    if (length(significant_genes) == 0) {
      shinyalert(
        title = "Aucun gène significatif",
        text  = "Aucun gène ne passe les seuils actuels. Ajustez les seuils logFC/p-value et relancez.",
        type  = "warning"
      )
      return(NULL)
    }
    
    # Récupération de l'OrgDb correspondant à l'espèce sélectionnée
    org_db <- species_to_orgdb[[input$species]]
    
    if (is.null(org_db)) {
      shinyalert(
        title = "Espèce non supportée",
        text  = paste0("Aucune base d'annotation disponible pour : ", input$species),
        type  = "error"
      )
      return(NULL)
    }
    
    showNotification(
      paste0("Calcul ORA en cours sur ", length(input$ora_go_ontology), 
             " ontologie(s)..."),
      duration = NULL,
      id       = "ora_running",
      type     = "message"
    )
    
    # Lancement du calcul
    ego_list <- tryCatch(
      run_ORA_go(
        gene_list = significant_genes,
        label     = "Enrichissement",
        org_db    = org_db,
        ontology  = input$ora_go_ontology,
        p_adj     = "BH",
        p_cutoff  = 0.05,
        key_type  = "SYMBOL"
      ),
      error = function(e) NULL
    )
    
    removeNotification("ora_running")
    
    # Vérification que le calcul a réussi
    if (is.null(ego_list)) {
      shinyalert(
        title = "Erreur de calcul",
        text  = "Le calcul ORA a échoué. Vérifiez l'espèce et le format des gènes.",
        type  = "error"
      )
      return(NULL)
    }
    
    valid_results <- !sapply(ego_list, is.null)
    
    if (!any(valid_results)) {
      shinyalert(
        title = "Aucun gène mappé",
        text  = paste0(
          "Aucun de vos gènes n'a pu être mappé à la base d'annotation ",
          "pour l'espèce '", input$species, "'.\n\n",
          "Vérifiez que l'espèce sélectionnée correspond bien à vos données ",
          "et que vos identifiants sont des symboles (ex: TP53, BRCA1)."
        ),
        type = "error"
      )
      return(NULL)
    }
    
    # Si certaines ontologies ont échoué, on les retire de la liste
    if (!all(valid_results)) {
      failed <- names(ego_list)[!valid_results]
      showNotification(
        paste0("Avertissement : aucun résultat pour ", paste(failed, collapse = ", ")),
        duration = 6,
        type     = "warning"
      )
    }
    
    showNotification("Calcul ORA terminé", duration = 4, type = "message")
    
    ego_list[valid_results]
  })
  
  #####=======================Mise à jour du sélecteur d'ontologie ORA==========
  
  # Après chaque calcul ORA, on restreint le selectInput aux ontologies calculées
  observeEvent(ora_results(), {
    req(ora_results())
    available <- names(ora_results())
    
    labels_map <- c(
      "BP" = "Processus Biologique (BP)",
      "CC" = "Composant Cellulaire (CC)",
      "MF" = "Fonction Moléculaire (MF)"
    )
    
    choices <- setNames(available, labels_map[available])
    
    updateSelectInput(
      session,
      inputId  = "ora_displayed_ontology",
      choices  = choices,
      selected = available[1]
    )
  })
  
  #####=======================Rendu du plot ORA================================
  
  output$ora_plot <- renderPlotly({
    
    # Cas 1 : pas encore lancé -> plot vide avec message
    if (is.null(input$deg_file) || is.null(raw_data()) || input$run_ora == 0) {
      return(
        plotly_empty(type = "scatter", mode = "markers") %>%
          layout(
            xaxis = list(visible = FALSE),
            yaxis = list(visible = FALSE),
            annotations = list(
              list(
                text = "Sélectionnez la ou les ontologie(s) puis cliquez sur « Lancer l'enrichissement ORA »",
                showarrow = FALSE,
                font = list(size = 16, color = "#666")
              )
            )
          )
      )
    }
    
    # Cas 2 : calcul lancé -> on récupère l'ego de l'ontologie affichée
    req(ora_results())
    req(input$ora_displayed_ontology)
    
    ego <- ora_results()[[input$ora_displayed_ontology]]
    if (is.null(ego)) return(NULL)
    
    # Titre : celui saisi par l'utilisateur, sinon un titre par défaut
    label <- if (is.null(input$ora_title) || input$ora_title == "") {
      paste0("Enrichissement GO-", input$ora_displayed_ontology)
    } else {
      input$ora_title
    }
    
    top_n <- input$ora_top_n_terms
    
    # Génère uniquement le plot demandé (plus rapide qu'une liste complète)
    selected_plot <- switch(input$ora_selected_plot,
                            "Dotplot"   = generate_dotplot(ego, label = label, top_n = top_n),
                            "Barplot"   = generate_barplot(ego, label = label, top_n = top_n),
                            "Cnetplot"  = generate_cnetplot(ego, label = label, top_n = top_n),
                            "Emapplot"  = generate_emapplot(ego, label = label, top_n = top_n),
                            "Goplot"    = generate_goplot(ego, label = label, top_n = top_n),
                            "Upsetplot" = generate_upsetplot(ego, label = label),
                            "Heatplot"  = generate_heatplot(ego, label = label, top_n = top_n)
    )
    
    if (is.null(selected_plot)) return(NULL)
    
    ggplotly(selected_plot, tooltip = c("x", "y", "text")) %>%
      layout(dragmode = "zoom", hovermode = "closest") %>%
      config(
        displayModeBar = input$ora_toolbox,
        modeBarButtonsToRemove = list("toImage"),
        displaylogo = FALSE
      ) %>%
      plotly::toWebGL()
  })
  
  #####=======================Statut affiché sous le bouton ORA================
  
  output$ora_status <- renderText({
    if (input$run_ora == 0) return("")
    
    results <- ora_results()
    if (is.null(results)) return("Aucun résultat disponible")
    
    n_ont <- length(results)
    paste0("Calcul terminé (", n_ont, " ontologie(s))")
  })
  
  #####=======================Téléchargement du plot ORA=======================
  
  output$downloadOra <- downloadHandler(
    filename = function() {
      paste0("Enrichissement_ORA_", input$ora_selected_plot, "_", 
             input$ora_displayed_ontology, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      ego <- ora_results()[[input$ora_displayed_ontology]]
      label <- if (is.null(input$ora_title) || input$ora_title == "") {
        paste0("Enrichissement GO-", input$ora_displayed_ontology)
      } else {
        input$ora_title
      }
      top_n <- input$ora_top_n_terms
      
      selected_plot <- switch(input$ora_selected_plot,
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
  
  
  ####============================ONGLET ENRICHISSEMENT GSEA====================
  
  #####=======================Préparation du ranking GSEA======================
  
  # GSEA prend un vecteur numérique nommé, trié décroissant 
  # (noms = gènes, valeurs = log2FC)
  gsea_ranked <- reactive({
    req(processed_data())
    
    df <- processed_data()
    
    # Retirer NA et dédupliquer (garder le |log2FC| max par gène)
    df <- df[!is.na(df$log2FC) & !is.na(df$GeneName), ]
    df <- df[order(-abs(df$log2FC)), ]
    df <- df[!duplicated(df$GeneName), ]
    
    ranked <- df$log2FC
    names(ranked) <- df$GeneName
    sort(ranked, decreasing = TRUE)
  })
  
  #####=======================Calcul GSEA======================================
  
  gsea_results <- eventReactive(input$run_gsea, {
    
    req(processed_data())
    req(length(input$gsea_go_ontology) > 0)
    
    ranked <- gsea_ranked()
    
    if (length(ranked) == 0) {
      shinyalert(
        title = "Aucun gène valide",
        text  = "Aucun gène valide pour construire le ranking GSEA.",
        type  = "warning"
      )
      return(NULL)
    }
    
    # Récupération de l'OrgDb correspondant à l'espèce sélectionnée
    org_db <- species_to_orgdb[[input$species]]
    
    if (is.null(org_db)) {
      shinyalert(
        title = "Espèce non supportée",
        text  = paste0("Aucune base d'annotation disponible pour : ", input$species),
        type  = "error"
      )
      return(NULL)
    }
    
    showNotification(
      paste0("Calcul GSEA en cours sur ", length(input$gsea_go_ontology), 
             " ontologie(s)... Peut prendre plusieurs minutes."),
      duration = NULL,
      id       = "gsea_running",
      type     = "message"
    )
    
    # Lancement du calcul GSEA
    gse_list <- tryCatch(
      run_gsea_go(
        ranked_gene_list = ranked,
        label            = "GSEA",
        org_db           = org_db,
        ontology         = input$gsea_go_ontology,
        p_adj            = "BH",
        p_cutoff         = 0.05,
        key_type         = "SYMBOL"
      ),
      error = function(e) NULL
    )
    
    removeNotification("gsea_running")
    
    if (is.null(gse_list)) {
      shinyalert(
        title = "Erreur de calcul",
        text  = "Le calcul GSEA a échoué. Vérifiez l'espèce et le format des gènes.",
        type  = "error"
      )
      return(NULL)
    }
    
    # Filtrage des ontologies vides (aucun terme enrichi)
    valid_results <- sapply(gse_list, function(g) !is.null(g) && nrow(as.data.frame(g)) > 0)
    
    if (!any(valid_results)) {
      shinyalert(
        title = "Aucun terme enrichi",
        text  = "Aucun terme GO n'est significativement enrichi avec ces paramètres.",
        type  = "warning"
      )
      return(NULL)
    }
    
    if (!all(valid_results)) {
      failed <- names(gse_list)[!valid_results]
      showNotification(
        paste0("Avertissement : aucun résultat pour ", paste(failed, collapse = ", ")),
        duration = 6,
        type     = "warning"
      )
    }
    
    showNotification("Calcul GSEA terminé", duration = 4, type = "message")
    
    gse_list[valid_results]
  })
  
  #####=======================Mise à jour du sélecteur d'ontologie GSEA========
  
  observeEvent(gsea_results(), {
    req(gsea_results())
    available <- names(gsea_results())
    
    labels_map <- c(
      "BP" = "Processus Biologique (BP)",
      "CC" = "Composant Cellulaire (CC)",
      "MF" = "Fonction Moléculaire (MF)"
    )
    
    choices <- setNames(available, labels_map[available])
    
    updateSelectInput(
      session,
      inputId  = "gsea_displayed_ontology",
      choices  = choices,
      selected = available[1]
    )
  })
  
  #####=======================Rendu du plot GSEA===============================
  
  output$gsea_plot <- renderPlotly({
    
    # Cas 1 : pas encore lancé -> plot vide avec message
    if (is.null(input$deg_file) || is.null(raw_data()) || input$run_gsea == 0) {
      return(
        plotly_empty(type = "scatter", mode = "markers") %>%
          layout(
            xaxis = list(visible = FALSE),
            yaxis = list(visible = FALSE),
            annotations = list(
              list(
                text = "Sélectionnez la ou les ontologie(s) puis cliquez sur « Lancer l'enrichissement GSEA »",
                showarrow = FALSE,
                font = list(size = 16, color = "#666")
              )
            )
          )
      )
    }
    
    # Cas 2 : calcul lancé -> on récupère le gse de l'ontologie affichée
    req(gsea_results())
    req(input$gsea_displayed_ontology)
    
    gse <- gsea_results()[[input$gsea_displayed_ontology]]
    if (is.null(gse)) return(NULL)
    
    label <- if (is.null(input$gsea_title) || input$gsea_title == "") {
      paste0("GSEA GO-", input$gsea_displayed_ontology)
    } else {
      input$gsea_title
    }
    
    top_n <- input$gsea_top_n_terms
    
    # Génère uniquement le plot demandé
    selected_plot <- switch(input$gsea_selected_plot,
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
      config(
        displayModeBar = input$gsea_toolbox,
        modeBarButtonsToRemove = list("toImage"),
        displaylogo = FALSE
      ) %>%
      plotly::toWebGL()
  })
  
  #####=======================Statut affiché sous le bouton GSEA===============
  
  output$gsea_status <- renderText({
    if (input$run_gsea == 0) return("")
    
    results <- gsea_results()
    if (is.null(results)) return("Aucun résultat disponible")
    
    n_ont <- length(results)
    paste0("Calcul terminé (", n_ont, " ontologie(s))")
  })
  
  #####=======================Téléchargement du plot GSEA======================
  
  output$downloadGsea <- downloadHandler(
    filename = function() {
      paste0("Enrichissement_GSEA_", input$gsea_selected_plot, "_", 
             input$gsea_displayed_ontology, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      gse <- gsea_results()[[input$gsea_displayed_ontology]]
      label <- if (is.null(input$gsea_title) || input$gsea_title == "") {
        paste0("GSEA GO-", input$gsea_displayed_ontology)
      } else {
        input$gsea_title
      }
      top_n <- input$gsea_top_n_terms
      
      selected_plot <- switch(input$gsea_selected_plot,
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