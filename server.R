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
  
  
  ####============================ONGLET ENRICHISSEMENT=========================
  
  #####=======================Calcul ORA (lancé sur clic bouton)================
  
  # eventReactive : ne se déclenche QUE quand l'utilisateur clique sur le bouton.
  # Le résultat est ensuite mis en cache : changer l'ontologie affichée ou le
  # type de plot ne relance PAS le calcul.
  
  #####=========================MAPPING ESPÈCE -> ORG.DB========================
  # Table de correspondance entre le nom de l'espèce affiché côté UI
  # et l'objet OrgDb à passer à clusterProfiler
  species_to_orgdb <- list(
    "Homo sapiens"            = org.Hs.eg.db::org.Hs.eg.db,
    "Mus musculus"            = org.Mm.eg.db::org.Mm.eg.db,
    "Drosophila melanogaster" = org.Dm.eg.db::org.Dm.eg.db
  )  
  
  # Indicateur d'état du calcul (réactif pour déclencher les updates UI)
  calculation_running <- reactiveVal(FALSE)
  
  # Déclencheur : dès qu'on clique sur le bouton, on signale "en cours"
  # La priorité haute garantit que ça s'exécute avant ora_results
  observeEvent(input$run_enrichment, {
    calculation_running(TRUE)
  }, priority = 10)
  
  ora_results <- eventReactive(input$run_enrichment, {
    
    req(input$enrichment_method == "ORA")
    req(processed_data())
    req(length(input$go_ontology) > 0)
    
    # Active l'état "calcul en cours"
    calculation_running(TRUE)
    on.exit(calculation_running(FALSE))
    
    # FORCE Shiny à pousser l'état au navigateur AVANT de bloquer
    session$sendCustomMessage("dummy", list())
    Sys.sleep(0.05)  # laisse 50ms au navigateur pour rafraîchir
    
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
      paste0("Calcul ORA en cours sur ", length(input$go_ontology), 
             " ontologie(s)..."),
      duration = NULL,
      id       = "ora_running",
      type     = "message"
    )
    
    # Lancement du calcul
    # withProgress force Shiny à mettre à jour l'UI avant le calcul
    ego_list <-
        tryCatch(
          run_ORA_go(
            gene_list = significant_genes,
            label     = "Enrichissement",
            org_db    = org_db,
            ontology  = input$go_ontology,
            p_adj     = "BH",
            p_cutoff  = 0.05,
            key_type  = "SYMBOL"
          ),
          error = function(e) NULL
        )
    
    
    removeNotification("ora_running")
    
    # Vérification que le mapping a réussi pour au moins une ontologie
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
      ego_list <- ego_list[valid_results]
    }
    
    showNotification("Calcul ORA terminé", duration = 4, type = "message")
    
    ego_list[valid_results]
  })
  
  #####=======================Mise à jour du sélecteur d'ontologie==============
  
  # Après chaque calcul ORA, on restreint le selectInput "displayed_ontology"
  # aux ontologies qui ont effectivement été calculées
  observeEvent(ora_results(), {
    available <- names(ora_results())
    
    # Labels en français pour l'affichage
    labels_map <- c(
      "BP" = "Processus Biologique (BP)",
      "CC" = "Composant Cellulaire (CC)",
      "MF" = "Fonction Moléculaire (MF)"
    )
    
    choices <- setNames(available, labels_map[available])
    
    updateSelectInput(
      session,
      inputId  = "displayed_ontology",
      choices  = choices,
      selected = available[1]
    )
  })
  
  #####=======================Génération des plots ORA==========================
  
  # Construit la liste des plots disponibles pour l'ontologie sélectionnée.
  # Cette étape est légère car elle réutilise l'objet ego déjà calculé.
  ora_plots <- reactive({
    req(ora_results())
    req(input$displayed_ontology)
    
    ego <- ora_results()[[input$displayed_ontology]]
    
    # Si l'ontologie sélectionnée n'a pas été calculée, on sort
    if (is.null(ego)) return(NULL)
    
    # Titre : celui saisi par l'utilisateur, sinon un titre par défaut
    label <- if (is.null(input$enrichment_title) || input$enrichment_title == "") {
      paste0("Enrichissement GO-", input$displayed_ontology)
    } else {
      input$enrichment_title
    }
    
    plot_ORA(ego, label = label, top_n = input$top_n_terms)
  })
  
  #####=======================Affichage de l'erreur enrichissement==============
  
  # Affiche l'image d'erreur tant qu'aucun calcul n'a été lancé OU si le résultat est vide
  output$show_enrichment_error <- reactive({
    is.null(input$deg_file) || 
      is.null(raw_data()) || 
      input$run_enrichment == 0
  })
  outputOptions(output, "show_enrichment_error", suspendWhenHidden = FALSE)
  
  output$enrichment_error_img <- renderImage({
    list(
      src         = "www/erreur_format.jpg",  # TODO : image dédiée à l'enrichissement
      contentType = "image/jpeg",
      width       = "70%",
      height      = "auto",
      alt         = "Aucun calcul lancé ou aucun résultat disponible"
    )
  }, deleteFile = FALSE)
  
  #####=======================Rendu du plot d'enrichissement====================
  
  output$enrichment_plot <- renderPlotly({
    
    if (is.null(input$deg_file) || is.null(raw_data()) || input$run_enrichment == 0) {
      return(
        plotly_empty(type = "scatter", mode = "markers") %>%
          layout(
            xaxis = list(visible = FALSE),
            yaxis = list(visible = FALSE),
            annotations = list(
              list(
                text = "Sélectionner votre méthode d'enrichissement, le(s) ontologie et cliquez sur « Lancer l'enrichissement »",
                showarrow = FALSE,
                font = list(size = 16, color = "#666")
              )
            )
          )
      )
    }
    
    req(ora_plots())
    
    plot_choice <- if (input$enrichment_method == "ORA") {
      input$selected_plot_ora
    } else {
      input$selected_plot_gsea
    }
    
    selected_plot <- ora_plots()[[plot_choice]]
    if (is.null(selected_plot)) return(NULL)
    
    ggplotly(selected_plot, tooltip = c("x", "y", "text")) %>%
      layout(dragmode = "zoom", hovermode = "closest") %>%
      config(
        displayModeBar = input$enrichment_toolbox,
        modeBarButtonsToRemove = list("toImage"),
        displaylogo = FALSE
      ) %>%
      plotly::toWebGL()
  })
  
  #####=======================Statut affiché sous le bouton=====================
  
  # Petit texte d'aide qui informe l'utilisateur de l'état du calcul
  output$enrichment_status <- renderText({
    if (input$run_enrichment == 0) {
      return("")  # rien affiché tant qu'aucun calcul lancé
    }
    
    if (calculation_running()) {
      return("")  # rien pendant le calcul (le spinner suffit)
    }
    
    results <- ora_results()
    if (is.null(results)) {
      return("Aucun résultat disponible")
    }
    
    n_ont <- length(results)
    paste0("Calcul terminé (", n_ont, " ontologie(s))")
  })
  
  #####=======================Téléchargement du plot d'enrichissement===========
  
  output$downloadEnrichment <- downloadHandler(
    filename = function() {
      plot_choice <- if (input$enrichment_method == "ORA") {
        input$selected_plot_ora
      } else {
        input$selected_plot_gsea
      }
      paste0("Enrichissement_", plot_choice, "_", 
             input$displayed_ontology, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      plot_choice <- if (input$enrichment_method == "ORA") {
        input$selected_plot_ora
      } else {
        input$selected_plot_gsea
      }
      selected_plot <- ora_plots()[[plot_choice]]
      ggsave(file, plot = selected_plot, width = 12, height = 8, dpi = 300)
    }
  )
  
}