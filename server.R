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
    
    #####=======================Accueil interactif===============================
    
    # Initialiser une valeur réactive pour savoir quelle bulle est sélectionnée
    selected_bubble <- reactiveVal(NULL)
    
    # Observer les clics sur les bulles
    observeEvent(input$btn_deg, {
        selected_bubble("deg")
    })
    
    observeEvent(input$btn_go, {
        selected_bubble("go")
    })
    
    observeEvent(input$btn_pathway, {
        selected_bubble("pathway")
    })
    
    # Générer l'UI de l'explication
    output$explication_accueil <- renderUI({
        
        req(selected_bubble())
        
        if (selected_bubble() == "deg") {
            
            div(
                class = "explication-box",
                
                h3("I) Inspection des Données (DEG)"),
                
                p("Notre application permet une représentation visuelle interactive de vos données d’expression génique via un Volcano Plot dynamique. Cet outil est indispensable pour cibler rapidement les gènes significativement surexprimés ou sous-exprimés."),
                
                p("Ajustez en temps réel vos seuils de P-value et de Log2 Fold Change pour affiner vos résultats, et croisez ces visualisations avec notre tableau de données interactif pour une inspection ciblée de vos gènes d’intérêt.")
            )
            
        } else if (selected_bubble() == "go") {
            
            div(
                class = "explication-box",
                
                h3("II) Enrichissement de Termes GO"),
                
                h4("A) Enrichissement ORA (Over-Representation Analysis)"),
                
                p("Grâce à l’intégration des bases de données de la Gene Ontology, identifiez efficacement les termes GO surreprésentés..."),
                
                h4("B) Enrichissement GSEA (Gene Set Enrichment Analysis)"),
                
                p("Allez plus loin dans l’interprétation en évaluant l’enrichissement sur l’ensemble de votre profil d’expression...")
            )
            
        } else if (selected_bubble() == "pathway") {
            
            div(
                class = "explication-box",
                
                h3("III) Analyses d'enrichissement de voies biologiques"),
                
                h4("A) Enrichissement KEGG"),
                
                p("Allez au-delà des fonctions isolées et cartographiez vos gènes significatifs directement sur les voies métaboliques et de signalisation cellulaires..."),
                
                h4("B) Enrichissement REACTOME"),
                
                p("Plongez dans un réseau de réactions moléculaires extrêmement détaillé et rigoureusement documenté...")
            )
        }
    })
    
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
    
    species_to_orgdb <- list(
        "Homo sapiens"            = org.Hs.eg.db::org.Hs.eg.db,
        "Mus musculus"            = org.Mm.eg.db::org.Mm.eg.db,
        "Drosophila melanogaster" = org.Dm.eg.db::org.Dm.eg.db
    )
    
    species_to_kegg <- c(
        "Homo sapiens"            = "hsa",
        "Mus musculus"            = "mmu",
        "Drosophila melanogaster" = "dme"
    )
    
    species_to_reactome <- c(
        "Homo sapiens"            = "human",
        "Mus musculus"            = "mouse",
        "Drosophila melanogaster" = "fly"
    )
    
    #####=======================Sélection des gènes par direction=================
    get_directional_genes <- function(df, direction) {
        
        df_sig <- df[as.character(df$Significance) != "Not significant", ]
        
        if (direction == "up") {
            
            return(
                df_sig$GeneName[
                    as.character(df_sig$Significance) == "Upregulated"
                ]
            )
            
        } else if (direction == "down") {
            
            return(
                df_sig$GeneName[
                    as.character(df_sig$Significance) == "Downregulated"
                ]
            )
        }
        
        df_sig$GeneName
    }
  
  ####============================ONGLET DEG====================================
  
  #####=======================Volcano plot======================================
  
  create_volcano <- reactive({
    req(processed_data())
    selected_row <- input$deg_table_rows_selected
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
    ggplotly(create_volcano(), tooltip = c("x", "y", "colour", "text")) %>%
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
    # Alerte si fichier non présent
    if (is.null(input$deg_file) || is.null(raw_data())) {
      shinyalert(title = "Fichier manquant", text = "Veuillez d'abord importer un fichier CSV", type = "error")
      return(NULL)
    }
    
    req(processed_data())
    req(length(input$go_ora_ontology) > 0)
    
    w <- make_waiter(id = "ui_go_ora_plot", message = paste0("Calcul GO ORA en cours sur ", length(input$go_ora_ontology), " ontologie(s)..."))
    w$show()
    on.exit(w$hide())
    
    significant_genes <- get_directional_genes(processed_data(), input$go_ora_direction)
    
    if (length(significant_genes) == 0) {
      shinyalert(title = "Aucun gène significatif", text = "Ajustez vos seuils ou la direction.", type = "warning")
      return(NULL)
    }
    
    org_db <- species_to_orgdb[[input$species]]
    
    showNotification("Calcul GO ORA en cours...", duration = NULL, id = "go_ora_running", type = "message")
    
    ego_list <- tryCatch(
      run_ORA_go(gene_list = significant_genes, label = "GO ORA", org_db = org_db, ontology = input$go_ora_ontology, p_adj = "BH", p_cutoff = 0.05, key_type = "SYMBOL"),
      error = function(e) NULL
    )
    
    removeNotification("go_ora_running")
    if (is.null(ego_list)) return(NULL)
    
    valid <- !sapply(ego_list, is.null)
    if (!any(valid)) {
      shinyalert(title = "Aucun gène mappé", text = paste0("Vérifiez l'espèce sélectionnée."), type = "error")
      return(NULL)
    }
    showNotification("Calcul GO ORA terminé", duration = 4, type = "message")
    ego_list[valid]
  })
  
  observeEvent(go_ora_results(), {
    req(go_ora_results())
    available <- names(go_ora_results())
    labels_map <- c("BP" = "Processus Biologique (BP)", "CC" = "Composant Cellulaire (CC)", "MF" = "Fonction Moléculaire (MF)")
    updateSelectInput(session, "go_ora_displayed_ontology", choices = setNames(available, labels_map[available]), selected = available[1])
  })
  
  #####=======================Plot GO ORA======================================
  
  go_ora_plot_obj <- reactive({
    if (is.null(input$deg_file) || is.null(raw_data()) || input$run_go_ora == 0) return(NULL)
    req(go_ora_results(), input$go_ora_displayed_ontology)
    
    ego <- go_ora_results()[[input$go_ora_displayed_ontology]]
    if (is.null(ego)) return(NULL)
    
    label <- if (is.null(input$go_ora_title) || input$go_ora_title == "") paste0("GO ORA - ", input$go_ora_displayed_ontology) else input$go_ora_title
    top_n <- input$go_ora_top_n_terms
    
    switch(input$go_ora_selected_plot,
           "Dotplot"   = generate_dotplot(ego, label = label, top_n = top_n),
           "Barplot"   = generate_barplot(ego, label = label, top_n = top_n),
           "Cnetplot"  = generate_cnetplot(ego, label = label, top_n = top_n),
           "Emapplot"  = generate_emapplot(ego, label = label, top_n = top_n),
           "Goplot"    = generate_goplot(ego, label = label, top_n = top_n),
           "Upsetplot" = generate_upsetplot(ego, label = label),
           "Heatplot"  = generate_heatplot(ego, label = label, top_n = top_n),
           "Manhattan" = generate_manhattan_plot(go_results = go_ora_results(), label = label, p_cutoff = 0.05, top_n_labels = min(top_n, 15))
    )
  })
  
  output$ui_go_ora_plot <- renderUI({
    if (is.null(input$deg_file) || is.null(raw_data()) || input$run_go_ora == 0) {
      return(plotlyOutput("go_ora_empty", height = "500px"))
    }
    # Ces graphiques sont rendus en Plot Statique natif car ggraph est cassé par ggplotly()
    static_plots <- c("Cnetplot", "Emapplot", "Goplot", "Upsetplot")
    if (input$go_ora_selected_plot %in% static_plots) {
      plotOutput("go_ora_plot_static", height = "500px")
    } else {
      plotlyOutput("go_ora_plot_interactive", height = "500px")
    }
  })
  
  output$go_ora_empty <- renderPlotly({
    plotly_empty(type = "scatter", mode = "markers") %>%
      layout(xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
             annotations = list(list(text = "Sélectionnez la ou les ontologie(s) puis cliquez sur « Lancer GO ORA »", showarrow = FALSE, font = list(size = 16, color = "#666"))))
  })
  
  output$go_ora_plot_static <- renderPlot({
    req(go_ora_plot_obj())
    go_ora_plot_obj()
  })
  
  output$go_ora_plot_interactive <- renderPlotly({
    req(go_ora_plot_obj())
    ggplotly(go_ora_plot_obj(), tooltip = c("x", "y", "text")) %>%
      layout(dragmode = "zoom", hovermode = "closest") %>%
      config(displayModeBar = input$go_ora_toolbox, modeBarButtonsToRemove = list("toImage"), displaylogo = FALSE) %>%
      plotly::toWebGL()
  })
  
  output$go_ora_status <- renderText({
    if (input$run_go_ora == 0) return("")
    if (is.null(go_ora_results())) return("Aucun résultat disponible")
    paste0("Calcul terminé (", length(go_ora_results()), " ontologie(s))")
  })
  
  output$downloadGoOra <- downloadHandler(
    filename = function() paste0("GO_ORA_", input$go_ora_selected_plot, "_", input$go_ora_displayed_ontology, "_", Sys.Date(), ".png"),
    content = function(file) {
      req(go_ora_plot_obj())
      ggsave(file, plot = go_ora_plot_obj(), width = 12, height = 8, dpi = 300)
    }
  )
  
  ####============================ONGLET GO GSEA================================
  
  go_gsea_results <- eventReactive(input$run_go_gsea, {
    if (is.null(input$deg_file) || is.null(raw_data())) {
      shinyalert(title = "Fichier manquant", text = "Veuillez d'abord importer un fichier CSV valide dans l'onglet Accueil ou DEG.", type = "error")
      return(NULL)
    }
    
    req(processed_data())
    req(length(input$go_gsea_ontology) > 0)
    
    w <- make_waiter(id = "ui_go_gsea_plot", message = "Calcul GO GSEA en cours... (plusieurs minutes possibles)")
    w$show()
    on.exit(w$hide())
    
    ranked <- gsea_ranked()
    org_db <- species_to_orgdb[[input$species]]
    
    showNotification("Calcul GO GSEA en cours...", duration = NULL, id = "go_gsea_running", type = "message")
    
    gse_list <- tryCatch(
      run_gsea_go(ranked_gene_list = ranked, label = "GO GSEA", org_db = org_db, ontology = input$go_gsea_ontology, p_adj = "BH", p_cutoff = 0.05, key_type = "SYMBOL"),
      error = function(e) NULL
    )
    
    removeNotification("go_gsea_running")
    if (is.null(gse_list)) return(NULL)
    
    valid <- sapply(gse_list, function(g) !is.null(g) && nrow(as.data.frame(g)) > 0)
    if (!any(valid)) {
      shinyalert(title = "Aucun terme enrichi", type = "warning")
      return(NULL)
    }
    showNotification("Calcul GO GSEA terminé", duration = 4, type = "message")
    gse_list[valid]
  })
  
  observeEvent(go_gsea_results(), {
    req(go_gsea_results())
    available <- names(go_gsea_results())
    labels_map <- c("BP" = "Processus Biologique (BP)", "CC" = "Composant Cellulaire (CC)", "MF" = "Fonction Moléculaire (MF)")
    updateSelectInput(session, "go_gsea_displayed_ontology", choices = setNames(available, labels_map[available]), selected = available[1])
  })
  
  filter_gsea_by_direction <- function(gse, direction) {
    if (is.null(gse) || nrow(as.data.frame(gse)) == 0) return(gse)
    if (direction == "both") return(gse)
    df <- gse@result
    keep <- if (direction == "up") df$NES > 0 else df$NES < 0
    if (!any(keep)) return(NULL)
    gse@result <- df[keep, , drop = FALSE]
    gse
  }
  
  go_gsea_plot_obj <- reactive({
    if (is.null(input$deg_file) || is.null(raw_data()) || input$run_go_gsea == 0) return(NULL)
    req(go_gsea_results(), input$go_gsea_displayed_ontology)
    
    gse <- go_gsea_results()[[input$go_gsea_displayed_ontology]]
    if (is.null(gse)) return(NULL)
    
    gse <- filter_gsea_by_direction(gse, input$go_gsea_direction)
    if (is.null(gse) || nrow(as.data.frame(gse)) == 0) return("EMPTY_DIR")
    
    label <- if (is.null(input$go_gsea_title) || input$go_gsea_title == "") paste0("GO GSEA - ", input$go_gsea_displayed_ontology) else input$go_gsea_title
    top_n <- input$go_gsea_top_n_terms
    
    switch(input$go_gsea_selected_plot,
           "Dotplot"   = generate_dotplot(gse, label = label, top_n = top_n),
           "Cnetplot"  = generate_cnetplot(gse, label = label, top_n = top_n),
           "Emapplot"  = generate_emapplot(gse, label = label, top_n = top_n),
           "Upsetplot" = generate_upsetplot(gse, label = label),
           "Heatplot"  = generate_heatplot(gse, label = label, top_n = top_n),
           "Ridgeplot" = generate_ridgeplot(gse, label = label, top_n = top_n),
           "GSEAplot2" = generate_gseaplot2(gse, gene_set_ids = 1:min(3, nrow(as.data.frame(gse)))),
           "GSEArank"  = generate_gsearank(gse, gene_set_id = 1),
           "Manhattan" = generate_manhattan_plot(go_results = go_gsea_results(), label = label, p_cutoff = 0.05, top_n_labels = min(top_n, 15))
    )
  })
  
  output$ui_go_gsea_plot <- renderUI({
    if (is.null(input$deg_file) || is.null(raw_data()) || input$run_go_gsea == 0) {
      return(plotlyOutput("go_gsea_empty", height = "500px"))
    }
    static_plots <- c("Cnetplot", "Emapplot", "Upsetplot", "GSEAplot2", "GSEArank")
    if (input$go_gsea_selected_plot %in% static_plots) {
      plotOutput("go_gsea_plot_static", height = "500px")
    } else {
      plotlyOutput("go_gsea_plot_interactive", height = "500px")
    }
  })
  
  output$go_gsea_empty <- renderPlotly({
    plotly_empty(type = "scatter", mode = "markers") %>%
      layout(xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
             annotations = list(list(text = "Sélectionnez la ou les ontologie(s) puis cliquez sur « Lancer GO GSEA »", showarrow = FALSE, font = list(size = 16, color = "#666"))))
  })
  
  output$go_gsea_plot_static <- renderPlot({
    p <- go_gsea_plot_obj()
    req(p)
    if (is.character(p) && p == "EMPTY_DIR") {
      return(ggplot() + annotate("text", x=0, y=0, label="Aucun pathway dans cette direction") + theme_void())
    }
    p
  })
  
  output$go_gsea_plot_interactive <- renderPlotly({
    p <- go_gsea_plot_obj()
    req(p)
    if (is.character(p) && p == "EMPTY_DIR") {
      return(plotly_empty() %>% layout(annotations = list(list(text = "Aucun pathway dans cette direction", showarrow=FALSE, font=list(size=16, color="#666")))))
    }
    ggplotly(p, tooltip = c("x", "y", "text")) %>%
      layout(dragmode = "zoom", hovermode = "closest") %>%
      config(displayModeBar = input$go_gsea_toolbox, modeBarButtonsToRemove = list("toImage"), displaylogo = FALSE) %>%
      plotly::toWebGL()
  })
  
  output$go_gsea_status <- renderText({
    if (input$run_go_gsea == 0) return("")
    if (is.null(go_gsea_results())) return("Aucun résultat disponible")
    paste0("Calcul terminé (", length(go_gsea_results()), " ontologie(s))")
  })
  
  output$downloadGoGsea <- downloadHandler(
    filename = function() paste0("GO_GSEA_", input$go_gsea_selected_plot, "_", input$go_gsea_displayed_ontology, "_", Sys.Date(), ".png"),
    content = function(file) {
      p <- go_gsea_plot_obj()
      req(p); if (is.character(p)) return(NULL)
      ggsave(file, plot = p, width = 12, height = 8, dpi = 300)
    }
  )
  
  ####============================ONGLET PATHWAY ORA============================
  
  pathway_ora_results <- eventReactive(input$run_pathway_ora, {
    if (is.null(input$deg_file) || is.null(raw_data())) {
      shinyalert(title = "Fichier manquant", text = "Veuillez d'abord importer un fichier CSV valide dans l'onglet Accueil ou DEG.", type = "error")
      return(NULL)
    }
    req(processed_data())
    req(length(input$pathway_ora_databases) > 0)
    
    w <- make_waiter(id = "ui_pathway_ora_plot", message = "Calcul Pathway ORA en cours...")
    w$show()
    on.exit(w$hide())
    
    significant_genes <- get_directional_genes(processed_data(), input$pathway_ora_direction)
    if (length(significant_genes) == 0) {
      shinyalert(title = "Aucun gène significatif", type = "warning")
      return(NULL)
    }
    
    org_db <- species_to_orgdb[[input$species]]
    entrez_ids <- convert_symbols_to_entrez(significant_genes, org_db)
    if (is.null(entrez_ids) || length(entrez_ids) == 0) {
      shinyalert(title = "Conversion impossible", text = "Aucun de vos gènes n'a pu être converti en ENTREZID.", type = "error")
      return(NULL)
    }
    
    showNotification("Calcul Pathway ORA en cours...", duration = NULL, id = "pathway_ora_running", type = "message")
    
    results <- list()
    if ("KEGG" %in% input$pathway_ora_databases) {
      results$KEGG <- tryCatch(run_ORA_KEGG(gene_list = entrez_ids, label = "KEGG ORA", organism = species_to_kegg[[input$species]]), error = function(e) NULL)
    }
    if ("Reactome" %in% input$pathway_ora_databases) {
      results$Reactome <- tryCatch(run_ORA_pathway(gene_list = entrez_ids, label = "Reactome ORA", organism = species_to_reactome[[input$species]]), error = function(e) NULL)
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
  
  observeEvent(pathway_ora_results(), {
    req(pathway_ora_results())
    available <- names(pathway_ora_results())
    updateSelectInput(session, "pathway_ora_displayed_db", choices = available, selected = available[1])
  })
  
  pathway_ora_plot_obj <- reactive({
    if (is.null(input$deg_file) || is.null(raw_data()) || input$run_pathway_ora == 0) return(NULL)
    req(pathway_ora_results(), input$pathway_ora_displayed_db)
    
    res <- pathway_ora_results()[[input$pathway_ora_displayed_db]]
    if (is.null(res)) return(NULL)
    
    label <- if (is.null(input$pathway_ora_title) || input$pathway_ora_title == "") paste0("Pathway ORA - ", input$pathway_ora_displayed_db) else input$pathway_ora_title
    top_n <- input$pathway_ora_top_n_terms
    
    switch(input$pathway_ora_selected_plot,
           "Dotplot"   = generate_dotplot(res, label = label, top_n = top_n),
           "Barplot"   = generate_barplot(res, label = label, top_n = top_n),
           "Cnetplot"  = generate_cnetplot(res, label = label, top_n = top_n),
           "Emapplot"  = generate_emapplot(res, label = label, top_n = top_n),
           "Upsetplot" = generate_upsetplot(res, label = label),
           "Heatplot"  = generate_heatplot(res, label = label, top_n = top_n)
    )
  })
  
  output$ui_pathway_ora_plot <- renderUI({
    if (is.null(input$deg_file) || is.null(raw_data()) || input$run_pathway_ora == 0) {
      return(plotlyOutput("pathway_ora_empty", height = "500px"))
    }
    static_plots <- c("Cnetplot", "Emapplot", "Upsetplot")
    if (input$pathway_ora_selected_plot %in% static_plots) {
      plotOutput("pathway_ora_plot_static", height = "500px")
    } else {
      plotlyOutput("pathway_ora_plot_interactive", height = "500px")
    }
  })
  
  output$pathway_ora_empty <- renderPlotly({
    plotly_empty() %>% layout(xaxis = list(visible = FALSE), yaxis = list(visible = FALSE), annotations = list(list(text = "Sélectionnez la ou les base(s) puis cliquez sur « Lancer Pathway ORA »", showarrow = FALSE, font = list(size = 16, color = "#666"))))
  })
  
  output$pathway_ora_plot_static <- renderPlot({
    req(pathway_ora_plot_obj())
    pathway_ora_plot_obj()
  })
  
  output$pathway_ora_plot_interactive <- renderPlotly({
    req(pathway_ora_plot_obj())
    ggplotly(pathway_ora_plot_obj(), tooltip = c("x", "y", "text")) %>%
      layout(dragmode = "zoom", hovermode = "closest") %>%
      config(displayModeBar = input$pathway_ora_toolbox, modeBarButtonsToRemove = list("toImage"), displaylogo = FALSE) %>%
      plotly::toWebGL()
  })
  
  output$pathway_ora_status <- renderText({
    if (input$run_pathway_ora == 0) return("")
    if (is.null(pathway_ora_results())) return("Aucun résultat disponible")
    paste0("Calcul terminé (", length(pathway_ora_results()), " base(s))")
  })
  
  output$downloadPathwayOra <- downloadHandler(
    filename = function() paste0("Pathway_ORA_", input$pathway_ora_selected_plot, "_", input$pathway_ora_displayed_db, "_", Sys.Date(), ".png"),
    content = function(file) {
      req(pathway_ora_plot_obj())
      ggsave(file, plot = pathway_ora_plot_obj(), width = 12, height = 8, dpi = 300)
    }
  )
  
  ####============================ONGLET PATHWAY GSEA===========================
  
  pathway_gsea_results <- eventReactive(input$run_pathway_gsea, {
    if (is.null(input$deg_file) || is.null(raw_data())) {
      shinyalert(title = "Fichier manquant", text = "Veuillez d'abord importer un fichier CSV valide dans l'onglet Accueil ou DEG.", type = "error")
      return(NULL)
    }
    req(processed_data())
    req(length(input$pathway_gsea_databases) > 0)
    
    w <- make_waiter(id = "ui_pathway_gsea_plot", message = "Calcul Pathway GSEA en cours...")
    w$show()
    on.exit(w$hide())
    
    ranked <- gsea_ranked()
    if (length(ranked) == 0) return(NULL)
    
    org_db <- species_to_orgdb[[input$species]]
    ranked_entrez <- convert_ranked_to_entrez(ranked, org_db)
    
    showNotification("Calcul Pathway GSEA en cours...", duration = NULL, id = "pathway_gsea_running", type = "message")
    
    results <- list()
    if ("KEGG" %in% input$pathway_gsea_databases) {
      results$KEGG <- tryCatch(run_gsea_kegg(ranked_gene_list = ranked_entrez, label = "KEGG GSEA", organism = species_to_kegg[[input$species]]), error = function(e) NULL)
    }
    if ("Reactome" %in% input$pathway_gsea_databases) {
      results$Reactome <- tryCatch(run_gsea_reactome(ranked_gene_list = ranked_entrez, label = "Reactome GSEA", organism = species_to_reactome[[input$species]]), error = function(e) NULL)
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
  
  observeEvent(pathway_gsea_results(), {
    req(pathway_gsea_results())
    available <- names(pathway_gsea_results())
    updateSelectInput(session, "pathway_gsea_displayed_db", choices = available, selected = available[1])
  })
  
  pathway_gsea_plot_obj <- reactive({
    if (is.null(input$deg_file) || is.null(raw_data()) || input$run_pathway_gsea == 0) return(NULL)
    req(pathway_gsea_results(), input$pathway_gsea_displayed_db)
    
    gse <- pathway_gsea_results()[[input$pathway_gsea_displayed_db]]
    if (is.null(gse)) return(NULL)
    
    gse <- filter_gsea_by_direction(gse, input$pathway_gsea_direction)
    if (is.null(gse) || nrow(as.data.frame(gse)) == 0) return("EMPTY_DIR")
    
    label <- if (is.null(input$pathway_gsea_title) || input$pathway_gsea_title == "") paste0("Pathway GSEA - ", input$pathway_gsea_displayed_db) else input$pathway_gsea_title
    top_n <- input$pathway_gsea_top_n_terms
    
    switch(input$pathway_gsea_selected_plot,
           "Dotplot"   = generate_dotplot(gse, label = label, top_n = top_n),
           "Cnetplot"  = generate_cnetplot(gse, label = label, top_n = top_n),
           "Emapplot"  = generate_emapplot(gse, label = label, top_n = top_n),
           "Upsetplot" = generate_upsetplot(gse, label = label),
           "Heatplot"  = generate_heatplot(gse, label = label, top_n = top_n),
           "Ridgeplot" = generate_ridgeplot(gse, label = label, top_n = top_n),
           "GSEAplot2" = generate_gseaplot2(gse, gene_set_ids = 1:min(3, nrow(as.data.frame(gse)))),
           "GSEArank"  = generate_gsearank(gse, gene_set_id = 1)
    )
  })
  
  output$ui_pathway_gsea_plot <- renderUI({
    if (is.null(input$deg_file) || is.null(raw_data()) || input$run_pathway_gsea == 0) {
      return(plotlyOutput("pathway_gsea_empty", height = "500px"))
    }
    static_plots <- c("Cnetplot", "Emapplot", "Upsetplot", "GSEAplot2", "GSEArank")
    if (input$pathway_gsea_selected_plot %in% static_plots) {
      plotOutput("pathway_gsea_plot_static", height = "500px")
    } else {
      plotlyOutput("pathway_gsea_plot_interactive", height = "500px")
    }
  })
  
  output$pathway_gsea_empty <- renderPlotly({
    plotly_empty() %>% layout(xaxis = list(visible = FALSE), yaxis = list(visible = FALSE), annotations = list(list(text = "Sélectionnez la ou les base(s) puis cliquez sur « Lancer Pathway GSEA »", showarrow = FALSE, font = list(size = 16, color = "#666"))))
  })
  
  output$pathway_gsea_plot_static <- renderPlot({
    p <- pathway_gsea_plot_obj()
    req(p)
    if (is.character(p) && p == "EMPTY_DIR") return(ggplot() + annotate("text", x=0, y=0, label="Aucun pathway dans cette direction") + theme_void())
    p
  })
  
  output$pathway_gsea_plot_interactive <- renderPlotly({
    p <- pathway_gsea_plot_obj()
    req(p)
    if (is.character(p) && p == "EMPTY_DIR") {
      return(plotly_empty() %>% layout(annotations = list(list(text = "Aucun pathway dans cette direction", showarrow=FALSE, font=list(size=16, color="#666")))))
    }
    ggplotly(p, tooltip = c("x", "y", "text")) %>%
      layout(dragmode = "zoom", hovermode = "closest") %>%
      config(displayModeBar = input$pathway_gsea_toolbox, modeBarButtonsToRemove = list("toImage"), displaylogo = FALSE) %>%
      plotly::toWebGL()
  })
  
  output$pathway_gsea_status <- renderText({
    if (input$run_pathway_gsea == 0) return("")
    if (is.null(pathway_gsea_results())) return("Aucun résultat disponible")
    paste0("Calcul terminé (", length(pathway_gsea_results()), " base(s))")
  })
  
  output$downloadPathwayGsea <- downloadHandler(
    filename = function() paste0("Pathway_GSEA_", input$pathway_gsea_selected_plot, "_", input$pathway_gsea_displayed_db, "_", Sys.Date(), ".png"),
    content = function(file) {
      p <- pathway_gsea_plot_obj()
      req(p); if (is.character(p)) return(NULL)
      ggsave(file, plot = p, width = 12, height = 8, dpi = 300)
    }
  )
}