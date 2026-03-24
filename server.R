#==================================================================================================
#Author : Noé Méderlet 
#contact : noe.mederlet@univ-rouen.fr
#github : https://github.com/Medernoe/Projet_omics
#organism : Master Bims M2, université de rouen 
#project : Création d'une application interactive dédiée à l'analyse de données transcriptomiques,
#développée dans le cadre d'un projet universitaire du Master 2 de Bioinformatique de l'Université de Rouen.
#==================================================================================================
source("global.R")

# Serveur
function(input, output, session) {
  
  #==================================================================================================
  # CHARGEMENT DES DONNÉES
  
  # Lecture du fichier CSV uploadé par l'utilisateur
  required_columns <- c("GeneName", "log2FC", "pval")
  
  data <- reactive({
    req(input$file)
    
    # Vérification de l’extension du fichier
    ext <- tools::file_ext(input$file$name)
    if (tolower(ext) != "csv") {
      shinyalert(
        title = "Format non valide",
        text = "Veuillez importer un fichier .csv uniquement.",
        type = "error"
      )
      return(NULL)
    }
    
    # Lecture sécurisée du CSV
    df <- tryCatch(
      {
        read.csv(input$file$datapath, sep = ";")
      },
      error = function(e) {
        shinyalert(
          title = "Erreur de lecture",
          text = "Impossible de lire le fichier. Vérifiez qu'il s'agit d’un CSV valide.",
          type = "error"
        )
        return(NULL)
      }
    )
    
    if (is.null(df)) return(NULL)
    
    # Vérification des colonnes obligatoires
    if (!all(required_columns %in% colnames(df))) {
      shinyalert(
        title = "Colonnes manquantes",
        text = paste0(
          "Le fichier doit contenir les colonnes suivantes : ",
          paste(required_columns, collapse = ", ")
        ),
        type = "error"
      )
      return(NULL)
    }
    
    df
  })
  
  # Fixe un seuil de latence avant la modification via les sliders (optimisation)
  seuil_FC_debounced <- debounce(reactive(input$seuil_FC), 300)
  seuil_pvalue_debounced <- debounce(reactive(input$seuil_pvalue), 300)
  
  # Cette fonction réactive applique les seuils pour classifier les gènes
  processed_data <- reactive({
    if (is.null(data())) {
      return(NULL)
    }
    df <- data()
    # Fonction custom : détermine la significativité (Upregulated/Downregulated/Not significant)
    significativity(df, 
                    log2FC_cutoff = seuil_FC_debounced(),
                    P_cutoff = seuil_pvalue_debounced())
  })
  
  
  
  #==================================================================================================
  # DEG 
  #==================================================================================================
  
  # ----------------------- VOLCANO PLOT -----------------------

  # Création du volcano plot avec ggplot
  create_volcano <- reactive({
    req(data())  # S'assure que les données existent avant de créer le plot
    
    df <- processed_data()
    selected_row <- input$data_rows_selected  # Récupère la ligne sélectionnée dans le tableau
    
    # Fonction custom : génère le volcano plot avec highlight optionnel
    plot_v <- plot_volcano(df,
                           log2FC_cutoff = seuil_FC_debounced(),
                           P_cutoff = seuil_pvalue_debounced(),
                           seuil_v = input$v,      # Affiche/masque lignes verticales
                           seuil_h = input$h,      # Affiche/masque ligne horizontale
                           title = input$title_DEG,
                           highlight_row = selected_row)
    
    return(plot_v)
  })
  
  # Condition pour afficher l'image d'erreur (pas de fichier chargé)
  # Cette fonction reactive retourne TRUE si aucun fichier n'est chargé ou si data n'a pas le bon format  
  output$show_volcano_error <- reactive({
    is.null(input$file) || is.null(data())
  })
  outputOptions(output, "show_volcano_error", suspendWhenHidden = FALSE)
  
  # Image d'erreur pour le volcano plot
  output$volcano_error_img <- renderImage({
    list(
      src = "www/erreur_format.jpg",  
      contentType = "image/jpeg",      
      width = "70%",                   
      height = "auto",                 
      alt = "Format de fichier attendu"
    )
  }, deleteFile = FALSE)  
  
  # Rendu du volcano plot interactif avec Plotly
  output$volcano_plot <- renderPlotly({
    req(input$file)  
    
    # Fonction custom de création du vplot 
    volcano <- create_volcano()
    
    # Conversion du ggplot en plotly 
    ggplotly(volcano, tooltip = c("x", "y", "colour")) %>%
      layout(
        dragmode = "zoom",      
        hovermode = "closest"   
      ) %>%
      config(
        # Affiche/masque la barre d'outils selon l'input
        displayModeBar = input$toolbox_DEG,  
        modeBarButtonsToAdd = list("drawrect", "eraseshape"),
        modeBarButtonsToRemove = list("toImage"),
        displaylogo = FALSE
      ) %>%
      # WebGL boost les performances 
      plotly::toWebGL() 
  })
  
  # Téléchargement du volcano plot en PNG
  output$downloadVolcano <- downloadHandler(
    filename = function() {
      paste("Volcano_plot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = create_volcano(), width = 12, height = 8, dpi = 300)
    }
  )
  
  # ----------------------- TABLEAU DES DONNÉES -----------------------

  # Condition pour afficher le message d'erreur (pas de fichier chargé)
  output$show_data_error <- reactive({
    is.null(input$file) || is.null(data())
  })
  outputOptions(output, "show_data_error", suspendWhenHidden = FALSE)
  
  # Texte d'erreur pour le tableau
  output$data_error_text <- renderText({
    "Veuillez charger un fichier CSV au format attendu pour explorer les données"
  })
  
  # Rendu du tableau interactif avec DataTables
  output$data <- renderDT({
    req(input$file)
    
    df <- processed_data()
    
    # Affiche le tableau
    datatable(
      df,
      selection = 'single',  
      options = list(
        pageLength = 10,      
        scrollX = TRUE,       
        deferRender = TRUE,   
        scroller = TRUE   
      )
    )
  })
  

  
  #==================================================================================================
  # Enrichissement 
  #==================================================================================================
  
  # ----------------------- ORA -----------------------
  
  # Cette fonction réactive calcule les GO terms pour BP, CC, MC
  processed_data_ORA <- reactive({
    if (is.null(data())) {
      return(NULL)
    }
    
    df <- processed_data()
    
    cat("\n=== DEBUG Significance ===\n")
    cat("Valeurs uniques:", unique(df$Significance), "\n")
    cat("Nombre total de gènes:", nrow(df), "\n")
    
    significant_genes <- df$GeneName[as.character(df$Significance) != "Not significant"]
    
    if (length(significant_genes) == 0) {
      return(NULL)
    }
    
    label <- ifelse(is.null(input$title_ORA) || input$title_ORA == "", 
                    "Enrichissement", 
                    input$title_ORA)

    ego_list <- run_go_enrichment(
      gene_list = significant_genes, 
      label = label,
      org_db = org.Hs.eg.db,  # TODO: adapter selon input$species plus tard
      ontology = c("BP", "CC", "MF"),
      p_adj = "BH",
      q_cutoff = 0.05,
      key_type = "SYMBOL"
    )
    
    cat("\n=== Résultats enrichissement ===\n")
    cat("BP:", nrow(as.data.frame(ego_list$BP)), "termes\n")
    cat("CC:", nrow(as.data.frame(ego_list$CC)), "termes\n")
    cat("MF:", nrow(as.data.frame(ego_list$MF)), "termes\n")
    
    return(ego_list)
  })
  #
  #
  
  # Création des plots pour ORA
  #
  #
  create_plot_ORA <- reactive({
    req(processed_data_ORA())  # S'assure que les données existent avant de créer le plot
    
    # recuperer le Go terme souhaiter par l'utilisateur 
    ego_list <- processed_data_ORA()
    
    # selection de l'ontologie choisie (BP, CC ou MF)
    selected_go <- input$GO  
    ego <- ego_list[[selected_go]]
    # recuperer le top_n terme souhaiter par l'utilisateur 
    top_n <- input$top_n_go
    # recuperer le titre souhaiter par l'utilisateur 
    label <- ifelse(is.null(input$title_ORA) || input$title_ORA == "", 
                    paste0("Enrichissement GO-", selected_go), 
                    input$title_ORA)
    
    # Fonction custom : génère les différents plot 
    list_plot_ora <- plot_ORA(ego, label = label, top_n = top_n)
    
    return(list_plot_ora)
  })
  #
  #
  
  
  
  # Condition pour afficher l'image d'erreur (pas de fichier chargé)
  # Cette fonction reactive retourne TRUE si aucun fichier n'est chargé ou si data n'a pas le bon format  

#modifier ici pour mettre l'erreur propre a ora plot enrichissement
  output$show_enrichment_error <- reactive({
    is.null(input$file) || is.null(data()) || is.null(processed_data_ORA())
  })
  outputOptions(output, "show_enrichment_error", suspendWhenHidden = FALSE)
  
  # Image d'erreur pour l'enrichissement
  output$enrichment_error_img <- renderImage({ 
    list(
      src = "www/erreur_format.jpg",  # TODO faire une image d'erreur
      contentType = "image/jpeg",      
      width = "70%",                   
      height = "auto",                 
      alt = "Aucun gène significatif ou format du fichier incorrect"
    )
  }, deleteFile = FALSE)
  
  
  # Rendu du ora plot interactif avec Plotly
  output$enrichment_plot <- renderPlotly({
    req(create_plot_ORA())  
    
    # Fonction custom de création du vplot 
    list_ora_plot <- create_plot_ORA()
    
    #choix du plot
    selected_plot <- list_ora_plot[[input$Choosen_plot]]
    
    # Vérification que le plot existe
    if (is.null(selected_plot)) {
      return(NULL)
    }
    
    # Conversion du ggplot en plotly 
    ggplotly(selected_plot ,tooltip = c("x", "y", "text")) %>%
      layout(
        dragmode = "zoom",      
        hovermode = "closest"   
      ) %>%
      config(
        # Affiche/masque la barre d'outils selon l'input
        displayModeBar = input$toolbox_ORA,  
        modeBarButtonsToAdd = list("drawrect", "eraseshape"),
        modeBarButtonsToRemove = list("toImage"),
        displaylogo = FALSE
      ) %>%
      # WebGL boost les performances 
      plotly::toWebGL() 
  })
  
  # Téléchargement du volcano plot en PNG
  output$downloadEnrichment <- downloadHandler(
    filename = function() {
      paste0("Enrichissement_", input$Choosen_plot, "_", 
             input$GO, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      # Récupère le plot sélectionné
      list_ora_plot <- create_plot_ORA()
      selected_plot <- list_ora_plot[[input$Choosen_plot]]
      # Save le plot
      ggsave(file, plot = selected_plot, width = 12, height = 8, dpi = 300)
    }
  )

}