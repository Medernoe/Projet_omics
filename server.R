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
  # VOLCANO PLOT

  # Création du volcano plot avec ggplot
  create_volcano <- reactive({
    req(processed_data())  # S'assure que les données existent avant de créer le plot
    
    df <- processed_data()
    selected_row <- input$data_rows_selected  # Récupère la ligne sélectionnée dans le tableau
    
    # Fonction custom : génère le volcano plot avec highlight optionnel
    plot_v <- plot_volcano(df,
                           log2FC_cutoff = seuil_FC_debounced(),
                           P_cutoff = seuil_pvalue_debounced(),
                           seuil_v = input$v,      # Affiche/masque lignes verticales
                           seuil_h = input$h,      # Affiche/masque ligne horizontale
                           title = input$Titre,
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
    p <- create_volcano()
    
    # Conversion du ggplot en plotly 
    ggplotly(p, tooltip = c("x", "y", "colour")) %>%
      layout(
        dragmode = "zoom",      
        hovermode = "closest"   
      ) %>%
      config(
        # Affiche/masque la barre d'outils selon l'input
        displayModeBar = input$toolbox,  
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
  
  #==================================================================================================
  # TABLEAU DES DONNÉES

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
}