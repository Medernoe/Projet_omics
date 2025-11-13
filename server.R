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
server <- function(input, output) {
  #==============================================================================    
  # Dataframe des pvalue et logFC
  data <- reactive({
    if (is.null(input$file)) {
      return(NULL)
    }
    read.csv(input$file$datapath, sep=";")
  })
  
  #==============================================================================  
  # Données traitées 
  processed_data <- reactive({
    if (is.null(data())) {
      return(NULL)
    }
    df <- data()
    significativity(df, 
                    log2FC_cutoff = input$seuil_FC,
                    P_cutoff = input$seuil_pvalue)
  })
  
  #==============================================================================  
  # Fonction réactive qui crée le plot
  create_volcano <- reactive({
    req(processed_data()) 
    
    df <- processed_data()
    selected_row <- input$data_rows_selected 
    
    # Utilise la fonction custom : plot_volcano
    plot_v <- plot_volcano(df,
                           log2FC_cutoff = input$seuil_FC,
                           P_cutoff = input$seuil_pvalue,
                           seuil_v = input$v,
                           seuil_h = input$h,
                           title = input$Titre,
                           highlight_row = selected_row)  
    
    return(plot_v)
  })
  
  #==============================================================================  
  # Sortie volcanoplot avec plotly pour interactivité
  output$volcano_plot <- renderPlotly({
    # Vérifier si le fichier est chargé
    validate(
      need(!is.null(input$file), 
           "Veuillez charger un fichier CSV au format attendu pour visualiser son Volcano Plot")
    )
    
    p <- create_volcano()
    
    # Convertir le plot en plotly
    ggplotly(p, tooltip = c("x", "y", "colour")) %>%
      layout(
        dragmode = "zoom",
        hovermode = "closest"
      ) %>%
      config(
        displayModeBar = input$toolbox, #permet de selectionner ou pas la toolbox
        modeBarButtonsToAdd = list("drawrect", "eraseshape"),
        modeBarButtonsToRemove = list("toImage"),
        displaylogo = FALSE
      )
  })
  
  #==============================================================================  
  # Téléchargement du volcano plot
  output$downloadVolcano <- downloadHandler(
    filename = function() {
      paste("Volcano_plot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = create_volcano(), width = 12, height = 8, dpi = 300)
    }
  )
  
  #==============================================================================  
  # Sortie tableau
  output$data <- renderDT({
    # Vérifier si le fichier est chargé
    validate(
      need(!is.null(input$file), 
           "Veuillez charger un fichier CSV au format attendu pour afficher le tableau des données")
    )
    
    df <- processed_data()
    
    # Affiche la table avec sélection activée
    datatable(
      df,
      selection = 'single',
      options = list(
        pageLength = 10,
        scrollX = TRUE
      )
    )
  })
}