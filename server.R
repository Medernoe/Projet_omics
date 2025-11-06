#==================================================================================================
#Author : Noé Méderlet 
#contact : noe.mederlet@univ-rouen.fr
#github : https://github.com/Medernoe/Projet_omics
#organism : Master Bims M2, université de rouen 
#project : Création d'une application interactive dédiée à l’analyse de données transcriptomiques,
#développée dans le cadre d’un projet universitaire du Master 2 de Bioinformatique de l’Université de Rouen.
#==================================================================================================
source("global.R")

# Serveur
server <- function(input, output) {
  #==============================================================================    
  #dataframe des pvalue et logFC
  data <- reactive({read.csv(input$file$datapath)})
  
  #==============================================================================  
  # Données traitées 
  processed_data <- reactive({
    req(data())
    df <- data()
    significativity(df, 
                    log2FC_cutoff = input$seuil_FC,
                    P_cutoff = input$seuil_pvalue)
  })
  
  #==============================================================================  
  # Fonction réactive qui crée le plot
  create_volcano <- reactive({
    # assure que les données existent
    req(processed_data())  
    
    df <- processed_data()
    # Récupère la ligne sélectionnée
    selected_row <- input$data_rows_selected 
    
    # DEBUG : Afficher dans la console R
    cat("Selected row:", selected_row, "\n")
    cat("Length:", length(selected_row), "\n")
    cat("Is NULL:", is.null(selected_row), "\n")
    
    #utilise la fonction custom : plot_volcano
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
  #sortie volcanoplot
  output$volcano_plot <- renderPlot({
    create_volcano()
  })
  
  #==============================================================================  
  #Téléchargement du volcano plot
  output$downloadVolcano <- downloadHandler(
    filename = function() {
      #creer le nom 
      paste("Volcano_plot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      # Utiliser la fonction réactive pour recreer le plot 
      ggsave(file, plot = create_volcano(), width = 12, height = 8, dpi = 300)
    }
  )
  
  #==============================================================================  
  #sortie tableau
  output$data <- renderDT({
    req(processed_data())
    
    df <- processed_data()
    
    #affiche la table avec sélection activée
    datatable(
      df,
      selection = 'single',  # Permet de sélectionner UNE seule ligne, à modifier pour prendre plusieur lignes
      options = list(
        pageLength = 10,
        scrollX = TRUE
      )
    )
  })
}