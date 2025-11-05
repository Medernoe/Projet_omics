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
  # Fonction réactive qui crée le plot
  create_volcano <- reactive({
    req(data())  # S'assure que les données existent
    
    df <- data()
    
    df <- significativity(df, 
                          log2FC_cutoff = input$seuil_FC,
                          P_cutoff = input$seuil_pvalue)
    
    plot_v <- plot_volcano(df,
                           log2FC_cutoff = input$seuil_FC,
                           P_cutoff = input$seuil_pvalue,
                           seuil_v = input$v,
                           seuil_h = input$h,
                           title = input$Titre)
    
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
      # Utiliser la fonction réactive
      ggsave(file, plot = create_volcano(), width = 12, height = 8, dpi = 300)
    }
  )
  
  #==============================================================================  
  #sortie tableau
  output$data <- renderDataTable({
    df <- data()
    
    #fonction custom : script fonctions.R => seuil de significativité 
    df <- significativity(df, 
                          log2FC_cutoff = input$seuil_FC,
                          P_cutoff = input$seuil_pvalue)
    #affiche la table
    datatable(df)
  })
}