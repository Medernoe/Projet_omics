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
  
  #dataframe des pvalue et logFC
  data <- reactive({read.csv(input$file$datapath)})
  
  #sortie volcanoplot
  output$volcano_plot <- renderPlot({
    
    #transforme l'objet réactive de data en df 
    df <- data()
    
    #fonction custom : script fonctions.R => seuil de significativité 
    df <- significativity(df, 
                          log2FC_cutoff = input$seuil_FC,
                          P_cutoff = input$seuil_pvalue)
    
    #fonction custom : script fonctions.R => produit un V plot 
    plot_volcano(df,
                 log2FC_cutoff = input$seuil_FC,
                 P_cutoff = input$seuil_pvalue,
                 seuil_v = input$v,
                 seuil_h = input$h,
                 title = input$Titre)
  })
  
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