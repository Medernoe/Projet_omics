#==================================================================================================
#Author : Noé Méderlet 
#contact : noe.mederlet@gmail.com
#github : 
#organism : Master Bims M2, université de rouen 
#project : 


#==================================================================================================
library(ggplot2)
library(DT)
source("fonctions.R")


# Serveur
server <- function(input, output) {
  
  #data
  data <- reactive({read.csv(input$file$datapath)})
  
  
  output$volcano_plot <- renderPlot({
    
    #transforme l'objet réactive de data en df 
    df <- data()
    
    df <- significativity(df, 
                          log2FC_cutoff = input$seuil_FC,
                          P_cutoff = input$seuil_pvalue)
    
    #fonction custom : script fonctions.R
    plot_volcano(df,
                 log2FC_cutoff = input$seuil_FC,
                 P_cutoff = input$seuil_pvalue,
                 seuil_v = input$v,
                 seuil_h = input$h,
                 title = "volcano plot")
  })
  
  output$data <- renderDataTable({
    df <- data()
    
    df <- significativity(df, 
                          log2FC_cutoff = input$seuil_FC,
                          P_cutoff = input$seuil_pvalue)
    
    datatable(df)
  })
}








