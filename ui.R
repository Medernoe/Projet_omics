library(shiny)
library(shinydashboard)
library(ggplot2)

ui <- dashboardPage(
  dashboardHeader(title = "Volcano Plot App"),
  
  dashboardSidebar(
    fileInput("file", "Choisir fichier DEG :"),
    sidebarMenu(
      menuItem("Home", tabName = "Home", icon = icon("house")),
      menuItem("DEG", tabName = "DEG", icon = icon("magnifying-glass-chart")),
      menuItem("Enrichissement", tabName = "Enrichissement", icon = icon("chart-column")),
      menuItem("Aide", tabName = "aide", icon = icon("question"))
    )
  ),
  
  dashboardBody(
    tabItems(
      # Onglet Home (tu peux le laisser vide ou mettre du contenu)
      tabItem(tabName = "Home",
              h2("Bienvenue")
      ),
      
      # Onglet DEG
      tabItem(tabName = "DEG",
              fluidRow(
                box(plotOutput("volcano_plot"), width = 8,
                    downloadButton("downloadVolcano", "Télécharger")),
                box(
                  title = "Seuils",
                  width = 4,
                  sliderInput("seuil_FC", "Log2 FC :", min = 0, max = 5, value = 1, step = 0.5, width = "100%"),
                  checkboxInput("v", "Activer ligne verticale (logFC)", value = TRUE),
                  sliderInput("seuil_pvalue", "P-value :", min = 0, max = 1, value = 1, width = "100%"),
                  checkboxInput("h", "Activer ligne horizontale (p-value)", value = TRUE)
                )
              ),
              fluidRow(
                box(
                  title = "Tableau des données",
                  width = 12,
                  dataTableOutput("data")
                )
              )
      ),
      
      # Onglet Enrichissement
      tabItem(tabName = "Enrichissement",
              h3("Contenu pour Enrichissement")
      ),
      
      # Onglet Aide
      tabItem(tabName = "aide",
              h3("Aide / Documentation")
      )
    )
  )
)

