source("global.R")

ui <- dashboardPage(
  #Entete
  dashboardHeader(title = ""),  # On met un espace vide pour éviter l'erreur


  #Sidebar   
  dashboardSidebar(
    fileInput("file", "Choisir fichier DEG :"),
    
    selectInput(
      inputId = "species",
      label = "Choisir une espèce :",
      choices = c("Homo sapiens", "Mus musculus", "Drosophila melanogaster"),
      selected = "Homo sapiens"
    ),
    
    sidebarMenu(
      menuItem("Home", tabName = "Home", icon = icon("house")),
      menuItem("DEG", tabName = "DEG", icon = icon("magnifying-glass-chart")),
      menuItem("Enrichissement", tabName = "Enrichissement", icon = icon("chart-column")),
      menuItem("Tuto", tabName = "Tuto", icon = icon("book-open")),
      menuItem("À propos", tabName = "aide", icon = icon("question"))
    )
  ),
  
  dashboardBody(
    # Injecter notre titre stylisé
    tags$head(
      tags$style(HTML("
        /* Masquer le titre de base */
        .main-header .logo {
          visibility: hidden;
        }
        
        /* Changer le fond de l'en-tête */
        .main-header {
          background-color: #001f3f;  /* bleu foncé */
          
        }

        /* Ajouter un nouveau titre personnalisé */
        .main-header::before {
          content: 'DƎGO';
          position: absolute;
          left: 45px;
          top: 1px;
          bottom: 1px;
          font-family: 'Courier New', monospace;
          font-size: 55px;
          font-weight: bold;
          color: white;
          display: flex;              /* permet d’ajuster verticalement le texte */
          align-items: center;        /* centre verticalement dans l’espace défini par top/bottom */
        }
      ")
      )
    ),
    
    
    tabItems(
      # Onglet Home 
      tabItem(
        tabName = "Home",
        fluidRow(
          column(
            width = 12,
            align = "center",
            
            # Logo
            tags$img(src = "logo.png", height = "300px"),
            br(), br(),
            
            # Titre principal
            h2("Bienvenue dans l’application interfacée DEGO !"),
            br(),
            
            # Texte de présentation
            p(strong("DEGO"), "pour", strong("DEG"), "(Differentially Expressed Genes) et", 
              strong("EGO"), "(Enrichment Gene Ontology), est une application interactive dédiée à l’analyse de données transcriptomiques."),
            p("Cette application a été développée dans le cadre d’un projet universitaire du",
              strong("Master 2 de Bioinformatique de l’Université de Rouen"), "."),
            br(),
            
            # Informations sur le créateur
            p(strong("Réalisée par :"), "Noé Mederlet"),
            p(strong("GitHub :"), tags$a(href = "https://github.com/tonLien", "Lien vers le dépôt", target = "_blank")),
            p(strong("Contact :"), tags$a(href = "mailto:ton.email@exemple.com", "ton.email@exemple.com")),
            br(),
            
            # Sections d'information
            p("Pour plus d’informations sur l’utilisation de l’application, consultez la section ", strong("« Tuto »"), "."),
            p("Pour en savoir davantage sur le projet, rendez-vous dans la section ", strong("« À propos »"), ".")
          )
        )
      ),
      
      
      
      # Onglet DEG
      tabItem(tabName = "DEG",
              fluidRow(
                box(plotOutput("volcano_plot"), width = 8,
                    downloadButton("downloadVolcano", "Télécharger")),
                
                box(
                  width = 4,
                  title = 'Titre', 
                  textInput("Title","Entrer un titre pour la figure :")
                ),
                
                box(
                  title = "Seuils",
                  width = 4,
                  sliderInput("seuil_FC", "Log2 FC :", min = 0, max = 5, value = 1, step = 0.5, width = "100%"),
                  checkboxInput("v", "Activer ligne verticale (logFC)", value = TRUE),
                  sliderInput("seuil_pvalue", "P-value :", min = 0, max = 1, value = 0.5, width = "100%"),
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
              h2("Coming soon")
      ),
      
      # Onglet Tuto
      tabItem(tabName = "Tuto",
              h2("Coming soon")
      ),
              
      # Onglet Aide
      tabItem(tabName = "aide",
              h3("Aide / Documentation")
      )
    )
  )
)

