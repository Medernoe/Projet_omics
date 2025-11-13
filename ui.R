#==================================================================================================
#Author : Noé Méderlet 
#contact : noe.mederlet@univ-rouen.fr
#github : https://github.com/Medernoe/Projet_omics
#organism : Master Bims M2, université de rouen 
#project : Création d'une application interactive dédiée à l'analyse de données transcriptomiques,
#développée dans le cadre d'un projet universitaire du Master 2 de Bioinformatique de l'Université de Rouen.
#==================================================================================================
source("global.R")

#==================================================================================================
# INTERFACE UTILISATEUR

dashboardPage(
  #================================================================================================
  # Header
  dashboardHeader(title = ""),  # Espace vide car le titre est écrit en CSS
  
  
  #================================================================================================
  # Sidebar
  dashboardSidebar(
    
    #-----Input fichier-----
    fileInput("file", "Choisir fichier DEG :"),
    
    #-----Sélection organisme-----
    selectInput(
      inputId = "species",
      label = "Choisir une espèce :",
      choices = c("Homo sapiens", "Mus musculus", "Drosophila melanogaster"),
      selected = "Homo sapiens"
    ),
    
    #-----Menu-----
    sidebarMenu(
      menuItem("Home", tabName = "Home", icon = icon("house")),
      menuItem("DEG", tabName = "DEG", icon = icon("magnifying-glass-chart")),
      menuItem("Enrichissement", tabName = "Enrichissement", icon = icon("chart-column")),
      menuItem("Tuto", tabName = "Tuto", icon = icon("book-open")),
      menuItem("À propos", tabName = "more", icon = icon("question"))
    )
  ),
  
  #================================================================================================
  # body 
  dashboardBody(
    
    #-----CSS externe-----
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
    ),
    
    #-----Contenu des onglets-----
    tabItems(
      #==============================================================================================
      # ONGLET : HOME
      tabItem(
        tabName = "Home",
        fluidRow(
          column(
            width = 12,
            align = "center",
            
            #-----Logo-----
            tags$img(src = "logo.png", height = "250px"),
            br(), br(),
            
            #-----Titre principal-----
            h2("Bienvenue dans l'application interfacée DEGO !"),
            br(),
            
            #-----Présentation-----
            p(strong("DEGO"), "pour", strong("DEG"), "(Differentially Expressed Genes) et", 
              strong("EGO"), "(Enrichment Gene Ontology), est une application interactive dédiée à l'analyse de données transcriptomiques."),
            p("Cette application a été développée dans le cadre d'un projet universitaire du",
              strong("Master 2 de Bioinformatique de l'Université de Rouen"), "."),
            br(),
            
            #-----Informations-----
            p(strong("Réalisée par :"), "Noé Mederlet"),
            p(strong("GitHub :"), tags$a(href = "https://github.com/Medernoe/Projet_omics", "Lien vers le dépôt", target = "_blank")),
            p(strong("Contact :"), tags$a(href = "noe.mederlet@univ-rouen.fr", "noe.mederlet@univ-rouen.fr")),
            br(),
            
            p("Pour plus d'informations sur l'utilisation de l'application, consultez la section ", strong("« Tuto »"), "."),
            p("Pour en savoir davantage sur le projet, rendez-vous dans la section ", strong("« À propos »"), ".")
          )
        )
      ),
      
      
      #==============================================================================================
      # ONGLET : DEG
      tabItem(
        tabName = "DEG",
        
        #--------------------------------------------------------------------------------------------
        # Volcano Plot + Paramètres seuils 
        fluidRow(
          
          #-----Colonne gauche : Volcano Plot-----
          column(
            width = 8,
            box(
              title = "Volcano Plot",
              width = 12,
              
              # Image d'erreur conditionnelle (affichée si pas de fichier)
              conditionalPanel(
                condition = "output.show_volcano_error",
                div(
                  style = "text-align: center; padding: 20px;",
                  p("Veuillez charger un fichier CSV au format attendu pour visualiser son Volcano Plot",
                    style = "font-size: 16px; color: #666; margin-bottom: 20px;"),
                  imageOutput("volcano_error_img", height = "450px", inline = TRUE)
                )
              ),
              
              # Plot conditionnel (affiché si fichier chargé)
              conditionalPanel(
                condition = "!output.show_volcano_error",
                plotlyOutput("volcano_plot", height = "500px")
              ),
              
              # Bouton de téléchargement
              downloadButton("downloadVolcano", "Télécharger")
            )
          ),
          
          #-----Colonne droite : Paramètres seuils-----
          column(
            width = 4,
            
            # Box : Titre du graphique
            box(
              title = "Titre",
              width = 12,
              textInput("Titre", "Entrer un titre pour la figure :")
            ),
            
            # Box : Activation ToolBox
            box(
              title = "ToolBox",
              width = 12,
              checkboxInput("toolbox", "Activer la barre d'outils", value = TRUE)
            ),
            
            # Box : Seuils de significativité
            box(
              title = "Seuils",
              width = 12,
              # Seuil Log2 Fold Change
              sliderInput("seuil_FC", "Log2 FC :", min = 0, max = 5, value = 1, step = 0.5, width = "100%"),
              checkboxInput("v", "Activer ligne verticale (logFC)", value = TRUE),
              # Seuil P-value
              sliderInput("seuil_pvalue", "P-value :", min = 0, max = 1, value = 0.05, width = "100%"),
              checkboxInput("h", "Activer ligne horizontale (p-value)", value = TRUE)
            )
          )
        ),
        
        #--------------------------------------------------------------------------------------------
        # Tableau des données
        fluidRow(
          box(
            title = "Tableau des données",
            width = 12,
            
            # Texte d'erreur conditionnel (affiché si pas de fichier)
            conditionalPanel(
              condition = "output.show_data_error",
              div(
                style = "text-align: center; padding: 50px;",
                h3(textOutput("data_error_text"), style = "color: #666;")
              )
            ),
            
            # Tableau conditionnel (affiché si fichier chargé)
            conditionalPanel(
              condition = "!output.show_data_error",
              DTOutput("data")
            )
          )
        )
      ),
      
      
      #==============================================================================================
      # ONGLET : ENRICHISSEMENT (vide)
      tabItem(
        tabName = "Enrichissement",
        h2("Coming soon")
      ),
      
      #==============================================================================================
      # ONGLET : TUTO (vide)
      tabItem(
        tabName = "Tuto",
        h2("Coming soon")
      ),
      
      #==============================================================================================
      # ONGLET : À PROPOS
      tabItem(
        tabName = "more",
        fluidRow(
          box(
            title = "À propos du projet",
            width = 12,
            # Charger le contenu HTML externe
            includeHTML("www/about.html")
          )
        )
      )
    )
  )
)