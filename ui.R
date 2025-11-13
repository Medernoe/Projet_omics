#==================================================================================================
#Author : Noé Méderlet 
#contact : noe.mederlet@univ-rouen.fr
#github : https://github.com/Medernoe/Projet_omics
#organism : Master Bims M2, université de rouen 
#project : Création d'une application interactive dédiée à l'analyse de données transcriptomiques,
#développée dans le cadre d'un projet universitaire du Master 2 de Bioinformatique de l'Université de Rouen.
#==================================================================================================
source("global.R")

ui <- dashboardPage(
  
  #Entete
  dashboardHeader(title = ""),  # Espace vide car le titre est ecrit en HTML CSS
  
  
  #Sidebar   
  dashboardSidebar(
    
    #boite de selection du fichier d'entrée
    fileInput("file", "Choisir fichier DEG :"),
    
    #boite de selection de l'organisme
    selectInput(
      inputId = "species",
      label = "Choisir une espèce :",
      choices = c("Homo sapiens", "Mus musculus", "Drosophila melanogaster"),
      selected = "Homo sapiens"
    ),
    
    #Sous menu avec logo
    sidebarMenu(
      menuItem("Home", tabName = "Home", icon = icon("house")),
      menuItem("DEG", tabName = "DEG", icon = icon("magnifying-glass-chart")),
      menuItem("Enrichissement", tabName = "Enrichissement", icon = icon("chart-column")),
      menuItem("Tuto", tabName = "Tuto", icon = icon("book-open")),
      menuItem("À propos", tabName = "more", icon = icon("question"))
    )
  ),
  
  #Corps
  dashboardBody(
    # Injecter le fichier CSS externe
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
    ),
    
    #Remplie les champs définis dans sidebarMenu()
    tabItems(
      
      # Onglet Home 
      tabItem(
        tabName = "Home",
        fluidRow(
          column(
            width = 12,
            align = "center",
            
            # Logo
            tags$img(src = "logo.png", height = "250px"),
            br(), br(),
            
            # Titre principal
            h2("Bienvenue dans l'application interfacée DEGO !"),
            br(),
            
            # Texte de présentation
            p(strong("DEGO"), "pour", strong("DEG"), "(Differentially Expressed Genes) et", 
              strong("EGO"), "(Enrichment Gene Ontology), est une application interactive dédiée à l'analyse de données transcriptomiques."),
            p("Cette application a été développée dans le cadre d'un projet universitaire du",
              strong("Master 2 de Bioinformatique de l'Université de Rouen"), "."),
            br(),
            
            # Informations sur le créateur
            p(strong("Réalisée par :"), "Noé Mederlet"),
            p(strong("GitHub :"), tags$a(href = "https://github.com/Medernoe/Projet_omics", "Lien vers le dépôt", target = "_blank")),
            p(strong("Contact :"), tags$a(href = "noe.mederlet@univ-rouen.fr", "noe.mederlet@univ-rouen.fr")),
            br(),
            
            # Sections d'information
            p("Pour plus d'informations sur l'utilisation de l'application, consultez la section ", strong("« Tuto »"), "."),
            p("Pour en savoir davantage sur le projet, rendez-vous dans la section ", strong("« À propos »"), ".")
          )
        )
      ),
      
      
      # Onglet DEG
      tabItem(
        tabName = "DEG",
        
        #Premiere ligne de la page 
        fluidRow(
          # Colonne gauche : le graphique
          column(
            width = 8,
            box(
              title = "Volcano Plot",
              width = 12,
              plotlyOutput("volcano_plot", height = "500px"),
              downloadButton("downloadVolcano", "Télécharger")
            )
          ),
          
          # Colonne droite : paramètres et titre
          column(
            width = 4,
            
            #Box demandant un titre pour la figure Volcano plot 
            box(
              title = "Titre",
              width = 12,
              textInput("Titre", "Entrer un titre pour la figure :")
            ),
            
            #box d'activation et de désactivation de la barre d'outils
            box(
              title = "ToolBox",
              width = 12,
              checkboxInput("toolbox", "Activer la barre d'outils", value = TRUE),
            ),
            
            #box pour fixer les seuils pvalue et logFC + checkbox pour l'ajout de trait fixé au seuil
            box(
              title = "Seuils",
              width = 12,
              #modifier ici max = 5 par max(log2FC)
              sliderInput("seuil_FC", "Log2 FC :", min = 0, max = 5, value = 1, step = 0.5, width = "100%"),
              checkboxInput("v", "Activer ligne verticale (logFC)", value = TRUE),
              sliderInput("seuil_pvalue", "P-value :", min = 0, max = 1, value = 0.05, width = "100%"),
              checkboxInput("h", "Activer ligne horizontale (p-value)", value = TRUE)
            )
          )
        ),
        
        #Seconde ligne sur la page 
        fluidRow(
          box(
            title = "Tableau des données",
            width = 12,
            DTOutput("data")  
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
      
      # Onglet A propos
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