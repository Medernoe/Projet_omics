#==================================================================================================
#Author : Noé Méderlet 
#contact : noe.mederlet@univ-rouen.fr
#github : https://github.com/Medernoe/Projet_omics
#organism : Master Bims M2, université de rouen 
#project : Création d'une application interactive dédiée à l’analyse de données transcriptomiques,
#développée dans le cadre d’un projet universitaire du Master 2 de Bioinformatique de l’Université de Rouen.
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
    # Injecter titre stylisé en HTML
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
          display: flex;              
          align-items: center;
        }
      ")
      )
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
            p(strong("GitHub :"), tags$a(href = "https://github.com/Medernoe/Projet_omics", "Lien vers le dépôt", target = "_blank")),
            p(strong("Contact :"), tags$a(href = "noe.mederlet@univ-rouen.fr", "noe.mederlet@univ-rouen.fr")),
            br(),
            
            # Sections d'information
            p("Pour plus d’informations sur l’utilisation de l’application, consultez la section ", strong("« Tuto »"), "."),
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
              plotOutput("volcano_plot", height = "500px"),
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
            
            #box pour fixer les seuils pvalue et logFC + checkbox pour l'ajout de trait fixé au seuil
            box(
              title = "Seuils",
              width = 12,
              sliderInput("seuil_FC", "Log2 FC :", min = 0, max = 5, value = 1, step = 0.5, width = "100%"),
              checkboxInput("v", "Activer ligne verticale (logFC)", value = TRUE),
              sliderInput("seuil_pvalue", "P-value :", min = 0, max = 1, value = 0.5, width = "100%"),
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
            # Utilisation de HTML
            HTML('
            <div class="panel-group" id="faqAccordion">
    
              <!-- Question 1 -->
              <div class="panel panel-default">
                <div class="panel-heading" style="background-color:#001f3f; color:white; cursor:pointer;" 
                      data-toggle="collapse" data-parent="#faqAccordion" href="#q1">
                  <h4 class="panel-title">Objectif ?</h4>
                </div>
                <div id="q1" class="panel-collapse collapse">
                  <div class="panel-body">
                    Ce projet est proposé dans le cadre d’un enseignement universitaire. 
                    Il a pour vocation l’approfondissement des compétences et de la compréhension 
                    des langages R et RShiny. L’application a pour but de permettre la 
                    génération de graphiques d’annotation fonctionnelle à partir de matrices 
                    en entrée contenant des valeurs de LogFC.
                  </div>
                </div>
              </div>
    
    
              <!-- Question 2 -->
              <div class="panel panel-default">
                <div class="panel-heading" style="background-color:#001f3f; color:white; cursor:pointer;" 
                      data-toggle="collapse" data-parent="#faqAccordion" href="#q2">
                  <h4 class="panel-title">Pourquoi un logo moche ?</h4>
                </div>
                <div id="q2" class="panel-collapse collapse">
                  <div class="panel-body">
                    Il est important de noter que je suis bioinformaticien et non graphiste. 
                    De plus le premier avis d’un utilisateur sur la qualité d’un logiciel vient du premier regard. 
                    Avec un logo un peu daté, les utilisateurs penseront peut-être qu’il s’agit 
                    d’une vieille interface mais paradoxalement, cela pourrait aussi leur inspirer 
                    confiance en la robustesse du projet. (C\'était mieux avant !)
                  </div>
                </div>
              </div>
    
    
              <!-- Question 3 -->
              <div class="panel panel-default">
                <div class="panel-heading" style="background-color:#001f3f; color:white; cursor:pointer;" 
                      data-toggle="collapse" data-parent="#faqAccordion" href="#q3">
                  <h4 class="panel-title">À venir ?</h4>
                </div>
                <div id="q3" class="panel-collapse collapse">
                  <div class="panel-body">
                    Réécriture des algorithmes d’enrichissement pour proposer des approches 
                    réactives aux jeux de données, amélioration globale de l’interface, 
                    ajout de nouveaux modules dans l’onglet « Enrichissement ». 
                    Un tutoriel et une documentation détaillée seront également ajoutés 
                    pour accompagner les utilisateurs et assurer un suivi du projet.
                  </div>
                </div>
              </div>
              
            </div>
          ')
          )
        )
      ) 
    )
  )
)

