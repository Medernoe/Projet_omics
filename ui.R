#===========================Presentation========================================
#Authors : Noé Méderlet, Mehdi Tachekort, Mathieu Cartier, Valentin Fourdreane
#contact : noe.mederlet@univ-rouen.fr, mehdi.tachekort@univ-rouen.fr, mathieu.cartier@univ-rouen.fr, valentin.fourdraine@univ-rouen.fr
#github : https://github.com/Medernoe/Projet_omics
#organism : Master Bims M2, université de rouen 
#project : Création d'une application interactive dédiée à l'analyse de données transcriptomiques,
#développée dans le cadre d'un projet universitaire du Master 2 de Bioinformatique de l'Université de Rouen.
#==============================================================================#

#=========Chargement du global.R et de shinyalter===============================
source("global.R")
useShinyalert()

#=========================INTERFACE UTILISATEUR=================================

#=============================Dashboard=========================================

dashboardPage(
  
  ####============================Header==========================================
  
  dashboardHeader(title = ""),  # Espace vide car le titre est écrit en CSS
  
  ####===========================Sidebar==========================================
  dashboardSidebar(
    
    #####=======================Input fichier=======================================
    fileInput("deg_file", "Choisir fichier DEG :"),
    
    #####======================Sélection organisme==================================
    selectInput(
      inputId = "species",
      label = "Choisir une espèce :",
      choices = c("Homo sapiens", "Mus musculus", "Drosophila melanogaster"),
      selected = "Homo sapiens"
    ),
    
    #####===========================Menu============================================
    sidebarMenu(
      menuItem("Accueil", tabName = "home", icon = icon("house")),
      menuItem("DEG", tabName = "deg", icon = icon("magnifying-glass-chart")),
      menuItem("Enrichissement", tabName = "enrichment", icon = icon("chart-column")),
      menuItem("Documentation", tabName = "documentation", icon = icon("book-open")),
      menuItem("À propos", tabName = "about", icon = icon("question"))
    )
  ),
  
  #===============================Body============================================
  dashboardBody(
    
    ####==========================CSS externe=======================================
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
    ),
    
    ####=====================Contenu des onglets====================================
    tabItems(
      #####===========================Accueil=======================================
      tabItem(
        tabName = "home",
        fluidRow(
          column(
            width = 12,
            align = "center",
            
            ######==========================Logo============================================
            tags$img(src = "logo.png", height = "250px"),
            br(), br(),
            
            ######=====================Titre principal======================================
            
            h2("Bienvenue dans l'application interfacée DEGO !"),
            br(),
            
            ######=====================Présentation=========================================
            p(strong("DEGO"), "pour", strong("DEG"), "(Differentially Expressed Genes) et", 
              strong("EGO"), "(Enrichment Gene Ontology), est une application interactive dédiée à l'analyse de données transcriptomiques."),
            p("Cette application a été développée dans le cadre d'un projet universitaire du",
              strong("Master 2 de Bioinformatique de l'Université de Rouen"), "."),
            br(),
            
            ######=====================Informations======================================
            p(strong("Réalisée par :"), "Noé Mederlet"),
            p(strong("GitHub :"), tags$a(href = "https://github.com/Medernoe/Projet_omics", "Lien vers le dépôt", target = "_blank")),
            p(strong("Contact :"), tags$a(href = "noe.mederlet@univ-rouen.fr", "noe.mederlet@univ-rouen.fr")),
            br(),
            
            p("Pour plus d'informations sur l'utilisation de l'application, consultez la section ", strong("« Tuto »"), "."),
            p("Pour en savoir davantage sur le projet, rendez-vous dans la section ", strong("« À propos »"), ".")
          )
        )
      ),
      
      
      #####==========================DEG==============================================
      tabItem(
        tabName = "deg",
        
        ######=====================Volcano Plot=========================================
        fluidRow(
          
          ########-----Colonne gauche : Volcano Plot-----
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
                withSpinner(plotlyOutput("volcano_plot", height = "500px"),
                            type = 0,        # Type de spinner (0 for custom 1 à 8)
                            color = "#3498db", # Couleur (bleu par défaut)
                            image = "sleepy-snorlax.gif",       
                            image.width = "50%",      
                            image.height = "50%",
                )
              ),
              
              # Bouton de téléchargement
              downloadButton("downloadVolcano", "Télécharger")
            )
          ),
          
          ########-----Colonne droite : Paramètres seuils-----
          column(
            width = 4,
            
            # Box : Titre du graphique
            box(
              title = "Titre",
              width = 12,
              textInput("deg_title", "Entrer un titre pour la figure :"),
            ),
            
            # Box : Activation ToolBox
            box(
              title = "ToolBox",
              width = 12,
              checkboxInput("deg_toolbox", "Activer la barre d'outils", value = TRUE),
            ),
            
            # Box : Seuils de significativité
            box(
              title = "Seuils",
              width = 12,
              # Seuil Log2 Fold Change
              sliderInput("fc_threshold", "Log2 FC :", min = 0, max = 5, value = 1, step = 0.5, width = "100%"),
              checkboxInput("show_vline", "Activer ligne verticale (logFC)", value = TRUE),
              # Seuil P-value
              sliderInput("pvalue_threshold", "P-value :", min = 0, max = 1, value = 0.05, width = "100%"),
              checkboxInput("show_hline", "Activer ligne horizontale (p-value)", value = TRUE)
            )
          )
        ),
        
        ######=========================Tableau de données===============================
        # Tableau des données
        fluidRow(
          box(
            title = "Tableau des données",
            width = 12,
            
            # Texte d'erreur conditionnel (affiché si pas de fichier)
            conditionalPanel(
              condition = "output.show_table_error",
              div(
                style = "text-align: center; padding: 50px;",
                h3(textOutput("table_error_text"), style = "color: #666;")
              )
            ),
            
            # Tableau conditionnel (affiché si fichier chargé)
            conditionalPanel(
              condition = "!output.show_table_error",
              DTOutput("deg_table")
            )
          )
        )
      ),
      
      #####==========================ENRICHISSEMENT===================================
      tabItem(
        tabName = "enrichment",
        
        ######-----Bandeau du haut : 3 box de contrôle d'analyse-----
        fluidRow(
          # Box 1 : Méthode d'enrichissement
          column(
            width = 4,
            box(
              title = "Méthode d'enrichissement",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              radioButtons(
                inputId = "enrichment_method",
                label = NULL,
                choices = c("ORA", "GSEA"),
                selected = "ORA",
                inline = TRUE
              )
            )
          ),
          
          # Box 2 : Ontologie GO (sélection multiple)
          column(
            width = 4,
            box(
              title = "Ontologie GO",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              checkboxGroupInput(
                inputId = "go_ontology",
                label = NULL,
                choices = c("Processus Biologique (BP)" = "BP",
                            "Composant Cellulaire (CC)" = "CC",
                            "Fonction Moléculaire (MF)" = "MF"),
                selected = "BP",
                inline = FALSE
              )
            )
          ),
          
          # Box 3 : Bouton de lancement
          column(
            width = 4,
            box(
              title = "Lancer l'analyse",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              div(
                style = "text-align: center; padding: 10px;",
                actionButton(
                  inputId = "run_enrichment",
                  label = "Lancer l'enrichissement",
                  icon = icon("play"),
                  class = "btn-success btn-lg",
                  style = "width: 100%;"
                ),
                br(), br(),
                # Statut affiché dynamiquement (à brancher côté server plus tard)
                div(
                  style = "font-size: 12px; color: #666; font-style: italic;",
                  textOutput("enrichment_status")
                )
              )
            )
          ),
          ######-----Plot + Paramètres de visualisation-----
          fluidRow(
            
            ########-----Colonne gauche (8/12) : Plot + contrôles d'affichage-----
            column(
              width = 8,
              
              # Box du plot
              box(
                title = "Enrichissement Plot",
                width = 12,
                
                # Image d'erreur conditionnelle (affichée si pas de fichier ou pas encore lancé)
                conditionalPanel(
                  condition = "output.show_enrichment_error",
                  div(
                    style = "text-align: center; padding: 20px;",
                    p("Veuillez charger un fichier CSV puis cliquer sur « Lancer l'enrichissement »",
                      style = "font-size: 16px; color: #666; margin-bottom: 20px;"),
                    imageOutput("enrichment_error_img", height = "450px", inline = TRUE)
                  )
                ),
                
                # Plot conditionnel
                conditionalPanel(
                  condition = "!output.show_enrichment_error",
                  withSpinner(plotlyOutput("enrichment_plot", height = "500px"),
                              type = 0,
                              image = "sleepy-snorlax.gif",
                              image.width = "50%",
                              image.height = "50%")
                ),
                
                downloadButton("downloadEnrichment", "Télécharger")
              ),
            ),
            
            ########-----Colonne droite (4/12) : Box unique avec les 3 sélections-----
            column(
              width = 4,
              
              box(
                title = "Paramètres du graphique",
                width = 12,
                
                # 1. Type de visualisation (dépend de la méthode ORA/GSEA)
                conditionalPanel(
                  condition = "input.enrichment_method == 'ORA'",
                  selectInput(
                    inputId = "selected_plot_ora",
                    label = "Type de visualisation :",
                    choices = c("Dotplot", "Barplot", "Cnetplot", "Emapplot",
                                "Goplot", "Upsetplot", "Heatplot"),
                    selected = "Dotplot"
                  )
                ),
                conditionalPanel(
                  condition = "input.enrichment_method == 'GSEA'",
                  selectInput(
                    inputId = "selected_plot_gsea",
                    label = "Type de visualisation :",
                    choices = c("Dotplot", "Cnetplot", "Emapplot", "Upsetplot",
                                "Heatplot", "Ridgeplot", "GSEAplot2", "GSEArank"),
                    selected = "Dotplot"
                  )
                ),
                
                # 2. Ontologie à afficher (mise à jour côté server après calcul)
                selectInput(
                  inputId = "displayed_ontology",
                  label = "Ontologie à afficher :",
                  choices = c("Processus Biologique (BP)" = "BP",
                              "Composant Cellulaire (CC)" = "CC",
                              "Fonction Moléculaire (MF)" = "MF"),
                  selected = "BP"
                ),
                
                # 3. Nombre de termes GO
                sliderInput(
                  inputId = "top_n_terms",
                  label = "Nombre de termes GO à afficher :",
                  min = 5, max = 50, value = 10, step = 5, width = "100%"
                )
              )
            ),
            column(
              width = 4,
              # Box d'affichage : titre et toolbox sous le plot
              box(
                title = "Affichage",
                width = 12,
                textInput("enrichment_title", "Titre de la figure :", placeholder = "Saisir un titre..."),
                checkboxInput("enrichment_toolbox", "Activer la barre d'outils", value = TRUE)
              )
            )
            
          ))
        ),
      #####=============================Documentation===========================
      tabItem(
        tabName = "documentation",
        h2("Coming soon")
      ),
      
      #####============================À PROPOS=======================================
      tabItem(
        tabName = "about",
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