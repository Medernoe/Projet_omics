#===========================Presentation========================================
#Authors : Noé Méderlet, Mehdi Tachekort, Mathieu Cartier, Valentin Fourdreane
#contact : noe.mederlet@univ-rouen.fr, mehdi.tachekort@univ-rouen.fr, mathieu.cartier@univ-rouen.fr, valentin.fourdraine@univ-rouen.fr
#github : https://github.com/Medernoe/Projet_omics
#organism : Master Bims M2, université de rouen 
#project : Création d'une application interactive dédiée à l'analyse de données transcriptomiques,
#développée dans le cadre d'un projet universitaire du Master 2 de Bioinformatique de l'Université de Rouen.
#==============================================================================#

#=========Chargement du global.R===============================================
source("global.R")

#=========================INTERFACE UTILISATEUR=================================

dashboardPage(
  
  ####============================Header==========================================
  
  # Le titre est en CSS
  dashboardHeader(title = ""),
  
  ####===========================Sidebar==========================================
  dashboardSidebar(
    
    ####==========================Accueil========================================
    sidebarMenu(
      id="menu_home",
      menuItem("Accueil", tabName = "home", icon = icon("house"))
    ),
    
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
      id = "menu_navigation",
      menuItem("DEG", tabName = "deg", icon = icon("magnifying-glass-chart")),
      
      # Sous-menu Enrichissement (4 sous-pages)
      menuItem("Enrichissement", icon = icon("chart-column"), startExpanded = FALSE,
               menuSubItem("GO ORA",       tabName = "go_ora",       icon = icon("chart-bar")),
               menuSubItem("GO GSEA",      tabName = "go_gsea",      icon = icon("chart-line")),
               menuSubItem("Pathway ORA",  tabName = "pathway_ora",  icon = icon("diagram-project")),
               menuSubItem("Pathway GSEA", tabName = "pathway_gsea", icon = icon("share-nodes"))
      ),
      
      menuItem("Documentation", tabName = "documentation", icon = icon("book-open")),
      menuItem("À propos", tabName = "about", icon = icon("question"))
    )
  ),
  
  #===============================Body============================================
  dashboardBody(
    
    ####==========================Initialisation waiter==========================
    use_waiter(),
    waiter_show_on_load(html = hamster_loader_server, color = "#009688"),
    
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
            
            ######=====================Informations=========================================
            p(strong("Réalisée par :"), "Noé Mederlet"),
            p(strong("GitHub :"), tags$a(href = "https://github.com/Medernoe/Projet_omics", "Lien vers le dépôt", target = "_blank")),
            p(strong("Contact :"), tags$a(href = "noe.mederlet@univ-rouen.fr", "noe.mederlet@univ-rouen.fr")),
            br(),
            
            p("Pour plus d'informations sur l'utilisation de l'application, consultez la section ", strong("« Documentation »"), "."),
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
              status = "primary",
              solidHeader = TRUE,
              
              conditionalPanel(
                condition = "output.show_volcano_error",
                div(
                  style = "text-align: center; padding: 20px;",
                  p("Veuillez charger un fichier CSV au format attendu pour visualiser son Volcano Plot",
                    style = "font-size: 16px; color: #666; margin-bottom: 20px;"),
                  imageOutput("volcano_error_img", height = "450px", inline = TRUE)
                )
              ),
              
              conditionalPanel(
                condition = "!output.show_volcano_error",
                plotlyOutput("volcano_plot", height = "500px")
              ),
              
              downloadButton("downloadVolcano", "Télécharger")
            )
          ),
          
          ########-----Colonne droite : Paramètres seuils-----
          column(
            width = 4,
            
            box(
              title = "Titre",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              textInput("deg_title", "Entrer un titre pour la figure :")
            ),
            
            box(
              title = "ToolBox",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              checkboxInput("deg_toolbox", "Activer la barre d'outils", value = TRUE)
            ),
            
            box(
              title = "Seuils",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              sliderInput("fc_threshold", "Log2 FC :", min = 0, max = 5, value = 1, step = 0.5, width = "100%"),
              checkboxInput("show_vline", "Activer ligne verticale (logFC)", value = TRUE),
              sliderInput("pvalue_threshold", "P-value :", min = 0, max = 1, value = 0.05, width = "100%"),
              checkboxInput("show_hline", "Activer ligne horizontale (p-value)", value = TRUE)
            )
          )
        ),
        
        ######=========================Tableau de données===============================
        fluidRow(
          box(
            title = "Tableau des données",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            
            conditionalPanel(
              condition = "output.show_table_error",
              div(
                style = "text-align: center; padding: 50px;",
                h3(textOutput("table_error_text"), style = "color: #666;")
              )
            ),
            
            conditionalPanel(
              condition = "!output.show_table_error",
              DTOutput("deg_table")
            )
          )
        )
      ),
      
      #####==========================GO ORA===========================================
      tabItem(
        tabName = "go_ora",
        
        fluidRow(
          column(
            width = 6,
            box(
              title = "Ontologie GO",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              checkboxGroupInput(
                inputId = "go_ora_ontology",
                label = NULL,
                choices = c("Processus Biologique (BP)" = "BP",
                            "Composant Cellulaire (CC)" = "CC",
                            "Fonction Moléculaire (MF)" = "MF"),
                selected = "BP",
                inline = FALSE
              )
            )
          ),
          
          column(
            width = 6,
            box(
              title = "Lancer l'analyse",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              div(
                style = "text-align: center; padding: 10px;",
                actionButton(
                  inputId = "run_go_ora",
                  label = "Lancer GO ORA",
                  icon = icon("play"),
                  class = "btn-teal btn-lg",
                  style = "width: 100%;"
                ),
                br(), br(),
                div(
                  style = "font-size: 12px; color: #666; font-style: italic;",
                  textOutput("go_ora_status")
                )
              )
            )
          )
        ),
        
        fluidRow(
          column(
            width = 8,
            
            box(
              title = "Enrichissement Plot (GO ORA)",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              plotlyOutput("go_ora_plot", height = "500px"),
              downloadButton("downloadGoOra", "Télécharger")
            ),
            
            box(
              title = "Affichage",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              textInput("go_ora_title", "Titre de la figure :", placeholder = "Saisir un titre..."),
              checkboxInput("go_ora_toolbox", "Activer la barre d'outils", value = TRUE)
            )
          ),
          
          column(
            width = 4,
            
            box(
              title = "Paramètres du graphique",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              
              selectInput(
                inputId = "go_ora_selected_plot",
                label = "Type de visualisation :",
                choices = c("Dotplot", "Barplot", "Cnetplot", "Emapplot",
                            "Goplot", "Upsetplot", "Heatplot", "Manhattan"),
                selected = "Dotplot"
              ),
              
              selectInput(
                inputId = "go_ora_displayed_ontology",
                label = "Ontologie à afficher :",
                choices = c("Processus Biologique (BP)" = "BP",
                            "Composant Cellulaire (CC)" = "CC",
                            "Fonction Moléculaire (MF)" = "MF"),
                selected = "BP"
              ),
              
              sliderInput(
                inputId = "go_ora_top_n_terms",
                label = "Nombre de termes GO à afficher :",
                min = 5, max = 50, value = 10, step = 5, width = "100%"
              )
            ),
            
            box(
              title = "Direction des gènes",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              radioButtons(
                inputId = "go_ora_direction",
                label = NULL,
                choices = c("Sur-exprimés" = "up",
                            "Sous-exprimés" = "down",
                            "Les deux" = "both"),
                selected = "both"
              )
            )
          )
        )
      ),
      
      #####==========================GO GSEA==========================================
      tabItem(
        tabName = "go_gsea",
        
        fluidRow(
          column(
            width = 6,
            box(
              title = "Ontologie GO",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              checkboxGroupInput(
                inputId = "go_gsea_ontology",
                label = NULL,
                choices = c("Processus Biologique (BP)" = "BP",
                            "Composant Cellulaire (CC)" = "CC",
                            "Fonction Moléculaire (MF)" = "MF"),
                selected = "BP",
                inline = FALSE
              )
            )
          ),
          
          column(
            width = 6,
            box(
              title = "Lancer l'analyse",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              div(
                style = "text-align: center; padding: 10px;",
                actionButton(
                  inputId = "run_go_gsea",
                  label = "Lancer GO GSEA",
                  icon = icon("play"),
                  class = "btn-teal btn-lg",
                  style = "width: 100%;"
                ),
                br(), br(),
                div(
                  style = "font-size: 12px; color: #666; font-style: italic;",
                  textOutput("go_gsea_status")
                )
              )
            )
          )
        ),
        
        fluidRow(
          column(
            width = 8,
            
            box(
              title = "Enrichissement Plot (GO GSEA)",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              plotlyOutput("go_gsea_plot", height = "500px"),
              downloadButton("downloadGoGsea", "Télécharger")
            ),
            
            box(
              title = "Affichage",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              textInput("go_gsea_title", "Titre de la figure :", placeholder = "Saisir un titre..."),
              checkboxInput("go_gsea_toolbox", "Activer la barre d'outils", value = TRUE)
            )
          ),
          
          column(
            width = 4,
            
            box(
              title = "Paramètres du graphique",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              
              selectInput(
                inputId = "go_gsea_selected_plot",
                label = "Type de visualisation :",
                choices = c("Dotplot", "Cnetplot", "Emapplot", "Upsetplot",
                            "Heatplot","GSEArank", "Manhattan"),
                selected = "Dotplot"
              ),
              
              selectInput(
                inputId = "go_gsea_displayed_ontology",
                label = "Ontologie à afficher :",
                choices = c("Processus Biologique (BP)" = "BP",
                            "Composant Cellulaire (CC)" = "CC",
                            "Fonction Moléculaire (MF)" = "MF"),
                selected = "BP"
              ),
              
              sliderInput(
                inputId = "go_gsea_top_n_terms",
                label = "Nombre de termes GO à afficher :",
                min = 5, max = 50, value = 10, step = 5, width = "100%"
              )
            ),
            
            box(
              title = "Direction des pathways (NES)",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              radioButtons(
                inputId = "go_gsea_direction",
                label = NULL,
                choices = c("Activés (NES > 0)" = "up",
                            "Réprimés (NES < 0)" = "down",
                            "Les deux" = "both"),
                selected = "both"
              )
            )
          )
        )
      ),
      
      #####==========================PATHWAY ORA======================================
      tabItem(
        tabName = "pathway_ora",
        
        fluidRow(
          column(
            width = 6,
            box(
              title = "Bases de données",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              checkboxGroupInput(
                inputId = "pathway_ora_databases",
                label = NULL,
                choices = c("KEGG" = "KEGG",
                            "Reactome" = "Reactome"),
                selected = "KEGG",
                inline = FALSE
              )
            )
          ),
          
          column(
            width = 6,
            box(
              title = "Lancer l'analyse",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              div(
                style = "text-align: center; padding: 10px;",
                actionButton(
                  inputId = "run_pathway_ora",
                  label = "Lancer Pathway ORA",
                  icon = icon("play"),
                  class = "btn-teal btn-lg",
                  style = "width: 100%;"
                ),
                br(), br(),
                div(
                  style = "font-size: 12px; color: #666; font-style: italic;",
                  textOutput("pathway_ora_status")
                )
              )
            )
          )
        ),
        
        fluidRow(
          column(
            width = 8,
            
            box(
              title = "Enrichissement Plot (Pathway ORA)",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              plotlyOutput("pathway_ora_plot", height = "500px"),
              downloadButton("downloadPathwayOra", "Télécharger")
            ),
            
            box(
              title = "Affichage",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              textInput("pathway_ora_title", "Titre de la figure :", placeholder = "Saisir un titre..."),
              checkboxInput("pathway_ora_toolbox", "Activer la barre d'outils", value = TRUE)
            )
          ),
          
          column(
            width = 4,
            
            box(
              title = "Paramètres du graphique",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              
              selectInput(
                inputId = "pathway_ora_selected_plot",
                label = "Type de visualisation :",
                choices = c("Dotplot", "Barplot", "Cnetplot", "Emapplot",
                            "Upsetplot", "Heatplot"),
                selected = "Dotplot"
              ),
              
              selectInput(
                inputId = "pathway_ora_displayed_db",
                label = "Base de données à afficher :",
                choices = c("KEGG" = "KEGG",
                            "Reactome" = "Reactome"),
                selected = "KEGG"
              ),
              
              sliderInput(
                inputId = "pathway_ora_top_n_terms",
                label = "Nombre de pathways à afficher :",
                min = 5, max = 50, value = 10, step = 5, width = "100%"
              )
            ),
            
            box(
              title = "Direction des gènes",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              radioButtons(
                inputId = "pathway_ora_direction",
                label = NULL,
                choices = c("Sur-exprimés" = "up",
                            "Sous-exprimés" = "down",
                            "Les deux" = "both"),
                selected = "both"
              )
            )
          )
        )
      ),
      
      #####==========================PATHWAY GSEA=====================================
      tabItem(
        tabName = "pathway_gsea",
        
        fluidRow(
          column(
            width = 6,
            box(
              title = "Bases de données",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              checkboxGroupInput(
                inputId = "pathway_gsea_databases",
                label = NULL,
                choices = c("KEGG" = "KEGG",
                            "Reactome" = "Reactome"),
                selected = "KEGG",
                inline = FALSE
              )
            )
          ),
          
          column(
            width = 6,
            box(
              title = "Lancer l'analyse",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              div(
                style = "text-align: center; padding: 10px;",
                actionButton(
                  inputId = "run_pathway_gsea",
                  label = "Lancer Pathway GSEA",
                  icon = icon("play"),
                  class = "btn-teal btn-lg",
                  style = "width: 100%;"
                ),
                br(), br(),
                div(
                  style = "font-size: 12px; color: #666; font-style: italic;",
                  textOutput("pathway_gsea_status")
                )
              )
            )
          )
        ),
        
        fluidRow(
          column(
            width = 8,
            
            box(
              title = "Enrichissement Plot (Pathway GSEA)",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              plotlyOutput("pathway_gsea_plot", height = "500px"),
              downloadButton("downloadPathwayGsea", "Télécharger")
            ),
            
            box(
              title = "Affichage",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              textInput("pathway_gsea_title", "Titre de la figure :", placeholder = "Saisir un titre..."),
              checkboxInput("pathway_gsea_toolbox", "Activer la barre d'outils", value = TRUE)
            )
          ),
          
          column(
            width = 4,
            
            box(
              title = "Paramètres du graphique",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              
              selectInput(
                inputId = "pathway_gsea_selected_plot",
                label = "Type de visualisation :",
                choices = c("Dotplot", "Cnetplot", "Emapplot", "Upsetplot",
                            "Heatplot", "GSEArank"),
                selected = "Dotplot"
              ),
              
              selectInput(
                inputId = "pathway_gsea_displayed_db",
                label = "Base de données à afficher :",
                choices = c("KEGG" = "KEGG",
                            "Reactome" = "Reactome"),
                selected = "KEGG"
              ),
              
              sliderInput(
                inputId = "pathway_gsea_top_n_terms",
                label = "Nombre de pathways à afficher :",
                min = 5, max = 50, value = 10, step = 5, width = "100%"
              )
            ),
            
            box(
              title = "Direction des pathways (NES)",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              radioButtons(
                inputId = "pathway_gsea_direction",
                label = NULL,
                choices = c("Activés (NES > 0)" = "up",
                            "Réprimés (NES < 0)" = "down",
                            "Les deux" = "both"),
                selected = "both"
              )
            )
          )
        )
      ),
      
      #####=============================Documentation=================================
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
            status = "primary",
            solidHeader = TRUE,
            includeHTML("www/about.html")
          )
        )
      )
    )
  )
)