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
  
  # Logo + texte "DEGO" dans le titre du header.
  # Si l'ensemble est tronqué, ajouter `titleWidth = 280` dans dashboardHeader().
  dashboardHeader(
    title = tags$a(
      href = "#",
      style = "display: flex; align-items: center; gap: 8px; text-decoration: none;",
      tags$img(
        src   = "logo-modified.png",
        alt   = "Logo DEGO",
        style = "height: 45px; margin-top: 2px; margin-bottom: 2px;"
      ),
      tags$span(
        "DEGO",
        style = "font-weight: 700; font-size: 20px; color: white; letter-spacing: 1px;"
      )
    )
  ),
  
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
      selected = "Mus musculus"
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
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
      
      # Élément audio invisible : musique d'attente jouée pendant les analyses.
      # `preload="auto"` permet au navigateur de précharger le fichier dès le
      # chargement de la page, pour qu'il soit prêt au moment du clic.
      tags$audio(
        id = "musique_attente",
        src = "musique_attente.mp3",
        preload = "auto",
        style = "display: none;"
      ),
      
      # JS handler : permet à server$sendCustomMessage("disableCheckboxes", ...)
      # de désactiver/réactiver visuellement des cases d'un checkboxGroupInput.
      # Verrouille les ontologies/bases déjà calculées sans les décocher.
      tags$script(HTML("
        Shiny.addCustomMessageHandler('disableCheckboxes', function(msg) {
          var group = document.getElementById(msg.group_id);
          if (!group) return;
          // On commence par tout réactiver dans le groupe
          group.querySelectorAll('input[type=\"checkbox\"]').forEach(function(cb) {
            cb.disabled = false;
            var lbl = cb.closest('label');
            if (lbl) {
              lbl.style.opacity = '';
              lbl.style.cursor = '';
              lbl.title = '';
            }
          });
          // Puis on désactive uniquement les valeurs reçues
          (msg.values || []).forEach(function(val) {
            var cb = group.querySelector('input[type=\"checkbox\"][value=\"' + val + '\"]');
            if (cb) {
              cb.disabled = true;
              var lbl = cb.closest('label');
              if (lbl) {
                lbl.style.opacity = '0.5';
                lbl.style.cursor = 'not-allowed';
                lbl.title = 'Déjà calculé — modifiez le fichier, l\\'espèce ou les seuils pour relancer';
              }
            }
          });
        });
        
        // Handlers play/pause pour la musique d'attente.
        // Le clic sur un bouton 'Lancer' compte comme interaction utilisateur,
        // donc le play() ne sera pas bloqué par la politique d'autoplay.
        // On garde une référence au timer pour pouvoir l'annuler si stopMusic
        // est appelé avant la fin des 20s (cas d'erreur côté serveur par exemple).
        var musicTimer = null;
        Shiny.addCustomMessageHandler('playMusic', function(msg) {
          var audio = document.getElementById('musique_attente');
          if (!audio) return;
          
          // Si la musique tourne déjà (clic répété), on ne la relance pas.
          if (!audio.paused && audio.currentTime > 0) return;
          
          audio.currentTime = 0;
          audio.volume = (msg && msg.volume != null) ? msg.volume : 0.5;
          var p = audio.play();
          if (p !== undefined) {
            p.catch(function(err) {
              console.warn('Lecture audio bloquée :', err);
            });
          }
          
          // Arrêt automatique après 20 secondes avec fade-out doux.
          if (musicTimer) clearTimeout(musicTimer);
          musicTimer = setTimeout(function() {
            var step = audio.volume / 12;
            var fade = setInterval(function() {
              if (audio.volume - step > 0) {
                audio.volume = audio.volume - step;
              } else {
                audio.pause();
                audio.currentTime = 0;
                audio.volume = 0.5;
                clearInterval(fade);
              }
            }, 50);
          }, 20000);
        });
        Shiny.addCustomMessageHandler('stopMusic', function() {
          var audio = document.getElementById('musique_attente');
          if (!audio) return;
          // Fade out doux sur 600ms pour éviter une coupure brutale
          var step = audio.volume / 12;
          var fade = setInterval(function() {
            if (audio.volume - step > 0) {
              audio.volume = audio.volume - step;
            } else {
              audio.pause();
              audio.currentTime = 0;
              audio.volume = 0.5;
              clearInterval(fade);
            }
          }, 50);
        });
      "))
    ),
    
    ####=====================Contenu des onglets====================================
    tabItems(
      
      #####===========================Accueil=======================================
      tabItem(
        tabName = "home",
        
        # 1. CSS Personnalisé
        tags$head(
          tags$style(HTML("
                  .logo-container { text-align: center; margin-bottom: 20px; }
                  .logo-container img { height: 80px; margin: 0 20px; }
                  
                  .bubble-container { display: flex; justify-content: space-around; flex-wrap: wrap; margin-bottom: 30px; }
                  
                  .bubble-btn {
                      width: 250px;
                      height: 250px;
                      border-radius: 50%;
                      border: 4px solid rgba(255, 241, 107, 0.8)!important;
                      background-color: white;
                      background-size: cover;
                      background-position: center;
                      box-shadow: 0 4px 8px rgba(0,0,0,0.1);
                      transition: transform 0.3s, box-shadow 0.3s;
                      display: flex;
                      flex-direction: column;
                      align-items: center;
                      justify-content: center;
                      text-align: center;
                      padding: 20px;
                      cursor: pointer;
                  }
                  .bubble-btn:hover {
                    transform: scale(1.05); /* Zoom au survol */
                    box-shadow: 0 8px 16px rgba(0,0,0,0.2);
                  }
                  .bubble-title {
                    font-weight: bold; color: #333; background: rgba(255,255,255,0.9);
                    padding: 8px; border-radius: 5px; margin-top: auto;
                  }
                  .explication-box {
                    background-color: #f9f9f9; border-left: 5px solid rgba(255, 241, 107, 0.8); padding: 20px;
                    font-size: 16px; border-radius: 5px; box-shadow: 0 2px 5px rgba(0,0,0,0.1);
                  }
              "))
        ),
        
        # 2. Section des Logos (Fac et Master)
        div(class = "logo-container",
            tags$img(src = "logouni.png", alt = "Logo Université de Rouen"),
            tags$img(src = "logomaster.png", alt = "Logo Master BIMS")
        ),
        
        # 3. Logo principal de l'application
        div(style = "text-align: center; margin-bottom: 25px;",
            tags$img(src = "logo-modified.png", alt = "Logo de l'Application", style = "max-height: 250px; max-width: 100%; object-fit: contain;")
        ),
        
        # 4. Description introductive traduite
        div(style = "text-align: center; max-width: 900px; margin: 0 auto 50px auto; padding: 25px; background-color: #ffffff; border-radius: 12px; box-shadow: 0 4px 15px rgba(0,0,0,0.05); border-top: 4px solid #009688;",
            tags$h3(style = "color: #009688; font-weight: bold; margin-bottom: 20px;", 
                    "Bienvenue sur Differential expression Gene|GO"),
            tags$p(style = "font-size: 16px; color: #444; line-height: 1.6;",
                   "Votre solution complète pour l'analyse de données génomiques globales, l'enrichissement de termes de la Gene Ontology (GO) et l'exploration des voies biologiques (pathways)."),
            tags$p(style = "font-size: 16px; color: #444; line-height: 1.6;",
                   "Conçue pour les chercheurs et les biologistes, notre application facilite le processus complexe d'analyse des données transcriptomiques. En identifiant les termes GO significatifs et en révélant les voies enrichies, elle vous aide à propulser vos découvertes scientifiques et à approfondir votre compréhension des fonctions géniques, des processus biologiques et des interactions moléculaires.")
        ),
        
        # 5. Section des 3 Bulles cliquables
        div(class = "bubble-container",
            actionButton("btn_deg", label = div(class="bubble-title", "Inspection & Analyse Différentielle"),
                         class = "bubble-btn", style = "background-image:linear-gradient(rgba(0,150,136,0.55), rgba(255, 241, 107, 0.8)),url('img_deg.png');"),
            actionButton("btn_go", label = div(class="bubble-title", "Enrichissement GO"),
                         class = "bubble-btn", style = "background-image:linear-gradient(rgba(0,150,136,0.55), rgba(255, 241, 107, 0.8)),url('img_go.png');"),
            actionButton("btn_pathway", label = div(class="bubble-title", "Enrichissement KEGG/Reactome"),
                         class = "bubble-btn", style = "background-image:linear-gradient(rgba(0,150,136,0.55), rgba(255, 241, 107, 0.8)),url('img_pathway.png');")
        ),
        
        # 6. Zone d'affichage dynamique de l'explication
        uiOutput("explication_accueil")
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
              uiOutput("ui_go_ora_plot"), # UI OUTPUT ICI
              downloadButton("downloadGoOra", "Télécharger")
            )
          ),
          
          column(
            width = 4,
            box(
              title = "Affichage",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              textInput("go_ora_title", "Titre de la figure :", placeholder = "Saisir un titre..."),
              checkboxInput("go_ora_toolbox", "Activer la barre d'outils", value = TRUE)
            ),
            
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
              title = "Niveaux d'expressions",
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
        ),
        
        ######=====================Tableau GO termes (ORA)===========================
        fluidRow(
          box(
            title = "Tableau des termes GO enrichis (ORA)",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            DTOutput("go_ora_table")
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
              uiOutput("ui_go_gsea_plot"), # UI OUTPUT ICI
              downloadButton("downloadGoGsea", "Télécharger")
            )
          ),
          
          column(
            width = 4,
            box(
              title = "Affichage",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              textInput("go_gsea_title", "Titre de la figure :", placeholder = "Saisir un titre..."),
              checkboxInput("go_gsea_toolbox", "Activer la barre d'outils", value = TRUE)
            ),
            
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
            )
            # Pas de box "Niveaux d'expressions" pour GSEA (volontairement retiré)
          )
        ),
        
        ######=====================Tableau GO termes (GSEA)==========================
        fluidRow(
          box(
            title = "Tableau des termes GO enrichis (GSEA)",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            DTOutput("go_gsea_table")
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
              uiOutput("ui_pathway_ora_plot"), # UI OUTPUT ICI
              downloadButton("downloadPathwayOra", "Télécharger")
            )
          ),
          
          column(
            width = 4,
            box(
              title = "Affichage",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              textInput("pathway_ora_title", "Titre de la figure :", placeholder = "Saisir un titre..."),
              checkboxInput("pathway_ora_toolbox", "Activer la barre d'outils", value = TRUE)
            ),
            
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
              title = "Niveaux d'expressions",
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
        ),
        
        ######=====================Tableau pathways (ORA)============================
        fluidRow(
          box(
            title = "Tableau des pathways enrichis (ORA)",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            DTOutput("pathway_ora_table")
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
              uiOutput("ui_pathway_gsea_plot"), # UI OUTPUT ICI
              downloadButton("downloadPathwayGsea", "Télécharger")
            )
          ),
          
          column(
            width = 4,
            box(
              title = "Affichage",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              textInput("pathway_gsea_title", "Titre de la figure :", placeholder = "Saisir un titre..."),
              checkboxInput("pathway_gsea_toolbox", "Activer la barre d'outils", value = TRUE)
            ),
            
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
            )
          )
        ),
        
        ######=====================Tableau pathways (GSEA)===========================
        fluidRow(
          box(
            title = "Tableau des pathways enrichis (GSEA)",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            DTOutput("pathway_gsea_table")
          )
        )
      ),
      
      #####=============================Documentation=================================
      tabItem(
        tabName = "documentation",
        fluidRow(
          column(
            width = 12,
            box(
              title = "Documentation Officielle & Guide d'Utilisation",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              
              div(style = "padding: 15px;",
                  
                  # Introduction
                  tags$p(style = "font-size: 16px; font-weight: bold; color: #2c3e50;",
                         "Bienvenue dans la documentation officielle de DEGO, une application interactive dédiée à l’analyse et à la visualisation de données transcriptomiques. Ce guide vous accompagnera étape par étape dans la préparation de vos données, l’exploration de l’expression différentielle et les analyses d’enrichissement fonctionnel."
                  ),
                  tags$hr(),
                  
                  # 1. Démarrage et Importation
                  tags$h3(icon("upload"), " 1. Démarrage et Importation des Données"),
                  tags$p("Pour garantir le bon fonctionnement de l’application, vos données doivent respecter un formatage précis avant d’être importées."),
                  
                  tags$h4(tags$strong("Préparation de votre fichier CSV")),
                  tags$p("Exportez vos données transcriptomiques (par exemple, issues de DESeq2 ou EdgeR) au format CSV avec un séparateur point-virgule (;). Votre fichier doit obligatoirement contenir les colonnes nommées exactement de la manière suivante :"),
                  tags$ul(
                    tags$li(tags$code("GeneName"), " : Le symbole officiel du gène (ex: TP53, BRCA1). ", tags$i("Attention : le format SYMBOL est requis pour que l’enrichissement fonctionne correctement.")),
                    tags$li(tags$code("log2FC"), " : Le Log2 Fold Change, indiquant le niveau de surexpression ou de sous-expression du gène."),
                    tags$li(tags$code("pval"), " : La p-value (ou p-value ajustée), représentant la significativité statistique de l’expression différentielle.")
                  ),
                  
                  tags$h4(tags$strong("Choix de l’organisme")),
                  tags$p("Dans le panneau latéral de gauche, utilisez le menu déroulant pour sélectionner l’organisme modèle correspondant à vos données. Les organismes actuellement supportés sont :"),
                  tags$ul(
                    tags$li(tags$i("Homo sapiens"), " (Humain)"),
                    tags$li(tags$i("Mus musculus"), " (Souris)"),
                    tags$li(tags$i("Drosophila melanogaster"), " (Mouche drosophile)")
                  ),
                  tags$p("Une fois l’organisme sélectionné, importez votre fichier CSV via le bouton ", tags$strong("Choisir fichier DEG"), ". L’application vérifiera automatiquement l’intégrité de vos données."),
                  tags$hr(),
                  
                  # 2. Inspection DEG
                  tags$h3(icon("magnifying-glass-chart"), " 2. Inspection des Données (Onglet DEG)"),
                  tags$p("Cette section vous permet d’explorer visuellement et tabulairement vos gènes différentiellement exprimés."),
                  
                  tags$h4(tags$strong("Volcano Plot Interactif")),
                  tags$p("Le Volcano Plot offre une vue globale de l’expression différentielle. Chaque point représente un gène."),
                  tags$ul(
                    tags$li(tags$strong("Seuils dynamiques : "), "Utilisez les curseurs situés à droite pour ajuster en temps réel le seuil de Log2 FC et le seuil de P-value. Les gènes significatifs seront automatiquement mis en évidence (rouge pour surexprimés, bleu pour sous-exprimés, gris pour non significatifs)."),
                    tags$li(tags$strong("Lignes de repère : "), "Vous pouvez activer ou désactiver l’affichage des lignes horizontales et verticales marquant vos seuils de coupure."),
                    tags$li(tags$strong("Interactivité : "), "Au survol d’un point, le nom du gène et ses valeurs exactes s’affichent."),
                    tags$li(tags$strong("Barre d’outils (ToolBox) : "), "Une boîte à outils activable permet de zoomer, sélectionner une zone ou réinitialiser la vue.")
                  ),
                  
                  tags$h4(tags$strong("Tableau de Données et Mise en Évidence")),
                  tags$p("Sous le graphique se trouve un tableau interactif listant vos gènes."),
                  tags$p(tags$strong(icon("lightbulb"), " Astuce : "), "Si vous cliquez sur une ligne spécifique du tableau, le gène correspondant sera mis en surbrillance (point violet et étiquette textuelle) directement sur le Volcano Plot."),
                  tags$hr(),
                  
                  # 3. Enrichissement ORA
                  tags$h3(icon("chart-pie"), " 3. Enrichissement ORA (Over-Representation Analysis)"),
                  tags$p("L’analyse ORA permet de déterminer si certains termes de la Gene Ontology (GO) sont statistiquement surreprésentés parmi vos gènes significativement différentiels (filtrés selon les seuils choisis dans l’onglet DEG)."),
                  
                  tags$h4(tags$strong("Paramétrage")),
                  tags$ul(
                    tags$li(tags$strong("Ontologie GO : "), "Sélectionnez une ou plusieurs catégories à analyser : Processus Biologique (BP), Composant Cellulaire (CC), Fonction Moléculaire (MF)."),
                    tags$li(tags$strong("Lancement : "), "Cliquez sur le bouton vert ", tags$strong("Lancer l’enrichissement ORA"), ". ", tags$i("Note : Le calcul peut prendre quelques secondes selon le nombre de gènes."))
                  ),
                  
                  tags$h4(tags$strong("Visualisation")),
                  tags$p("Une fois le calcul terminé, vous pouvez explorer les résultats via un menu déroulant proposant plusieurs représentations graphiques :"),
                  tags$ul(
                    tags$li(tags$strong("Dotplot : "), "Affiche les termes enrichis sous forme de bulles (taille = nombre de gènes, couleur = significativité)."),
                    tags$li(tags$strong("Barplot / Goplot : "), "Diagrammes classiques des termes les plus représentés."),
                    tags$li(tags$strong("Cnetplot : "), "Réseau reliant les termes GO aux gènes associés, permettant de voir les gènes partagés entre différentes voies."),
                    tags$li(tags$strong("Emapplot : "), "Carte de similarité regroupant les termes GO fonctionnellement proches."),
                    tags$li(tags$strong("Upsetplot : "), "Visualise les intersections complexes de gènes entre différents termes GO."),
                    tags$li(tags$strong("Heatplot : "), "Heatmap reliant les gènes aux termes GO.")
                  ),
                  tags$p("Vous pouvez ajuster le nombre de termes à afficher via le curseur dédié et télécharger la figure générée en haute résolution."),
                  tags$hr(),
                  
                  # 4. Enrichissement GSEA
                  tags$h3(icon("chart-line"), " 4. Enrichissement GSEA (Gene Set Enrichment Analysis)"),
                  tags$p("Contrairement à l’ORA qui se base sur un seuil strict, l’analyse GSEA prend en compte l’ensemble de vos gènes, classés par ordre décroissant selon leur Log2FC. Cela permet de détecter des variations d’expression subtiles mais coordonnées au sein d’une même voie biologique."),
                  
                  tags$h4(tags$strong("Paramétrage et Lancement")),
                  tags$p("Le fonctionnement est similaire à l’onglet ORA. Sélectionnez vos ontologies (BP, CC, MF) et lancez l’analyse. ", tags$i("L’algorithme GSEA étant plus lourd, le temps de calcul peut s’étendre à quelques minutes.")),
                  
                  tags$h4(tags$strong("Visualisation Spécifique")),
                  tags$p("En plus des graphiques communs avec l’ORA (Dotplot, Cnetplot, Emapplot, etc.), le module GSEA propose des visualisations expertes :"),
                  tags$ul(
                    tags$li(tags$strong("Ridgeplot : "), "Montre la distribution des Log2FC pour les gènes appartenant aux termes GO enrichis."),
                    tags$li(tags$strong("GSEAplot2 : "), "Affiche le score d’enrichissement classique (“running score”) pour plusieurs voies simultanément."),
                    tags$li(tags$strong("GSEArank : "), "Permet d’isoler le graphique de score (enrichment plot) pour une voie spécifique.")
                  ),
                  tags$hr(),
                  
                  # 5. Structure Technique
                  tags$h3(icon("cogs"), " 5. Structure Technique de l’Application"),
                  tags$p("L’application DEGO a été développée en R dans le cadre du Master 2 Bioinformatique de l’Université de Rouen. Elle repose sur une architecture modulaire pour faciliter sa maintenance et sa scalabilité."),
                  
                  tags$h4(tags$strong("Fichiers Sources")),
                  tags$ul(
                    tags$li(tags$code("ui.R"), " / ", tags$code("server.R"), " : Gèrent respectivement l’interface utilisateur et la logique réactive du serveur."),
                    tags$li(tags$code("global.R"), " : Initialise l’environnement, charge les bibliothèques et source les dépendances."),
                    tags$li(tags$code("fonctions.R"), " : Contient toute la logique métier abstraite (fonctions de nettoyage, calculs statistiques, algorithmes GSEA/ORA et création des objets graphiques).")
                  ),
                  
                  tags$h4(tags$strong("Dépendances et Librairies Utilisées")),
                  tags$p("Le bon fonctionnement de DEGO repose sur les packages suivants :"),
                  
                  tags$h5(icon("r-project"), " Packages CRAN (Interface, Traitement & Visualisation) :"),
                  tags$p(
                    tags$code("shiny"), " ", tags$code("shinydashboard"), " ", tags$code("waiter"), " ", 
                    tags$code("ggplot2"), " ", tags$code("DT"), " ", tags$code("plotly"), " ", 
                    tags$code("shinyalert"), " ", tags$code("ggarchery"), " ", tags$code("qqman"), " ", 
                    tags$code("dplyr")
                  ),
                  
                  tags$h5(icon("dna"), " Packages Bioconductor (Analyse Transcriptomique & Ontologies) :"),
                  tags$p(
                    tags$code("clusterProfiler"), " ", tags$code("org.Hs.eg.db"), " ", 
                    tags$code("org.Mm.eg.db"), " ", tags$code("org.Dm.eg.db"), " ", 
                    tags$code("enrichplot"), " ", tags$code("ReactomePA"), " ", tags$code("DOSE")
                  )
              )
            )
          )
        )
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
            
            div(class = "about_description", style = "padding: 15px;",
                
                # --- Logo centré ---
                div(style = "text-align: center; margin-bottom: 30px;",
                    tags$img(src = "logo-modified.png", height = "150px")
                ),
                
                tags$h3("Le projet DEGO"),
                tags$p("DEGO (Differentially Expressed Genes & Enrichment Gene Ontology) est une application interactive que nous avons développée pour faciliter l'analyse de données transcriptomiques. Elle permet d'explorer facilement les gènes différentiellement exprimés (DEG) et de réaliser des analyses d'enrichissement (GO, KEGG, Reactome) via des approches ORA et GSEA."),
                
                tags$h3("Code source et documentation"),
                tags$p("Le projet est entièrement open-source. Vous pouvez retrouver le code, la documentation détaillée et des jeux de données d'exemple directement sur notre dépôt GitHub :"),
                tags$a(href = "https://github.com/Medernoe/Projet_omics", target = "_blank", "https://github.com/Medernoe/Projet_omics"),
                
                tags$h3("L'équipe"),
                tags$p("Cette application a été pensée et développée par notre groupe de 4 étudiants en Master 2 Bioinformatique et Modélisation (BIMS) à l'Université de Rouen Normandie :"),
                tags$ul(
                  tags$li(tags$strong("Noé Méderlet : "), tags$a(href="mailto:noe.mederlet@univ-rouen.fr", "noe.mederlet@univ-rouen.fr")),
                  tags$li(tags$strong("Mehdi Tachekort : "), tags$a(href="mailto:mehdi.tachekort@univ-rouen.fr", "mehdi.tachekort@univ-rouen.fr")),
                  tags$li(tags$strong("Mathieu Cartier : "), tags$a(href="mailto:mathieu.cartier@univ-rouen.fr", "mathieu.cartier@univ-rouen.fr")),
                  tags$li(tags$strong("Valentin Fourdraine : "), tags$a(href="mailto:valentin.fourdraine@univ-rouen.fr", "valentin.fourdraine@univ-rouen.fr"))
                ),
                
                # --- Section Demande d'ajout d'espèces ---
                tags$h3("Besoin d'une autre espèce ?"),
                tags$p("L'application intègre actuellement les annotations pour l'Humain, la Souris et la Drosophile. Si vous travaillez sur un organisme différent et que vous souhaitez utiliser notre outil, n'hésitez pas à nous envoyer un e-mail ! Nous serons ravis d'ajouter votre espèce d'intérêt lors d'une prochaine mise à jour."),
                
                tags$h3("Détails techniques"),
                tags$p("L'application est développée en R et utilise les bibliothèques ", tags$code("shiny"), " et ", tags$code("shinydashboard"), " pour son interface utilisateur. Pour la fluidité de l'expérience, nous utilisons également des modules comme ", tags$code("shinyalert"), " et ", tags$code("waiter"), " pour gérer les temps de chargement et les notifications."),
                tags$p("Concernant le traitement des données et l'analyse bioinformatique, l'application s'appuie sur des packages de référence de la communauté Bioconductor. L'analyse d'enrichissement (ORA et GSEA) est notamment réalisée via ", tags$code("clusterProfiler"), " et ", tags$code("ReactomePA"), ". Les conversions d'identifiants entre espèces et l'accès aux ontologies sont gérés de manière dynamique par les bases de données d'annotation spécifiques (", tags$code("org.Hs.eg.db"), ", ", tags$code("org.Mm.eg.db"), ", et ", tags$code("org.Dm.eg.db"), ")."),
                tags$p("Enfin, les visualisations graphiques reposent fortement sur ", tags$code("enrichplot"), " pour la construction des graphes biologiques (Cnetplot, Emapplot, etc.), combiné à ", tags$code("ggplot2"), " pour les représentations plus poussées, comme nos versions personnalisées du Volcano plot et du Manhattan plot. L'interactivité des figures et des tableaux d'exploration est quant à elle assurée par ", tags$code("plotly"), " et ", tags$code("DT"), "."),
                
                # --- Citations et Références ---
                tags$h3("Références"),
                tags$p("Nous nous sommes appuyés sur plusieurs publications scientifiques et méthodes de référence pour le développement de cette application :"),
                tags$ul(
                  tags$li(
                    tags$strong("Over-representation analysis (ORA) : "), "Nguyen, T. M., et al. (2019). ", 
                    tags$em("Over-representation analysis: a comprehensive review."), 
                    " BMC Bioinformatics. ", 
                    tags$a(href="https://doi.org/10.1186/s12859-019-2710-2", target="_blank", "[DOI: 10.1186/s12859-019-2710-2]")
                  ),
                  tags$li(
                    tags$strong("Méthode GSEA : "), "Subramanian A, et al. (2005). ", 
                    tags$em("Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles."), 
                    " PNAS. ", 
                    tags$a(href="https://doi.org/10.1073/pnas.0506580102", target="_blank", "[DOI: 10.1073/pnas.0506580102]")
                  ),
                  tags$li(
                    tags$strong("clusterProfiler 4.0 : "), "Wu T, et al. (2021). ", 
                    tags$em("clusterProfiler 4.0: A universal enrichment tool for interpreting omics data."), 
                    " The Innovation. ", 
                    tags$a(href="https://doi.org/10.1016/j.xinn.2021.100141", target="_blank", "[DOI: 10.1016/j.xinn.2021.100141]")
                  ),
                  tags$li(
                    tags$strong("ReactomePA : "), "Yu G & He QY (2016). ", 
                    tags$em("ReactomePA: an R/Bioconductor package for reactome pathway analysis and visualization."), 
                    " Molecular BioSystems. ", 
                    tags$a(href="https://doi.org/10.1039/C5MB00663E", target="_blank", "[DOI: 10.1039/C5MB00663E]")
                  )
                ),
                
                tags$h3("Remerciements"),
                tags$p("Ce projet a été réalisé dans le cadre de notre formation. Un grand merci à nos encadrantes, Hélène Dauchel et Solène Pety, pour leur aide et leurs conseils tout au long du développement !")
            )
          )
        )
      )
    )
  )
)