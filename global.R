#===========================CHARGEMENT DES LIBRAIRIES===========================
# Packages CRAN
library(shiny)
library(shinydashboard)
library(waiter)
library(ggplot2)
library(DT)
library(plotly)
library(shinyalert)
library(ggarchery)
library(qqman)
library(dplyr)

# Packages Bioconductor
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Dm.eg.db)
library(enrichplot)
library(ReactomePA)
library(DOSE)

#===========================CHARGEMENT DES FONCTIONS============================
source("fonctions.R")

####============================WAITER HELPERS===============================
hamster_loader_server <- HTML('
<div aria-label="Hamster loader" role="img" class="wheel-and-hamster">
  <div class="wheel"></div>
  <div class="hamster">
    <div class="hamster__body">
      <div class="hamster__head">
        <div class="hamster__ear"></div>
        <div class="hamster__eye"></div>
        <div class="hamster__nose"></div>
      </div>
      <div class="hamster__limb hamster__limb--fr"></div>
      <div class="hamster__limb hamster__limb--fl"></div>
      <div class="hamster__limb hamster__limb--br"></div>
      <div class="hamster__limb hamster__limb--bl"></div>
      <div class="hamster__tail"></div>
    </div>
  </div>
  <div class="spoke"></div>
</div>
')

# Fabrique un Waiter qui s'affiche pendant un calcul lourd
# - id : ID de l'output sur lequel afficher le loader
#        Si NULL, le loader est plein écran
# - message : message texte affiché sous le hamster
make_waiter <- function(id = NULL, message = "Calcul en cours...") {
  Waiter$new(
    id    = id,
    html  = tagList(
      hamster_loader_server,
      tags$h4(message, style = "color: #009688; margin-top: 20px; text-align: center;")
    ),
    color = "rgba(255, 255, 255, 0.85)" 
  )
}