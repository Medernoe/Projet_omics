#===========================INSTALLATION AUTOMATIQUE============================
# Ce bloc vérifie que tous les packages nécessaires sont installés.
# - Les packages CRAN sont installés via install.packages()
# - Les packages Bioconductor sont installés via BiocManager::install()
# Au premier lancement, l'installation peut prendre 15-30 minutes.

# ---- 1. Liste des packages CRAN ----
cran_packages <- c(
  "shiny",
  "shinydashboard",
  "waiter",
  "ggplot2",
  "DT",
  "plotly",
  "shinyalert",
  "ggarchery",
  "qqman",
  "dplyr"
)

# ---- 2. Liste des packages Bioconductor ----
bioc_packages <- c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  "org.Dm.eg.db",
  "enrichplot",
  "ReactomePA",
  "DOSE"
)

# ---- 3. Fonction d'installation et de chargement ----
install_and_load <- function(packages, installer = install.packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste0("📦 Installation de '", pkg, "' en cours..."))
      installer(pkg)
    }
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  }
}

# ---- 4. Bootstrap : on s'assure que BiocManager est dispo ----
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  message("📦 Installation de BiocManager (gestionnaire Bioconductor)...")
  install.packages("BiocManager")
}

# ---- 5. Installation + chargement ----
message("🔍 Vérification des packages CRAN...")
install_and_load(cran_packages, installer = install.packages)

message("🔍 Vérification des packages Bioconductor...")
install_and_load(bioc_packages, installer = BiocManager::install)

message("✅ Tous les packages sont chargés.")

#===========================FIN INSTALLATION====================================

#===========================CHARGEMENT DES FONCTIONS============================
source("fonctions.R")

####============================WAITER HELPERS===============================

# Le HTML du hamster loader (dupliqué depuis ui.R car le server n'y a pas accès)
# Si tu veux éviter la duplication, déplace cette variable dans global.R
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
# - id : ID de l'output sur lequel afficher le loader (ex: "go_ora_plot")
#        Si NULL, le loader est plein écran
# - message : message texte affiché sous le hamster
make_waiter <- function(id = NULL, message = "Calcul en cours...") {
  Waiter$new(
    id    = id,
    html  = tagList(
      hamster_loader_server,
      tags$h4(message, style = "color: #009688; margin-top: 20px; text-align: center;")
    ),
    color = "rgba(255, 255, 255, 0.85)"  # fond semi-transparent
  )
}