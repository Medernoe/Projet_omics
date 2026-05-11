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