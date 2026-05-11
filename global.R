# global.R

# ── Liste des packages requis ──────────────────────────────────────────────────
packages <- c(
  # Shiny & UI
  "shiny",
  "shinydashboard",
  "waiter",
  "shinyalert",
  
  # Visualisation
  "ggplot2",
  "plotly",
  "ggarchery",
  "enrichplot",
  
  # Tableaux
  "DT",
  
  # Enrichissement
  "clusterProfiler",
  
  # Annotation (Bioconductor)
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  "org.Dm.eg.db"
)

# ── Fonction d'installation intelligente ──────────────────────────────────────
install_if_missing <- function(pkgs) {
  
  # Sépare CRAN et Bioconductor
  bioc_pkgs <- c("clusterProfiler", "enrichplot", "org.Hs.eg.db", "org.Mm.eg.db", "org.Dm.eg.db")
  cran_pkgs  <- setdiff(pkgs, bioc_pkgs)
  
  # Packages non installés
  missing_cran <- cran_pkgs[!cran_pkgs %in% installed.packages()[, "Package"]]
  missing_bioc <- bioc_pkgs[!bioc_pkgs %in% installed.packages()[, "Package"]]
  
  # Installation CRAN
  if (length(missing_cran) > 0) {
    message("📦 Installation CRAN : ", paste(missing_cran, collapse = ", "))
    install.packages(missing_cran, dependencies = TRUE)
  }
  
  # Installation Bioconductor
  if (length(missing_bioc) > 0) {
    message("🧬 Installation Bioconductor : ", paste(missing_bioc, collapse = ", "))
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(missing_bioc, ask = FALSE)
  }
  
  message("✅ Tous les packages sont disponibles.")
}

install_if_missing(packages)

# ── Chargement ─────────────────────────────────────────────────────────────────
invisible(lapply(packages, library, character.only = TRUE))

source("fonctions.R")