# 1. Image de base
FROM rocker/shiny:4.4.2

# 2. Dépendances système Linux
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    cmake \
    libglpk-dev \
    libuv1-dev \
    curl \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libfontconfig1-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# 3. Création du dossier de l'app
WORKDIR /app

# 4. Installation de BiocManager (sans forcer le repo)
RUN R -e "install.packages('BiocManager')"

# 5. Installation groupée de TOUS les packages via BiocManager
# BiocManager va piocher dans le bon snapshot CRAN de Rocker et dans Bioconductor 3.20
RUN R -e "options(warn=2); BiocManager::install(c( \
    'shinydashboard', 'waiter', 'ggplot2', 'DT', 'plotly', \
    'shinyalert', 'ggarchery', 'qqman', 'dplyr', \
    'clusterProfiler', 'org.Hs.eg.db', 'org.Mm.eg.db', 'org.Dm.eg.db', \
    'enrichplot', 'ReactomePA', 'DOSE', 'ggtree' \
    ), ask = FALSE, update = FALSE)"

# 6. Copie des fichiers de votre application
COPY . /app

# 7. Exposition du port et lancement
EXPOSE 3838
CMD ["R", "-e", "shiny::runApp('/app', host = '0.0.0.0', port = 3838)"]