# Image par défaut pré construite de R + shiny
FROM rocker/shiny:4.4.2

# Dépendances système Linux
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

# Création du dossier de l'app
WORKDIR /app

# Installation de BiocManager (sans forcer le repo)
RUN R -e "install.packages('BiocManager')"

# Installation groupée de tous les packages via BiocManager
#option pour afficher les erreurs et stopper la construction de l'image
RUN R -e "options(warn=2); BiocManager::install(c( \
    'shinydashboard', 'waiter', 'ggplot2', 'DT', 'plotly', \
    'shinyalert', 'ggarchery', 'qqman', 'dplyr', \
    'clusterProfiler', 'org.Hs.eg.db', 'org.Mm.eg.db', 'org.Dm.eg.db', \
    'enrichplot', 'ReactomePA', 'DOSE', 'ggtree' \
    ), ask = FALSE, update = FALSE)"

# Copie des fichiers de l'application
COPY . /app

# Exposition du port et lancement
EXPOSE 3838
CMD ["R", "-e", "shiny::runApp('/app', host = '0.0.0.0', port = 3838)"]