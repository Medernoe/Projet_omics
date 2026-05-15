FROM rocker/shiny:4.3.1

# Dépendances système
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# On copie le lockfile et le dossier renv en premier
COPY renv.lock /app/renv.lock
COPY renv/ /app/renv/

# Creation du .Rprofile dans l'image 
RUN echo 'source("renv/activate.R")' > .Rprofile

# Installation de renv et restauration
RUN R -e "install.packages('renv', repos='https://cloud.r-project.org')"
RUN R -e "renv::restore()"

# Copie du reste des fichiers (UI, Server, Global, etc.)
COPY . /app

EXPOSE 3838
CMD ["R", "-e", "shiny::runApp('/app', host = '0.0.0.0', port = 3838)"]