FROM bioconductor/bioconductor_docker:RELEASE_3_18

# Installation de Shiny Server
RUN apt-get update && apt-get install -y \
    gdebi-core \
    && wget https://download3.rstudio.org/ubuntu-18.04/x86_64/shiny-server-1.5.20.1002-amd64.deb \
    && gdebi -n shiny-server-1.5.20.1002-amd64.deb \
    && rm shiny-server-1.5.20.1002-amd64.deb

# Installation des dépendances système supplémentaires
RUN apt-get update && apt-get install -y \
    libssl-dev \
    libcurl4-gnutls-dev \
    libxml2-dev \
    libgit2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Installation des packages CRAN
RUN R -e "install.packages(c( \
    'shiny', \
    'shinydashboard', \
    'shinyjs', \
    'shinyalert', \
    'fresh', \
    'ggplot2', \
    'readxl', \
    'data.table', \
    'DT', \
    'shinydisconnect'), \
    repos='https://cloud.r-project.org/')"

# Installation explicite des packages Bioconductor
RUN R -e 'BiocManager::install(c("ReactomePA", "clusterProfiler", "pathview"), update=FALSE, ask=FALSE)'

RUN R -e 'BiocManager::install(c("org.At.eg.db", "org.Bt.eg.db", "org.Cf.eg.db", "org.Gg.eg.db", "org.EcK12.eg.db", "org.Dm.eg.db", "org.Hs.eg.db", "org.Mm.eg.db", "org.Ss.eg.db", "org.Rn.eg.db", "org.Ce.eg.db", "org.Xl.eg.db", "org.Sc.sgd.db", "org.Dr.eg.db"))'

# Vérification de l'installation des packages Bioconductor
RUN R -e 'library(clusterProfiler); library(ReactomePA); library(pathview); cat("Tous les packages Bioconductor sont installés correctement\n")'

# Configuration du répertoire de travail
WORKDIR /srv/shiny-server/HEATraN

# Création de l'arborescence de l'application
RUN mkdir -p bin data www

# Copie de l'application respectant la structure du projet
COPY ./app.R ./
COPY ./bin/ ./bin/
COPY ./bin/fun/ ./bin/fun/
COPY ./data/ ./data/
COPY ./www/ ./www/

# Définition des permissions
RUN chmod -R 755 .

# Configuration pour afficher les erreurs détaillées
RUN echo "sanitize_errors false;" >> /etc/shiny-server/shiny-server.conf

# Exposition du port pour Shiny
EXPOSE 3838

# Commande de démarrage
CMD ["/usr/bin/shiny-server"]
