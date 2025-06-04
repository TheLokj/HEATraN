FROM bioconductor/bioconductor_docker:RELEASE_3_18

# Installing Shiny Server
RUN apt-get update && apt-get install -y \
    gdebi-core \
    && wget https://download3.rstudio.org/ubuntu-18.04/x86_64/shiny-server-1.5.20.1002-amd64.deb \
    && gdebi -n shiny-server-1.5.20.1002-amd64.deb \
    && rm shiny-server-1.5.20.1002-amd64.deb

# Installing additional system dependencies
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
    tini \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Installing CRAN packages
RUN R -e "install.packages(c( \
    'shiny', \
    'shinydashboard', \
    'shinydashboardPlus', \
    'shinyjs', \
    'shinyalert', \
    'shinycssloaders', \
    'fresh', \
    'ggplot2', \
    'ggtext', \
    'ggridges', \
    'ggupset', \
    'readxl', \
    'data.table', \
    'DT', \
    'ini', \
    'plyr', \
    'knitr' \
  ), repos='https://cloud.r-project.org/')"


# Installation of Bioconductor packages
RUN R -e 'BiocManager::install( \
  c("ReactomePA", "clusterProfiler", "pathview", \
    "enrichplot", "topGO", "GO.db", \
    "org.At.tair.db", "org.Bt.eg.db", "org.Cf.eg.db", "org.Gg.eg.db", \
    "org.EcK12.eg.db", "org.Dm.eg.db", "org.Hs.eg.db", "org.Mm.eg.db", \
    "org.Ss.eg.db", "org.Rn.eg.db", "org.Ce.eg.db", "org.Xl.eg.db", \
    "org.Sc.sgd.db", "org.Dr.eg.db" \
  ), \
  update=FALSE, ask=FALSE \
)'

# Checking the installation of Bioconductor packages
RUN R -e 'library(clusterProfiler); library(ReactomePA); library(pathview); cat("Tous les packages Bioconductor sont installés correctement\n")'

# Configuring the working directory
WORKDIR /srv/shiny-server/HEATraN

# Creating the application tree
RUN mkdir -p bin data www

# Copy of the application respecting the project structure
COPY conf.ini ./
COPY app.R ./
COPY bin/ ./bin/
COPY bin/fun/ ./bin/fun/
COPY data/ ./data/
COPY www/ ./www/

# Setting permissions
RUN chmod -R a+rX ./

# Configuration to display detailed errors (debug)
RUN echo "sanitize_errors false;" >> /etc/shiny-server/shiny-server.conf

EXPOSE 3838

CMD ["Rscript", "./app.R"]
