# HEATraN

## About

HEATraN, litteraly ***H**yper-**E**xpression **A**nalysis **T**ool **ra**mpantly developed in **N**ormandy* is a bioinformatics analysis tool dedicated to transcriptomic analysis. 
It was developed as part of a student project in the Bioinformatics Master of Rouen Normandy University.

<img src='./www/logo.png' width='256' alt='HEATraN logo' style='display:block;margin-left: auto;margin-right: auto;'>

## How to use it

You must first get the repository. You can either download from GitHub by clicking on `Code>Download ZIP`, or enter directly into your terminal: 

`git clone https://github.com/TheLokj/HEATraN.git`

Then, HEATraN can be launched using several ways. 

### Online

*WIP*

### Locally

The two below methods allow you to access the tool by connecting to the local address http://localhost:3838/HEATraN/.

#### Docker

The safest option is to build the HEATraN image using the dockerfile provided. 
This ensures that all the necessary dependencies and libraries are installed, in the correct versions and compatible with each other. 

`docker build -t heatran-app .`

You can then launch it by typing:

`docker run -p 3838:3838 heatran-app R -e "shiny::runApp('/srv/shiny-server/HEATraN/app.R', host='0.0.0.0', port=3838)"`

#### R & Rstudio

The simplest option is to launch directly the app with R, from a terminal or directly within Rstudio.
This requires you to have installed the necessary R packages and associated software libraries first (see below for more details). 

Then, in Rstudio, simply go to the root of the project, open the `app.R` script and click on ‘Run App’ in the top right-hand corner.

From a terminal, you can enter : 

`R -e "shiny:runApp('app.r', port=3838')"`

#### Requirements for manual installation 

##### Bioconductor Packages

| **Library** | **Version** |
| :-- | :-- |
| R | 4.3.3 |
| shiny | 1.10.0 |
| shinyalert | 3.1.0 |
| shinycssloaders | 1.1.0 |
| shinydashboard | 0.7.3 |
| shinydisconnect | 0.1.1 |
| shinyjs | 2.1.0 |
| fresh | 0.2.1 |
| ggplot2 | 3.5.2 |
| ggtext | 0.1.2 |
| readxl | 1.4.5 |
| data.table | 1.17.0 |
| DT | 0.33 |

##### Bioconductor Packages

| **Library** | **Version** |
| :-- | :-- |
| BiocManager | 1.30.23 |
| BiocVersion | 3.18.1 |
| ReactomePA | 1.46.0 |
| clusterProfiler | 4.10.1 |
| pathview | 1.42.0 |
| org.At.tair.db | 3.18.0 |
| org.Bt.eg.db | 3.18.0 |
| org.Cf.eg.db | 3.18.0 |
| org.Gg.eg.db | 3.18.0 |
| org.EcK12.eg.db | 3.18.0 |
| org.Dm.eg.db | 3.18.0 |
| org.Hs.eg.db | 3.18.0 |
| org.Mm.eg.db | 3.18.0 |
| org.Ss.eg.db | 3.18.0 |
| org.Rn.eg.db | 3.18.0 |
| org.Ce.eg.db | 3.18.0 |
| org.Xl.eg.db | 3.18.0 |
| org.Sc.sgd.db | 3.18.0 |
| org.Dr.eg.db | 3.18.0 |

If one of these packages does not install, make sure that all the required software libraries are installed:

`sudo apt-get install libssl-dev libcurl4-gnutls-dev libxml2-dev libgit2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev`

## Repository structure

```
HEATraN/
├── app.R                # Main script
├── dockerfile           # Docker file to build Docker
├── config.ini           # Configuration file
├── bin/                 # Directory containing main scripts
│   ├── server.R         # Shiny Server script
│   ├── ui.R             # Shiny UI script
│   └── fun/             # Directory containing subscripts
│         └── pathway.R  # R script managing pathways
├── data/                # Data directory
│   └── example.tsv      # Differentially expressed results example file
└── www/                 # Additional graphics resources directory
    ├── logo.png         # HEATraN logo 
    ├── template.Rmd     # Template for report export
    ├── doc.HTML         # Documentation n°1 
    ├── stat.HTML        # Documentation n°2
    └── style.css        # HEATraN custom CSS stylesheet
```

## Citation

To cite this project, please refers to:

`Daher R., Naid El Djoudi L., Lesage L., Dauchel H. (2025). HEATraN. [https://github.com/TheLokj/HEATraN].`

Please also remember to cite the authors of the packages used by HEATraN and listed above.