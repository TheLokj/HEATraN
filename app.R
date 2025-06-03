# Developed by DAHER Rayan, NAIT EL DJOUDI Lamia & LESAGE Louison (@thelokj).
# rayan.daher@univ-rouen.fr
# lamia.nait-el-djoudi@univ-rouen.fr
# louison.lesage@univ-rouen.fr
# Students at Rouen Normandy University
# Master of Bioinformatics, class M2.2 BIMS 2026 
# Last updated : 03/06/2025
# HEATraN version 1.0.0

library(shiny)
library(shinydashboard) 
library(shinyjs) 
library(shinycssloaders) 
library(shinyalert) 
library(fresh) 
library(ggplot2)
library(ggtext)
library(readxl) 
library(data.table)
library(DT)
library(ini)
library(fs)
library(knitr)
library(enrichplot)
library(shinydashboardPlus) 
library(GOSemSim)
library(topGO)
library(GO.db)
library(plyr)
library(clusterProfiler)
library(ReactomePA)
library(pathview)
library(DOSE)
library(enrichplot)
library(dplyr)

shinyAppDir("./bin/")

