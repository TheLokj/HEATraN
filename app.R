# Developed by DAHER Rayan, NAIT EL DJOUDI Lamia & LESAGE Louison (@thelokj).
# rayan.daher@univ-rouen.fr
# lamia.nait-el-djoudi@univ-rouen.fr
# louison.lesage@univ-rouen.fr
# Students at Rouen Normandy University
# Master of Bioinformatics, class M2.2 BIMS 2026 
# Last updated : 04/06/2025
# HEATraN version 1.0.0

library(shiny)
library(shinydashboard) 
library(shinydashboardPlus) 
library(shinyjs) 
library(shinycssloaders) 
library(shinyalert) 
library(plyr)
library(dplyr)
library(data.table)
library(DT)
library(readxl) 
library(ini)
library(fresh) 
library(ggplot2)
library(ggtext)
library(topGO)
library(GO.db)
library(clusterProfiler)
library(ReactomePA)
library(pathview)
library(enrichplot)
library(knitr)

shinyAppDir("./bin/")

