# Developed by DAHER Rayan, NAIT EL DJOUDI Lamia & LESAGE Louison (@thelokj).
# rayan.daher@univ-rouen.fr
# lamia.nait-el-djoudi@univ-rouen.fr
# louison.lesage@univ-rouen.fr
# Students at Rouen Normandy University
# Master of Bioinformatics, class M2.2 BIMS 2026 
# Last updated : 05/06/2025
# HEATraN version 1.0.0

pkgs <- c("shiny", "shinydashboard", "shinydashboardPlus", "shinyjs", "shinycssloaders", "shinyalert", "shinydashboardPlus",
          "plyr", "dplyr", "data.table", "DT", "readxl", "ini", "fresh",
          "ggplot2", "ggtext", "topGO", "GO.db", "grid", "Rgraphviz",
          "clusterProfiler", "ReactomePA", "pathview", "enrichplot", "knitr")

invisible(suppressWarnings(suppressPackageStartupMessages(lapply(pkgs, function(pkg) {
        require(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts  = FALSE)
      }))))

source("bin/ui.R", local=F)
source("bin/server.R", local=F)
source("bin/fun/pathway.R", local=F)

app <- shinyApp(
  ui, server,
  onStart = function() {
    cat("                                                                             \n")
    cat(" DAHER Rayan, NAIT EL DJOUDI Lamia & LESAGE Louison present...               \n")
    cat("                                                                             \n")
    cat(":::    ::: ::::::::::     ::: ::::::::::: :::::::::      :::     ::::    ::: \n")
    cat(":+:    :+: :+:          :+: :+:   :+:     :+:    :+:   :+: :+:   :+:+:   :+: \n")
    cat("+:+    +:+ +:+         +:+   +:+  +:+     +:+    +:+  +:+   +:+  :+:+:+  +:+ \n")
    cat("+#++:++#++ +#++:++#   +#++:++#++: +#+     +#++:++#:  +#++:++#++: +#+ +:+ +#+ \n")
    cat("+#+    +#+ +#+        +#+     +#+ +#+     +#+    +#+ +#+     +#+ +#+  +#+#+# \n")
    cat("#+#    #+# #+#        #+#     #+# #+#     #+#    #+# #+#     #+# #+#   #+#+# \n")
    cat("###    ### ########## ###     ### ###     ###    ### ###     ### ###    #### \n")
    cat("                                                                       v1.0.0\n")
    cat("                                                                             \n")
    cat("Connect yourself on below link to access to the app from your favorite browser!\n")
  },
  onStop(function() {
    if (config$FILE$clear_cache == "True") {
      if (file.exists("./out")) {
        if ((dir_size("./out") / (1024^2)) > config$FILE$max_cache_mb) {
          message("Cache exceeding the authorised limit: cleaning cacache...")
          files_to_delete <- list.files(path = "./out", full.names = T)
          unlink(files_to_delete)
        }
      }
    }
    message("Thanks for using HEATraN!")
  })
)

runApp(app, host='0.0.0.0', port=3838)

