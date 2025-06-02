# Developed by DAHER Rayan, NAIT EL DJOUDI Lamia & LESAGE Louison (@thelokj).
# rayan.daher@univ-rouen.fr
# lamia.nait-el-djoudi@univ-rouen.fr
# louison.lesage@univ-rouen.fr
# Students at Rouen Normandy University
# Master of Bioinformatics, class M2.2 BIMS 2026 
# Last updated : 13/05/2025
# HEATraN version 0.3.0

setwd("../")
source("./bin/fun/pathway.R")

config <- read.ini("./conf.ini")

app_palette = colorRampPalette(c("#ef940b", "#7e3535"))
options(enrichplot.colours = c("#ef940b", "#7e3535"))

emptyTable <- data.frame(Gene=NA, Log2FC=NA, p_value=NA)
emptyTableGo <- data.frame(GO=NA, Description=NA, p_value=NA, q_value=NA)
emptyTable2 <- data.frame(Pathway=NA, p_value=NA, q_value=NA)
brushInfo <- reactiveVal(NULL)
sigGO <- reactiveVal(NULL)
sigGOresult <- reactiveVal(NULL)
pathwayORAEnrichment <- reactiveVal(NULL)
pathwayGSEAEnrichment <- reactiveVal(NULL)
goOraResults <- reactiveVal(NULL)
goGseaResults <- reactiveVal(NULL)
selectionMode <- reactiveVal("None")
preprocessedData <- reactiveVal(emptyTable)
requiredNames <- c(config$DATA$gene_name, config$DATA$gene_id, config$DATA$basemean, config$DATA$log2FC,  config$DATA$pvalue, config$DATA$adjusted_pvalue)

emptyPlot <- ggplot() + 
  geom_richtext(aes(x = 0.5, y = 0.5, 
                    label = "<b style='font-family:Mulish,sans-serif;font-size:20pt;'>Data required for visual representation</b>"), 
                size = 6, color = "red3", fill = NA, label.color = NA) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = NA))

emptyGseaPlot <- ggplot() + 
  geom_richtext(aes(x = 0.5, y = 0.5, 
                    label = "<b style='font-family:Mulish,sans-serif;font-size:20pt;'>Please select at least one enriched pathway</b>"), 
                size = 6, color = "red3", fill = NA, label.color = NA) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = NA))

emptyPathwayPlot <- ggplot() + 
  geom_richtext(aes(x = 0.5, y = 0.5, 
                    label = "<b style='font-family:Mulish,sans-serif;font-size:20pt;'>Please select an enriched data</b>"), 
                size = 6, color = "red3", fill = NA, label.color = NA) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = NA))

tooSmallTreePlot <- ggplot() + 
  geom_richtext(aes(x = 0.5, y = 0.5, 
                    label = "<b style='font-family:Mulish,sans-serif;font-size:20pt;'>Tree plot requires at least 2 significant data</b>"), 
                size = 6, color = "red3", fill = NA, label.color = NA) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = NA))

function(input, output, session) {
  
  hideTab(inputId = "go_tabs", target = "Results (ORA)", session = session)
  hideTab(inputId = "go_tabs", target = "Results (GSEA)", session = session)
  
  hideTab(inputId = "pathway_tabs", target = "Results (ORA)", session = session)
  hideTab(inputId = "pathway_tabs", target = "Results (GSEA)", session = session)
  hideTab(inputId = "pathway_tabs", target = "View Pathway", session = session)
  
  # -----------------------------------------
  # Data processing
  # -----------------------------------------
  
  # Reactive function controlling the imported data
  importedData <- reactive ({
    message("Importing data")
    # While there is no input, return NULL
    if(is.null(input$tableInput)){
      return(NULL)}
    # After the importation, print the loaded table
    else{
      inputPath <- (input$tableInput)[1,4]
      format <- tail(strsplit(inputPath, ".", fixed=T)[[1]], 1)
      # Check the file format
      if(format=="csv" || format=="tsv"){
        table = fread(inputPath, sep="auto", h=T)
      } else if(format=="xlsx" | format=="xls"){
        # Inform the user for the excel files
        shinyalert("Excel file imported", "The first datasheet of this file will be imported.", confirmButtonCol = "#7e3535")
        table <- read_excel(inputPath, sheet=1)
      }
      # Print an alert if the user select another file
      else{
        shinyalert("Wrong format!", "Format accepted : .csv, .tsv, .xls and .xlsx", type = "error", confirmButtonCol = "#7e3535")
        return(NULL)
      }
      # Check if the table is correctly build (with 6 columns)
      if(dim(table)[2]>=6){
        return(table)
      } else {
        shinyalert("Incorrect table!", "The table must at least be built of six columns", type = "error", confirmButtonCol = "#7e3535")
        return(NULL)
      }
    }
  }) |> bindEvent(input$tableInput, ignoreNULL=F, ignoreInit=T)
  
  # Observe importation in order to preprocess data
  observe ({
    message("Preprocessing data")
    dataToPreprocess <- importedData()
    # If the required columns aren't found, 
    print(requiredNames)
    if (FALSE %in% c(requiredNames %in% colnames(dataToPreprocess))){
      shinyalert(html = TRUE, "Please select the variables", "Wrong column name", type = "info", confirmButtonCol = "#7e3535",
                 text = tagList(
                   selectInput("GeneNameCol", "Gene Name", colnames(dataToPreprocess)),
                   selectInput("GeneIDCol", "Gene ID", colnames(dataToPreprocess)),
                   selectInput("BaseMeanCol", "Base mean", colnames(dataToPreprocess)),
                   selectInput("Log2FCCol", "Log2(FoldChange)", colnames(dataToPreprocess)),
                   selectInput("pvalCol", "p-value", colnames(dataToPreprocess)),
                   selectInput("padjCol", "Adjusted p-value", colnames(dataToPreprocess)),
                   checkboxInput("saveColNames", "use these column names in the future"),
                   HTML(paste("Next time, use these names to import direcly the table : <i>", paste(requiredNames, collapse=", "), "</i><br>")),
                 ),
                 callbackR = function(finished) { 
                   if(finished) {
                     if (length(unique(c(input$GeneNameCol, input$GeneIDCol, input$BaseMeanCol,input$Log2FCCol,input$pvalCol, input$padjCol)))==6){
                       if (input$saveColNames==T) { 
                         config$DATA$gene_name = input$GeneNameCol
                         config$DATA$gene_id = input$GeneIDCol
                         config$DATA$basemean = input$BaseMeanCol
                         config$DATA$log2FC = input$Log2FCCol
                         config$DATA$pvalue = input$pvalCol
                         config$DATA$adjusted_pvalue = input$padjCol
                         write.ini(config, "./conf.ini")
                       }
                       dataToPreprocess <- importedData()
                       colnames(dataToPreprocess)[colnames(dataToPreprocess)==input$GeneNameCol] = "GeneName"
                       colnames(dataToPreprocess)[colnames(dataToPreprocess)==input$GeneIDCol] = "GeneID"
                       colnames(dataToPreprocess)[colnames(dataToPreprocess)==input$BaseMeanCol] = "baseMean"
                       colnames(dataToPreprocess)[colnames(dataToPreprocess)==input$Log2FCCol] = "Log2FC"
                       colnames(dataToPreprocess)[colnames(dataToPreprocess)==input$pvalCol] = "pval"
                       colnames(dataToPreprocess)[colnames(dataToPreprocess)==input$padjCol] = "padj"
                       dataToPreprocess$minuslog10 <- -log(dataToPreprocess$pval)
                       preprocessedData(dataToPreprocess)
                       selectionMode("Sliders")
                     } else {
                       shinyalert("Incorrect choice!", "Each column must be unique.", type = "error", confirmButtonCol = "#7e3535")
                     }
                   } else {return(FALSE)}})
    } else {
      colnames(dataToPreprocess)[colnames(dataToPreprocess)==config$DATA$gene_name] = "GeneName"
      colnames(dataToPreprocess)[colnames(dataToPreprocess)==config$DATA$gene_id] = "GeneID"
      colnames(dataToPreprocess)[colnames(dataToPreprocess)==config$DATA$base_mean] = "baseMean"
      colnames(dataToPreprocess)[colnames(dataToPreprocess)==config$DATA$log2FC] = "Log2FC"
      colnames(dataToPreprocess)[colnames(dataToPreprocess)==config$DATA$pvalue] = "pval"
      colnames(dataToPreprocess)[colnames(dataToPreprocess)==config$DATA$adjusted_pvalue] = "padj"
      dataToPreprocess$minuslog10 <- -log(dataToPreprocess$pval)
      preprocessedData(dataToPreprocess)
      selectionMode("Sliders")}}) |> bindEvent(importedData())
  
  # Reactive function containing the selected points
  processedData <- reactive({
    message("Processing data")
    if (!isTRUE(all.equal(emptyTable, preprocessedData()))){
      df <- preprocessedData()
      updateSliderInput(session,'Log2FC',max=ceiling(max(abs(df$Log2FC))))
      # Adapt the plot limit when user zoom in it
      if (Zoom()$Zoomed==T){
        df <- df[df$Log2FC>=Zoom()$coords[1]&df$Log2FC<=Zoom()$coords[2]&(-log(df$padj))>=Zoom()$coords[3]&(-log(df$padj))<=Zoom()$coords[4],]
      }
      # Update the selected points according to the selection mode
      if (selectionMode() == "Brush") {
        df$selected <- ifelse(df$GeneID%in%brushedPoints(df, brushInfo()())$GeneID, "TRUE", "FALSE")
      } else if (selectionMode() == "Sliders"){
        df$selected <- ifelse((df$Log2FC>input$Log2FC&df$padj<input$pval)|(df$Log2FC<(-input$Log2FC)&df$padj<input$pval), "TRUE", "FALSE")
      }
      updateNumericInput(session, "export_GeneNumber", value=nrow(df), min=0, max=nrow(df))
      return(df)
    } else {
      return(emptyTable)
    }
  })
  
  output$goGseaPlot2 <- renderPlot({
    if (!is.null(goGseaResults()) && length(input$goGSEA) > 0 && input$goGSEA[1] != "None") {
      data = goGseaResults()
      plot = gseaplot2(data, as.numeric(input$goGSEA), 
                       pvalue_table = TRUE,
                       color = c("#E495A5", "#86B875", "#7DB0DD"))
      return(plot)
    } else {
      return(emptyGseaPlot)
    }
  })
  
  
    # Reactive function building the plot
    plot <- reactive({
      message("Loading plot")
      plot = ggplot(processedData(), aes(x=Log2FC, y=minuslog10, col=selected)) +
        geom_point() +
        scale_color_manual(values=c("FALSE" = "#384246", "TRUE" = "#E69F00")) +
        theme(legend.position = "none") +
        xlab("log2(FoldChange)") +
        ylab("-log10(p-value ajusted)") +
        theme(text = element_text(size = 14))    
      return(plot)})
    
    pathwayGSEADotPlot <- reactive({
      if (!is.null(pathwayGSEAEnrichment())) {
        data = pathwayGSEAEnrichment()$enrichment
        result_df <- data@result
        result_df[result_df$p.adjust <= input$qval, ]
        data@result <- result_df
        plot = dotplot(data, showCategory=30) + ggtitle("Dotplot" )  
        return(plot)
      } else {
        return(emptyPlot)
      }}
      )
    
    pathwayGseaPlot <- reactive({
      if ((!is.null(input$pathwayGSEA)) & (length(input$pathwayGSEA)>0)) {
        data = pathwayGSEAEnrichment()$enrichment
        plot = gseaplot2(data, as.numeric(input$pathwayGSEA)) 
        return(plot)
      } else {
        return(emptyGseaPlot)
      }
     })
    
    pathwayGseaUpsetPlot <- reactive({
      if (!is.null(pathwayGSEAEnrichment())) {
        plot = upsetplot(pathwayGSEAEnrichment()$enrichment, col="#7e3535") 
        return(plot)
      } else {
        return(emptyPlot)
      }
    })
    
    pathwayORAUpsetPlot <- reactive({
      if (!is.null(pathwayORAEnrichment())) {
        plot = upsetplot(pathwayORAEnrichment()$enrichment) 
        return(plot)
      } else {
        return(emptyPlot)
      }
    })
    
    pathwayORADotPlot <- reactive({
      if (!is.null(pathwayORAEnrichment())) {
        data = pathwayORAEnrichment()$enrichment
        result_df <- data@result
        result_df[result_df$p.adjust <= input$qval, ]
        data@result <- result_df
        plot = dotplot(data, showCategory=30) + ggtitle("Dotplot" )  
        return(plot)
      } else {
        return(emptyPlot)
      }}
    )
    
    pathwayORATreePlot <- reactive({
      # 1) récupérer le résultat ORA (une liste contenant $enrichment)
      res <- pathwayORAEnrichment()
      if (is.null(res)) {
        return(emptyPlot)
      }
      
      # 2) extraire le data.frame des résultats et compter les feuilles
      enr    <- res$enrichment
      n_leaf <- nrow(enr)
      
      # 3) gérer les tout petits cas (0 ou 1 terme)
      if (n_leaf <= 1) {
        return(tooSmallTreePlot)
      }
      
      # 4) calculer les similarités entre termes
      dataR <- pairwise_termsim(
        setReadable(
          enr,
          orgs[orgs$organism == input$species, "db"],
          keyType = "ENTREZID"
        )
      )
      
      # 5) définir k = n_leaf - 1 pour respecter 1 <= k <= n_leaf-1
      k_dyn <- n_leaf - 1
      #    (optionnel : plafonner k à une valeur max, p.ex. min(k_dyn, 5))
      
      # 6) tracer l’arbre avec k sur-mesure
      treeplot(
        dataR,
        k = k_dyn
      )
    })
    
    pathwayORACnetPlot <- reactive({
      if (!is.null(pathwayORAEnrichment())) {
        data = pathwayORAEnrichment()$enrichment
        dataR = setReadable(data, orgs[orgs$organism==input$species, "db"], 'ENTREZID')
        plot = cnetplot(dataR)
        return(plot)
      } else {
        return(emptyPlot)
      }
    })
    
    pathwayORAEmapPlot <- reactive({
      if (!is.null(pathwayORAEnrichment())) {
        data = pathwayORAEnrichment()$enrichment
        dataR = setReadable(data, orgs[orgs$organism==input$species, "db"], 'ENTREZID')
        plot = emapplot(pairwise_termsim(dataR))
        return(plot)
      } else {
        return(emptyPlot)
      }
    })
    
    pathwayGSEAEmapPlot <- reactive({
      if (!is.null(pathwayGSEAEnrichment())) {
        data = pathwayGSEAEnrichment()$enrichment
        dataR = setReadable(data, orgs[orgs$organism==input$species, "db"], 'ENTREZID')
        plot = emapplot(pairwise_termsim(dataR))
        return(plot)
      } else {
        return(emptyPlot)
      }
    })
    
    pathwayGSEATreePlot <- reactive({
      
      # 1) on récupère la liste
      res <- pathwayGSEAEnrichment()
      if (is.null(res)) {
        return(emptyPlot)
      }
      
      # 2) on extrait le tableau et le vecteur foldChange
      enr    <- res$enrichment     # data.frame des pathways
      geneList <- res$geneList     # named numeric vector
      
      # 3) on compte le nombre de termes
      n_leaf <- nrow(enr)
      if (n_leaf <= 1) {
        return(tooSmallTreePlot)
      }
      
      # 4) préparation du pairwise_termsim
      dataR <- pairwise_termsim(
        setReadable(
          enr,
          orgs[orgs$organism == input$species, "db"],
          keyType = "ENTREZID"
        )
      )
      
      # 5) ajustement dynamique de k
      k_dyn <- n_leaf - 1
      # si vous voulez plafonner, par ex. à max 5 clusters :
      # k_dyn <- min(n_leaf - 1, 5)
      
      # 6) on trace
      treeplot(
        dataR,
        foldChange = geneList,
        k          = k_dyn
      )
    })
    
    
    pathwayGSEACnetPlot <- reactive({
      if (!is.null(pathwayGSEAEnrichment())) {
        data = pathwayGSEAEnrichment()$enrichment
        dataR = setReadable(data, orgs[orgs$organism==input$species, "db"], 'ENTREZID')
        plot = cnetplot(dataR, foldChange=data@geneList)
        return(plot)
      } else {
        return(emptyPlot)
      }
    })
    
    observeEvent(input$pathwayGSEA, {
      if (!is.null(pathwayGSEAEnrichment()) & (input$pathwayGSEA[1]=="None")) {
        data = pathwayGSEAEnrichment()$enrichment
        if (pathwayEnrichment()$parameters$DB=="KEGG"){
          updateSelectInput(session ,"pathwayGSEA", choices = setNames(1:length(data$ID), data$Description))
        }
        else if (pathwayGSEAEnrichment()$parameters$DB=="Reactome"){
          updateSelectInput(session ,"pathwayGSEA", choices = setNames(1:length(data$Description), data$Description))
        }}
    })
    
    observe({
      if (!is.null(goOraResults())) { 
        showTab(inputId = "go_tabs", target = "Results (ORA)", session = session)
      } else {
        hideTab(inputId = "go_tabs", target = "Results (ORA)", session = session)
      }
      if (!is.null(goGseaResults())) { 
        showTab(inputId = "go_tabs", target = "Results (GSEA)", session = session)
      } else {
        hideTab(inputId = "go_tabs", target = "Results (GSEA)", session = session)
      }
    })
    
    # -----------------------------------------
    # User event
    # -----------------------------------------
      
    # Reactive function controling the selection mode 
    # As shiny do not allow to disable downloadButton, disable it with shinyJS
    # Disable also other buttons with shinyJS to standardize rendering
    
    observeEvent(input$q_val, {
      config$STAT$q_val = input$q_val
      write.ini(config, "./conf.ini")
    })
    
    observeEvent(input$ajust_method, {
      config$STAT$adjust_method = input$ajust_method
      write.ini(config, "./conf.ini")
    })
    
    #--- Whole Data Inspection ---#
    
    observeEvent(selectionMode(), {
      message(paste("Selection mode:", selectionMode()))
      if (selectionMode() != "None") {
        enable('Download')
        enable('DownloadTable')
        enable('SelectAll')
        enable('ResetButton')
        removeClass("Download", "disabled-button")
        removeClass("DownloadTable", "disabled-button")
        removeClass("SelectAll", "disabled-button")
        removeClass("ResetButton", "disabled-button")
        if (selectionMode() == "Sliders"){
          enable('pval')
          enable('Log2FC')
          if (Zoom()$Zoomed==F){
            disable('ZoomButton')
            addClass("ZoomButton", "disabled-button")
          } else {
            updateActionButton(session, "ZoomButton", label = "Unzoom selection", icon=icon('zoom-out', lib='glyphicon'))
          }
        } else if (selectionMode() == "Brush"){
          disable('pval')
          disable('Log2FC')
          enable('ZoomButton')
          removeClass("ZoomButton", "disabled-button")
          updateActionButton(session, "ZoomButton", label = "Zoom in selection", icon=icon('zoom-in', lib='glyphicon'))
        }
      } else {
        disable('pval')
        disable('Log2FC')
        disable('SelectAll')
        disable('ZoomButton')
        disable('ResetButton')
        addClass("SelectAll", "disabled-button")
        addClass("ZoomButton", "disabled-button")
        addClass("ResetButton", "disabled-button")
        addClass("Download", "disabled-button")
        addClass("DownloadTable", "disabled-button")
      }})
  
  # Download button event
  output$Download <- downloadHandler(
    filename <- function() { paste(unlist(strsplit(input$tableInput[,1], ".", fixed=T))[1], "_HEATraNplot.pdf") },
    content <- function(file) {
      ggsave(file, plot(), units = "mm", height=210, width=297)
    }
  )
  # Export table event
  output$DownloadTable <- downloadHandler(
    filename <- function() { paste(unlist(strsplit(input$tableInput[,1], ".", fixed=T))[1], "_HEATraNtable.csv") },
    content <- function(file) {
      write.csv(processedData()[,-c("selected","minuslog10")], file)
    }
  )
  
  # Reset button event
  observe({
    updateSliderInput(session,'Log2FC',value = 1)
    updateSliderInput(session,'pval',value = 0.05)
    if (selectionMode() == "Brush"){
      selectionMode("Sliders")
      session$resetBrush("plot_brush")
      brushInfo(NULL)
      updateActionButton(session, "ZoomButton", label = "Unzoom selection", icon=icon('zoom-out', lib='glyphicon'))
      if (Zoom()$Zoomed==F){
        disable('ZoomButton')
        addClass("ZoomButton", "disabled-button")
      }
    }
  }) |> bindEvent(input$ResetButton)
  
  # Select All button event
  observe({
    updateSliderInput(session,'Log2FC',value = 0)
    updateSliderInput(session,'pval',value = 1)
    if (selectionMode() == "Brush"){
      selectionMode("Sliders")
      session$resetBrush("plot_brush")
      brushInfo(NULL)
      updateActionButton(session, "ZoomButton", label = "Unzoom selection", icon=icon('zoom-out', lib='glyphicon'))
      if (Zoom()$Zoomed==F){
        disable('ZoomButton')
        addClass("ZoomButton", "disabled-button")
      }
    }
  }) |> bindEvent(input$SelectAll)
  
  # Brush use
  observe({
    message("Brush use")
    selectionMode("Brush")
    brushInfo(reactiveVal(input$plot_brush))
  }) |> bindEvent(input$plot_brush)
  
  # Zoom button event
  Zoom = reactive({
    if (selectionMode() == "Brush"){
      message("Zoom")
      coords <- c(brushInfo()()$xmin, brushInfo()()$xmax, brushInfo()()$ymin, brushInfo()()$ymax)
      session$resetBrush("plot_brush")
      brushInfo(NULL)
      selectionMode("Sliders")
      return(list(Zoomed=T, coords=coords))}
    else {
      message("Unzoom")
      disable('ZoomButton')
      addClass("ZoomButton", "disabled-button")
      return(list(Zoomed=F, coords=NULL))
    }
  }) |> bindEvent(input$ZoomButton, ignoreNULL=F)
  
  observe({
    if (!is.null(pathwayORAEnrichment()) || !is.null(pathwayGSEAEnrichment())) {
      if (!is.null(pathwayORAEnrichment())) {
        showTab(inputId = "pathway_tabs", target = "Results (ORA)", session = session)
        showTab(inputId = "pathway_tabs", target = "View Pathway", session = session)
      } 
      if (!is.null(pathwayGSEAEnrichment())) {
        showTab(inputId = "pathway_tabs", target = "Results (GSEA)", session = session)
        showTab(inputId = "pathway_tabs", target = "View Pathway", session = session)
      }
    }
  })

 #--- GO Enrichment ---#
  
  observeEvent(input$go_analysisButton, {
    req(preprocessedData())
    
    withProgress(message = 'Running GO analysis...', {
      df <- preprocessedData()
      
      # List of genes for the universe
      original_gene_list <- df$Log2FC
      names(original_gene_list) <- df$GeneID
      gene_list <- sort(na.omit(original_gene_list), decreasing = T)
      
      # Gene selection by adjusted p-value threshold
      sig_genes_df <- df[df$padj < input$go_pval,]
      genes <- sig_genes_df$Log2FC
      names(genes) <- sig_genes_df$GeneID
      genes <- na.omit(genes)
      
      up_genes <- names(genes)[genes > log2(input$go_oraFC)]
      down_genes <- names(genes)[genes < -log2(input$go_oraFC)]
      
      # Selecting the GO ontology and organism
      ontology <- input$inputGO
      go_organism <- orgs[orgs$organism==input$species, "db"]
      universe <- names(gene_list)
      library(go_organism, character.only = T)
      
      if ("Gene Set Enrichment Analysis (GSEA)" %in% input$go_analysisMethodChoice) {
        setProgress(message("Computing GSEA enrichment..."), value=0.5)
        gene_list <- df$Log2FC
        names(gene_list) <- df$GeneID
        gene_list <- sort(na.omit(gene_list), decreasing = T)
        gsea_result <- tryCatch({gseGO(geneList = gene_list, OrgDb = get(go_organism), ont = ontology, keyType = "ENSEMBL", pvalueCutoff = input$go_pval, verbose = F,
                                       pAdjustMethod = config$STAT$adjust_method)}, error = function(e) {return(NULL)})
        if (is.null(gsea_result)){
          shinyalert("Database issue", text="Check the species and compatibility of the IDs (if they are not in the correct format, please use an online converter). If the problem is still unresolved, check your significance threshold, which may be too restrictive.", type = "error")
          return()
        }
        if (nrow(gsea_result)>0) {
          goGseaResults(gsea_result) 
        } else {
          goGseaResults(NULL)
        }
        
      } 
      
      if ("Over Representation Analysis (ORA)" %in% input$go_analysisMethodChoice) {
        setProgress(message("Computing ORA enrichment..."), value=0.75)
        analyze_up <- "Over expressed DEG" %in% input$go_oraChoice
        analyze_down <- "Under expressed DEG" %in% input$go_oraChoice
        
        if (!analyze_up && !analyze_down) {
          shinyalert("Error", "Select at least one type of gene expression (Over or Under expressed DEG).", type = "error")
          return(NULL)
        }
        
        # Déterminer quels gènes analyser selon la sélection
        if (analyze_up && analyze_down) {
          # Les deux : utiliser l'union des gènes up et down
          genes_to_analyze <- union(up_genes, down_genes)
          analysis_type <- "both"
        } else if (analyze_up) {
          # Seulement up-regulated
          genes_to_analyze <- up_genes
          analysis_type <- "up"
        } else {
          # Seulement down-regulated
          genes_to_analyze <- down_genes
          analysis_type <- "down"
        }
        
        # Effectuer une seule analyse au lieu de trois
        if (length(genes_to_analyze) > 0) {
          result_go <- tryCatch({enrichGO(
            gene = genes_to_analyze,
            universe = universe,
            OrgDb = go_organism,
            keyType = "ENSEMBL",
            readable = T,
            ont = ontology,
            pvalueCutoff = input$go_pval,
            pAdjustMethod = config$STAT$adjust_method,
            qvalueCutoff = as.numeric(config$STAT$q_val)
          )}, error = function(e) {return(NULL)})
          if (is.null(result_go)){
            shinyalert("Unable to enrich", text="Check the species and compatibility of the IDs (if they are not in the correct format, please use an online converter). If the problem is still unresolved, check your significance threshold, which may be too restrictive.", type = "error")
            return()
          }
          if (nrow(result_go)>0){
            goOraResults(list(result = result_go, type = analysis_type))
            
            allGenesBinary <- as.integer(names(gene_list) %in% genes_to_analyze)
            names(allGenesBinary) <- names(gene_list)
            
            GOdata <- new("topGOdata",
                          ontology = ontology,
                          allGenes = allGenesBinary,
                          geneSelectionFun = function(x)(x < 0.05),
                          annot = annFUN.org,
                          mapping = go_organism,
                          ID = "ensembl")
            
            sigGO(GOdata)
            sigGOresult(runTest(sigGO(), algorithm = "classic", statistic = "fisher"))
          }
          else {
            goOraResults(NULL)
          }
        } else {
          goOraResults(NULL)
        }

        }
      })
    })
  
  observeEvent(goGseaResults(), {
    if (!is.null(goGseaResults())) {
      data = goGseaResults()
      updateSelectInput(session, "goGSEA", 
                        choices = setNames(1:length(data@result$Description), 
                                           data@result$Description))
    }
  })
  
    #--- Pathway Enrichment ---#
    
    observe({
      if (is.na(preprocessedData()[1,1])){
        shinyalert("Bad input", "Analysis require data!", type = "error", confirmButtonCol = "#7e3535")
      }
      else if ("ORA" %in% input$analysisMethodChoice & is.null(input$oraChoice)){
        shinyalert("Incomplete selection!", "You should select an interest for ORA method", type = "error", confirmButtonCol = "#7e3535")
      }
      else {
        withProgress(message = "Pathway enrichment...", {
            message("Starting Pathway Enrichment...")
            df = as.data.frame(preprocessedData())
            if ("ORA" %in% input$analysisMethodChoice) {
              setProgress(message = "Computing ORA enrichment...")
              
              # 1. Run enrichment and store in a temp var
              res <- pathway(
                df,
                organism   = input$species,
                DB         = input$dbPathwaychoice,
                analysis   = "ORA",
                oraInterest= input$oraChoice,
                pval       = input$pvalPathway,
                threshold = input$pathway_oraFC
              )
              # 2. Test whether any enriched terms were returned
              if (is.null(res) || (is.null(res$enrichment)) ||  nrow(res$enrichment@result) == 0) {
                # 3a. No enriched terms: set to NULL
                shinyalert("No results for ORA", text="Check your significance and fold-change threshold, which may be too restrictive.", type = "error")
                pathwayORAEnrichment(NULL)
              } else {
                # 3b. Enriched terms found: pass result and show tab
                pathwayORAEnrichment(res)
                showTab(inputId = "pathway_tabs", target = "View Pathway", session = session)
              }
            }
            
            if ("GSEA" %in% input$analysisMethodChoice) {
              setProgress(message = "Computing GSEA enrichment...")
              
              res <- pathway(
                df,
                organism   = input$species,
                DB         = input$dbPathwaychoice,
                analysis   = "GSEA",
                oraInterest= input$oraChoice,
                pval       = input$pvalPathway
              )
              
              if (is.null(res) || (is.null(res$enrichment)) ||  nrow(res$enrichment@result) == 0) {
                shinyalert("No results for GSEA", text="Check your significance threshold, which may be too restrictive.", type = "error")
                pathwayGSEAEnrichment(NULL)
              } else {
                pathwayGSEAEnrichment(res)
                showTab(inputId = "pathway_tabs", target = "View Pathway", session = session)
              }
            }
        })
    }}) |> bindEvent(input$analysisPathwayButton)
    
    observeEvent(pathwayGSEAEnrichment(), {
      data = pathwayGSEAEnrichment()$enrichment
      if (pathwayGSEAEnrichment()$parameters$DB=="KEGG"){
        updateSelectInput(session ,"pathwayGSEA", choices = setNames(1:length(data$ID), data$Description))
        } else if (pathwayGSEAEnrichment()$parameters$DB=="Reactome"){
        updateSelectInput(session ,"pathwayGSEA", choices = setNames(1:length(data$Description), data$Description))
        }
    }) 
  
  # -----------------------------------------
  # UI output 
  # -----------------------------------------

  #--- Whole Data Inspection ---#
    
  output$InfoSelect <- renderUI({
    if (is.null(preprocessedData) || is.na(preprocessedData()[1,1])){
      HTML(paste("<b style='color:	#FF0000'>Data required for selection</b><br/><br/>", sep=""))
    }
    else {
      if (is.null(input$plot_brush)){
        HTML(paste("<b>Current selection</b><br/><i>p-value</i>: [0 ; ", input$pval,"]", "    <br/>    ", "<i>log2(FoldChange)</i>: ", "[-", input$Log2FC, " ; ", input$Log2FC,"]<br/><br/>", sep=""))
      }
      else {
        HTML(paste("<b>Current selection</b><br/><i>p-value</i>: [", format(exp(-input$plot_brush$ymax), scientific=T, digits=3), " ; " , format(exp(-input$plot_brush$ymin), scientific=T, digits=3),"]", "     <br/>    ", "<i>log2(FoldChange)</i>: ", "[", round(input$plot_brush$xmin, 3), " ; ", round(input$plot_brush$xmax, 3),"]<br/><br/>", sep=""))
      }}
  })
  
  output$table <- DT::renderDT({
    message("Rendering datatable")
    if (!is.na(preprocessedData()[1,1])){
      processedData()[processedData()$selected==T,-c("selected","minuslog10")]}
    else {
      processedData()
    }})
  
  output$volcanoPlot <- renderPlot ({
    message("Rendering volcano plot")
    if (!is.na(preprocessedData()[1,1])){
      plot()}
    else {
      return(emptyPlot)
    }
  })
  
  #--- GO Enrichment ---#
  
  output$goBarplot <- renderPlot({
    if (!is.null(goOraResults()$result)) {
      barplot(goOraResults()$result)
    } else {
      emptyPlot
    }
  })
  
  output$goDotplot <- renderPlot({
    if (!is.null(goOraResults()$result)) {
      dotplot(goOraResults()$result)
    } else {
      emptyPlot
    }
  })
  
  output$goNetplot <- renderPlot({
    if (!is.null(goOraResults()$result)) {
      go_enrich <- pairwise_termsim(goOraResults()$result)
      cnetplot(go_enrich, layout = "kk", showCategory = 15)
    } else {
      emptyPlot
    }
  })
  
  output$goEmapPlot <- renderPlot({
    if (!(is.null(goGseaResults()))) {
      go_enrich <- pairwise_termsim(goOraResults()$result)
      emapplot(go_enrich, layout = "kk", showCategory = 15)
    } else { emptyPlot }
  })
  
  output$goTable <- DT::renderDT({
    if (!(is.null(goOraResults()$result))) {
      df = as.data.frame(goOraResults()$result)
      for (i in 1:length(df$ID)) {
        df$ID[i] = paste('<a href="https://www.ebi.ac.uk/QuickGO/term/', df$ID[i], '" target="_blank">', df$ID[i], '</a>', sep="")
      }
      if (nrow(df) != 0) {
        datatable(df, escape = F)
      }
      else {
        emptyTableGo}
    } else {
      emptyTableGo
    }
  })
  
  output$goUpsetplot <- renderPlot({
    if (!(is.null(goOraResults()$result))) {
      upsetplot(goOraResults()$result)
    } else {
      emptyPlot
    }
    
  })
  
  output$goGseaUpsetplot <- renderPlot({
    if (!(is.null(goGseaResults()))) {
      upsetplot(goGseaResults())
    } else { emptyPlot }
  })
  
  output$goGseaNetplot <- renderPlot({
    if (!(is.null(goGseaResults()))) {
      go_enrich <- pairwise_termsim(goGseaResults())
      cnetplot(go_enrich, layout = "kk", showCategory = 15)
    } else { emptyPlot }
  })
  
  output$goGseaEmapPlot <- renderPlot({
    if (!(is.null(goGseaResults()))) {
      go_enrich <- pairwise_termsim(goGseaResults())
      emapplot(go_enrich, layout = "kk", showCategory = 15)
    } else { emptyPlot }
  })
  
  output$goGseaDotplot <- renderPlot({
    if (!(is.null(goGseaResults()))) {
      dotplot(goGseaResults())
    } else { emptyPlot }
  })
  
  output$goGseaEnrichmentPlot <- renderPlot({
    if (!(is.null(goGseaResults()))) {
      gseaplot2(goGseaResults(), by = "all", geneSetID = input$goGSEA)
      } else { emptyPlot }
    })
  
  output$goGseaRidgeplot <- renderPlot({
    if (!(is.null(goGseaResults()))) {
      ridgeplot(goGseaResults(), showCategory = 13)
    } else { emptyPlot }
  })
  
  output$goGseaTable <- DT::renderDT({
    if (!(is.null(goGseaResults()))) {
      df = as.data.frame(goGseaResults())
      for (i in 1:length(df$ID)) {
        df$ID[i] = paste('<a href="https://www.ebi.ac.uk/QuickGO/term/', df$ID[i], '" target="_blank">', df$ID[i], '</a>', sep="")
      }
      if (nrow(df) != 0) {
        datatable(df, escape = F)
      }
      else {
        emptyTableGo
        }
    } else {
      emptyTableGo
    }
  })
  
  output$gosigOfnodesplot <- renderPlot({
    if (!is.null(sigGO())) {
    showSigOfNodes(sigGO(), score(sigGOresult()), firstSigNodes = 3)
    } else {
      emptyPlot
    }
    
  })
  
  #--- Pathway Enrichment ---#
    
    output$pathway_gsea_table <- DT::renderDT ({
      message("Rendering Pathway DataTable")
      if (!is.null(pathwayGSEAEnrichment())){
        df = as.data.frame(pathwayGSEAEnrichment()$enrichment)
        for (i in 1:length(df$ID)) {
          if (pathwayGSEAEnrichment()$parameters$DB=="KEGG"){
            df$ID[i] = paste('<a href="https://www.kegg.jp/entry/', df$ID[i], '" target="_blank">', df$ID[i], '</a>', sep="")
          } else if (pathwayGSEAEnrichment()$parameters$DB=="Reactome"){
            df$ID[i] = paste('<a href="https://reactome.org/PathwayBrowser/#/', df$ID[i], '" target="_blank">', df$ID[i], '</a>', sep="")
          }
        }
        if (nrow(df) != 0) {
          datatable(df, escape = F)
        }
        else {
          emptyTable2}
       }
      else {
        emptyTable2
      }
    })
    
    output$pathway_ora_table <- DT::renderDT ({
      message("Rendering Pathway DataTable")
      if (!is.null(pathwayORAEnrichment())){
        df = as.data.frame(pathwayORAEnrichment()$enrichment)
        for (i in 1:length(df$ID)) {
          if (pathwayORAEnrichment()$parameters$DB=="KEGG"){
            df$ID[i] = paste('<a href="https://www.genome.jp/pathway/', df$ID[i], '" target="_blank">', df$ID[i], '</a>', sep="")
          } else if (pathwayORAEnrichment()$parameters$DB=="Reactome"){
            df$ID[i] = paste('<a href="https://reactome.org/PathwayBrowser/#/', df$ID[i], '" target="_blank">', df$ID[i], '</a>', sep="")
          }
        }
        if (nrow(df) != 0) {
          datatable(df, escape = F)
        }
        else {
          emptyTable2}
      }
      else {
        emptyTable2
      }
    })
    
    output$pathway_gsea_upsetplot <- renderPlot ({
      message("Rendering Pathway Upset plot")
      if (!is.null(pathwayGseaUpsetPlot())){
        pathwayGseaUpsetPlot()}
      else {
      }
    })
    
    output$pathway_ora_upsetplot <- renderPlot ({
      message("Rendering Pathway Upset plot")
      if (!is.null(pathwayORAUpsetPlot())){
        pathwayORAUpsetPlot()}
      else {
      }
    })
    
    output$pathway_gsea_dotplot <- renderPlot ({
      message("Rendering Pathway Dot plot")
      if (!is.null(pathwayGSEADotPlot())){
        pathwayGSEADotPlot()}
      else {
      }
    })
    
    output$pathway_ora_dotplot <- renderPlot ({
      message("Rendering Pathway Dot plot")
      if (!is.null(pathwayORADotPlot())){
        pathwayORADotPlot()}
      else {
      }
    })
      
    output$pathway_ora_treeplot <- renderPlot ({
      message("Rendering Pathway Dot plot")
      if (!is.null(pathwayORADotPlot())){
        pathwayORATreePlot()}
      else {
      }
    })
    
    output$pathway_ora_cnetplot <- renderPlot ({
      message("Rendering Pathway Gene-Concept Network plot")
      if (!is.null(pathwayORACnetPlot())){
        pathwayORACnetPlot()}
      else {
      }
    })  
    
    output$pathway_ora_emapplot <- renderPlot ({
      message("Rendering Pathway Enrichment map plot")
      if (!is.null(pathwayORAEmapPlot())){
        pathwayORAEmapPlot()}
      else {
      }
    })  
    
    output$pathway_gsea_emapplot <- renderPlot ({
      message("Rendering Pathway Enrichment map plot")
      if (!is.null(pathwayGSEAEmapPlot())){
        pathwayGSEAEmapPlot()}
      else {
      }
    })  
    
    output$pathway_gsea_treeplot <- renderPlot ({
      message("Rendering Pathway Dot plot")
      if (!is.null(pathwayGSEATreePlot())){
        pathwayGSEATreePlot()}
      else {
      }
    })
    
    output$pathway_gsea_cnetplot <- renderPlot ({
      message("Rendering Pathway Gene-Concept Network plot")
      if (!is.null(pathwayGSEACnetPlot())){
        pathwayGSEACnetPlot()}
      else {
      }
    })  
    
    output$pathway_gsea_enrichplot <- renderPlot ({
      message("Rendering Pathway GSEA plot")
      if (!is.null(pathwayGseaPlot()) & pathwayGSEAEnrichment()$parameters$analysis=="GSEA"){
        pathwayGseaPlot()}
      else {
      }
  })
    
  # -----------------------------------------
  # Save and export functions
  # -----------------------------------------

    output$exportOptions <- renderUI({
      
      # 1) Condition générale : on n'affiche l'UI que si on a bien des données
      if ((!is.null(preprocessedData()) && !is.na(preprocessedData()[1,1])) 
          || (!is.null(goOraResults()$result)) 
          || (!is.null(goGseaResults()))
          || (!is.null(pathwayORAEnrichment())) 
          || (!is.null(pathwayGSEAEnrichment()))) {
        
        # 2) Construction du vecteur choices
        choices <- character(0)
        if (!is.null(preprocessedData()) && !is.na(preprocessedData()[1,1])) {
          choices <- "Whole Data Inspection (VolcanoPlot + Table)"
        }
        if (!is.null(goOraResults()$result))    choices <- c(choices, "GO ORA")
        if (!is.null(goGseaResults()))          choices <- c(choices, "GO GSEA")
        if (!is.null(pathwayORAEnrichment()))   choices <- c(choices, "Pathway ORA")
        if (!is.null(pathwayGSEAEnrichment()))  choices <- c(choices, "Pathway GSEA")
        
        # 3) Construction dynamique de la liste d'éléments à afficher
        ui_elems <- list()
        
        if (length(choices) == 1 && choices == "Whole Data Inspection (VolcanoPlot + Table)") {
          # Cas où on n'a QUE la Whole Data Inspection
          ui_elems <- tagList(
            ui_elems,
            HTML("<em>No enrichment was carried out. 
              The current export will only contain Whole Data Inspection.</em><br/>")
          )
        } else {
          # Cas général où on affiche les checkboxes
          ui_elems <- tagList(
            ui_elems,
            checkboxGroupInput("export_choices", "Analysis to export:", choices = choices)
          )
        }
        
        # 4) On ajoute toujours le bouton de téléchargement
        ui_elems <- tagList(
          ui_elems,
          downloadButton("exportReport", "Download report")
        )
        
        # 5) On retourne la totalité
        ui_elems
        
      } else {
        # Pas de données : message informatif
        HTML("<em>Only already-launched analyses can be exported.</em><br/>")
      }
    })

    
    
    output$exportReport <- downloadHandler(
      filename = function() {
        paste0("HEATraN_results_", Sys.Date(), ".html")
      },
      content = function(file) {
        tempReport <- file.path(tempdir(), "template.Rmd")
        file.copy("www/template.Rmd", tempReport, overwrite = TRUE)
        
        # Rendre le document R Markdown
        rmarkdown::render(
          input = tempReport,
          output_file = file,
          params = list(
            # General parameters
            nGenes = ifelse(!is.null(input$nGenesExport), input$nGenesExport, 20),
            organismInfo = input$species,
            
            # Basic outputs - Capture des objets réactifs
            InfoSelect = tryCatch({
              if (!is.null(preprocessedData()) && !is.na(preprocessedData()[1,1])) {
                if (is.null(input$plot_brush)) {
                  paste("Current selection - p-value: [0 ; ", input$pval,"] log2(FoldChange): [-", input$Log2FC, " ; ", input$Log2FC,"]", sep="")
                } else {
                  paste("Current selection - p-value: [", format(exp(-input$plot_brush$ymax), scientific=T, digits=3), " ; " , format(exp(-input$plot_brush$ymin), scientific=T, digits=3),"] log2(FoldChange): [", round(input$plot_brush$xmin, 3), " ; ", round(input$plot_brush$xmax, 3),"]", sep="")
                }
              } else NULL
            }, error = function(e) NULL),
            
            volcanoPlot = tryCatch({
              if (!is.null(preprocessedData()) && !is.na(preprocessedData()[1,1])) {
                p <- plot()  
                p  
              } else NULL
            }, error = function(e) NULL),
            
            table = tryCatch({
              if (!is.null(preprocessedData()) && !is.na(preprocessedData()[1,1])) {
                data <- processedData()
                if (!is.null(data) && nrow(data) > 0) {
                  data[data$selected == TRUE, -c("selected", "minuslog10")]
                } else NULL
              } else NULL
            }, error = function(e) NULL),
            
            # GO ORA outputs - Capture des fonctions réactives
            goBarplot = tryCatch({
              if (!is.null(goOraResults()) && !is.null(goOraResults()$result)) {
                barplot(goOraResults()$result)
              } else NULL
            }, error = function(e) NULL),
            
            goDotplot = tryCatch({
              if (!is.null(goOraResults()) && !is.null(goOraResults()$result)) {
                dotplot(goOraResults()$result)
              } else NULL
            }, error = function(e) NULL),
            
            goEmapPlot = tryCatch({
              if (!is.null(goOraResults()) && !is.null(goOraResults()$result)) {
                go_enrich <- pairwise_termsim(goOraResults()$result)
                emapplot(go_enrich, layout = "kk", showCategory = 15)
              } else NULL
            }, error = function(e) NULL),
            
            goNetplot = tryCatch({
              if (!is.null(goOraResults()) && !is.null(goOraResults()$result)) {
                go_enrich <- pairwise_termsim(goOraResults()$result)
                cnetplot(go_enrich, layout = "kk", showCategory = 15)
              } else NULL
            }, error = function(e) NULL),
            
            goTable = tryCatch({
              if (!is.null(goOraResults()) && !is.null(goOraResults()$result)) {
                as.data.frame(goOraResults()$result)
              } else NULL
            }, error = function(e) NULL),
            
            goUpsetplot = tryCatch({
              if (!is.null(goOraResults()) && !is.null(goOraResults()$result)) {
                upsetplot(goOraResults()$result)
              } else NULL
            }, error = function(e) NULL),
            
            gosigOfnodesplot = tryCatch({
              if (!is.null(sigGO()) && !is.null(sigGOresult())) {
                img_file <- tempfile(fileext = ".png")
                png(filename = img_file, width = 800, height = 600)
                showSigOfNodes(sigGO(), score(sigGOresult()), firstSigNodes = 3)
                dev.off()
                list(src = img_file, contentType = "image/png", width = "100%")
              } else {
                NULL
              }
            }, error = function(e) {
              NULL
            }),
            
            # GO GSEA outputs
            goGseaUpsetplot = tryCatch({
              if (!is.null(goGseaResults())) {
                upsetplot(goGseaResults())
              } else NULL
            }, error = function(e) NULL),
            
            goGseaNetplot = tryCatch({
              if (!is.null(goGseaResults())) {
                go_enrich <- pairwise_termsim(goGseaResults())
                cnetplot(go_enrich, layout = "kk", showCategory = 15)
              } else NULL
            }, error = function(e) NULL),
            
            goGseaEmapPlot = tryCatch({
              if (!is.null(goGseaResults())) {
                go_enrich <- pairwise_termsim(goGseaResults())
                emapplot(go_enrich, layout = "kk", showCategory = 15)
              } else NULL
            }, error = function(e) NULL),
            
            goGseaDotplot = tryCatch({
              if (!is.null(goGseaResults())) {
                dotplot(goGseaResults())
              } else NULL
            }, error = function(e) NULL),
            
            goGseaRidgeplot = tryCatch({
              if (!is.null(goGseaResults())) {
                ridgeplot(goGseaResults(), showCategory = 13)
              } else NULL
            }, error = function(e) NULL),
            
            goGseaTable = tryCatch({
              if (!is.null(goGseaResults())) {
                as.data.frame(goGseaResults())
              } else NULL
            }, error = function(e) NULL),
            
            # Pathway ORA outputs
            pathway_ora_table = tryCatch({
              if (!is.null(pathwayORAEnrichment())) {
                as.data.frame(pathwayORAEnrichment()$enrichment)
              } else NULL
            }, error = function(e) NULL),
            
            pathway_ora_upsetplot = tryCatch({
              if (!is.null(pathwayORAEnrichment())) {
                p <- pathwayORAUpsetPlot()
                p
              } else NULL
            }, error = function(e) NULL),
            
            pathway_ora_emapplot = tryCatch({
              if (!is.null(pathwayORAEmapPlot())) {
                p <- pathwayORAEmapPlot()
                p
              } else NULL
            }, error = function(e) NULL),
            
            pathway_ora_dotplot = tryCatch({
              if (!is.null(pathwayORAEnrichment())) {
                p <- pathwayORADotPlot()
                p
              } else NULL
            }, error = function(e) NULL),
            
            pathway_ora_treeplot = tryCatch({
              if (!is.null(pathwayORAEnrichment())) {
                p <- pathwayORATreePlot()
                p
              } else NULL
            }, error = function(e) NULL),
            
            pathway_ora_cnetplot = tryCatch({
              if (!is.null(pathwayORAEnrichment())) {
                p <- pathwayORACnetPlot()
                p
              } else NULL
            }, error = function(e) NULL),
            
            # Pathway GSEA outputs
            pathway_gsea_table = tryCatch({
              if (!is.null(pathwayGSEAEnrichment())) {
                as.data.frame(pathwayGSEAEnrichment()$enrichment)
              } else NULL
            }, error = function(e) NULL),
            
            pathway_gsea_upsetplot = tryCatch({
              if (!is.null(pathwayGSEAEnrichment())) {
                p <- pathwayGseaUpsetPlot()
                p
              } else NULL
            }, error = function(e) NULL),
            
            pathway_gsea_dotplot = tryCatch({
              if (!is.null(pathwayGSEAEnrichment())) {
                p <- pathwayGSEADotPlot()
                p
              } else NULL
            }, error = function(e) NULL),
            
            pathway_gsea_treeplot = tryCatch({
              if (!is.null(pathwayGSEAEnrichment())) {
                p <- pathwayGSEATreePlot()
                p
              } else NULL
            }, error = function(e) NULL),
            
            pathway_gsea_cnetplot = tryCatch({
              if (!is.null(pathwayGSEAEnrichment())) {
                p <- pathwayGSEACnetPlot()
                p
              } else NULL
            }, error = function(e) NULL),
            
            pathway_gsea_emapplot = tryCatch({
              if (!is.null(pathwayGSEAEmapPlot())) {
                p <- pathwayGSEAEmapPlot()
                p
              } else NULL
            }, error = function(e) NULL),
            
            # Complex data objects pour génération de plots individuels
            goOraResults = goOraResults(),
            goGseaResults = goGseaResults(),
            pathwayGSEAEnrichment = pathwayGSEAEnrichment(),
            pathwayORAEnrichment = pathwayORAEnrichment(),
            preprocessedData = preprocessedData(),
            processedData = processedData(),
            export_choices     = input$export_choices,
            includeGSEAPlots   = input$includeGSEAPlots
            
            # Export sections flags
            #exportSections = input$exportSections
          ),
          knit_root_dir = normalizePath(getwd()),
          envir = new.env(parent = globalenv())
        )
      })
    
  
  # Functions launched when the application is closed
  onStop(function() {
    if (config$FILE$clear_cache == "True") {
      if ((dir_size("./out") / (1024^2)) > config$FILE$max_cache_mb) {
        message("Cache exceeding the authorised limit: cleaning cache...")
        files_to_delete <- list.files(path = "./out", full.names = T)
        unlink(files_to_delete)
      }
    }
    message("Thanks for using HEATraN!")
  })
  
  combinedPathwayResults <- reactive({
    ora_results <- pathwayORAEnrichment()
    gsea_results <- pathwayGSEAEnrichment()
    
    common_cols <- c("ID", "Description", "pvalue", "p.adjust")
    
    # Si les deux analyses sont disponibles
    if (!is.null(ora_results) && !is.null(gsea_results)) {
      
      ora_data <- ora_results$enrichment@result[, common_cols, drop = FALSE]
      gsea_data <- gsea_results$enrichment@result[, common_cols, drop = FALSE]
      
      # Ajouter une colonne pour identifier le type d'analyse
      ora_data$analysis_type <- "ORA"
      gsea_data$analysis_type <- "GSEA"
      
      combined_data <- rbind(ora_data, gsea_data)
      
      # Enlever les doublons basés sur l'ID du pathway
      #combined_data <- combined_data[!duplicated(combined_data$ID), ]
      
      if (ora_results$parameters$DB == "KEGG") {
        choices <- setNames(combined_data$ID, 
                            paste0(combined_data$Description, " (", combined_data$analysis_type, ")"))
      } else { # Reactome
        choices <- setNames(combined_data$Description, 
                            paste0(combined_data$Description, " (", combined_data$analysis_type, ")"))
      }
      
      return(list(
        choices = choices,
        combined_data = combined_data,
        ora_results = ora_results,
        gsea_results = gsea_results
      ))
      
    } else if (!is.null(ora_results)) {
      # Seulement ORA disponible
      ora_data <- ora_results$enrichment@result[, common_cols, drop = FALSE]
      ora_data$analysis_type <- "ORA"
      
      if (ora_results$parameters$DB == "KEGG") {
        choices <- setNames(ora_data$ID, paste0(ora_data$Description, " (ORA)"))
      } else {
        choices <- setNames(ora_data$Description, paste0(ora_data$Description, " (ORA)"))
      }
      
      return(list(
        choices = choices,
        combined_data = ora_data,
        ora_results = ora_results,
        gsea_results = NULL
      ))
      
    } else if (!is.null(gsea_results)) {
      # Seulement GSEA disponible
      gsea_data <- gsea_results$enrichment@result[, common_cols, drop = FALSE]
      gsea_data$analysis_type <- "GSEA"
      
      if (gsea_results$parameters$DB == "KEGG") {
        choices <- setNames(gsea_data$ID, paste0(gsea_data$Description, " (GSEA)"))
      } else {
        choices <- setNames(gsea_data$Description, paste0(gsea_data$Description, " (GSEA)"))
      }
      
      return(list(
        choices = choices,
        combined_data = gsea_data,
        ora_results = NULL,
        gsea_results = gsea_results
      ))
    }
    
    return(NULL)
  })
  
  # Mettre à jour les choix de pathways
  output$pathwayImage <- renderImage({
    combined_results <- combinedPathwayResults()
    selected_pathway <- input$pathway
    
    if (!is.null(combined_results) && selected_pathway!="") {
      
      pathway_data <- combined_results$combined_data

      if (nrow(pathway_data) > 0) {
        if ((!is.null(combined_results$ora_results$parameters$DB) && combined_results$ora_results$parameters$DB == "KEGG")
            || (!is.null(combined_results$gsea_results$parameters$DB) && combined_results$gsea_results$parameters$DB == "KEGG")) {
          pathway_row <- which(pathway_data$ID == selected_pathway)
        } else {
          pathway_row <- which(pathway_data$Description == selected_pathway)
        }
        
        if (length(pathway_row) > 0) {
          analysis_type <- pathway_data$analysis_type[pathway_row[1]]
          combined_paths = c(combined_results$ora_results$pathway_images, combined_results$gsea_results$pathway_images)

        if (analysis_type == "ORA" && !is.null(combined_results$ora_results$pathway_images)) {
            image_path <- combined_paths[[pathway_row[1]]]
          } else if (analysis_type == "GSEA" && !is.null(combined_results$gsea_results$pathway_images)) {
            image_path <- combined_paths[[pathway_row[1]]]
          } else {
            image_path <- NULL
          }
          
          if (!is.null(image_path) && file.exists(image_path)) {
            return(list(src = image_path,
                        contentType = "image/png",
                        width = "100%",
                        alt = paste("Pathway:", selected_pathway, "(", analysis_type, ")")))
          } else {
            print(image_path)
          }
        }
      }
    } else if (!is.null(combined_results)) {
      ggsave("./out/empty.png", plot = emptyPathwayPlot)
      return(list(src = "./out/empty.png", contentType = 'image/png', width="100%", alt="Please select a pathway"))
    }
    
    # Image par défaut si aucun pathway n'est trouvé
    return(list(src = "", alt = "No pathway image available"))
  }, deleteFile = FALSE)
  
  observeEvent(combinedPathwayResults(), {
    combined_results <- combinedPathwayResults()
    
    if (!is.null(combined_results)) {
      # Afficher l'onglet "View Pathway"
      showTab(inputId = "pathway_tabs", target = "View Pathway", session = session)
      
      # Mettre à jour les choix
      updateSelectInput(session, "pathway", 
                        choices = combined_results$choices,
                        selected = names(combined_results$choices)[1])
    }
  })
  

}

