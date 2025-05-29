# Developed by LESAGE Louison (@thelokj).
# louison.lesage@univ-rouen.fr
# Student at Rouen Normandy University
# University project 2024-2025
# Last updated : 18/11/2024
# HEATraN version 0.2.0-a.5

emptyTable <- data.frame(Gene=NA, Log2FC=NA, p_value=NA)
brushInfo <- reactiveVal(NULL)
selectionMode <- reactiveVal("None")
preprocessedData <- reactiveVal(emptyTable)
requiredNames <- c("GeneName", "GeneID", "baseMean", "Log2FC", "pval", "padj")

function(input, output, session) {
  
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
    if (FALSE %in% c(requiredNames %in% colnames(dataToPreprocess))){
      shinyalert(html = TRUE, "Please select the variables", "Wrong column name", type = "info", confirmButtonCol = "#7e3535",
                 text = tagList(
                   selectInput("GeneNameCol", "Gene Name", colnames(dataToPreprocess)),
                   selectInput("GeneIDCol", "Gene ID", colnames(dataToPreprocess)),
                   selectInput("BaseMeanCol", "Base mean", colnames(dataToPreprocess)),
                   selectInput("Log2FCCol", "Log2(FoldChange)", colnames(dataToPreprocess)),
                   selectInput("pvalCol", "p-value", colnames(dataToPreprocess)),
                   selectInput("padjCol", "Adjusted p-value", colnames(dataToPreprocess)),
                   HTML("Next time, use these names to import direcly the table : <i>GeneName, GeneID, baseMean, Log2FC, pval, padj</i><br>"),
                 ),
                 callbackR = function(finished) {
                   if(finished) {
                     if (length(unique(c(input$GeneNameCol, input$GeneIDCol, input$BaseMeanCol,input$Log2FCCol,input$pvalCol, input$padjCol)))==6){
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
        df <- df[df$Log2FC>=Zoom()$coords[1]&df$Log2FC<=Zoom()$coords[2]&(-log(df$pval))>=Zoom()$coords[3]&(-log(df$pval))<=Zoom()$coords[4],]
      }
      # Update the selected points according to the selection mode
      if (selectionMode() == "Brush") {
        df$selected <- ifelse(df$GeneID%in%brushedPoints(df, brushInfo()())$GeneID, "TRUE", "FALSE")
      } else if (selectionMode() == "Sliders"){
        df$selected <- ifelse((df$Log2FC>input$Log2FC&df$pval<input$pval)|(df$Log2FC<(-input$Log2FC)&df$pval<input$pval), "TRUE", "FALSE")
      }
      return(df)
    } else {
      return(emptyTable)
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
      ylab("-log10(p-value)") +
      theme(text = element_text(size = 14))
    return(plot)})
  
  # -----------------------------------------
  # User event
  # -----------------------------------------
  
  # Reactive function controling the selection mode
  # As shiny do not allow to disable downloadButton, disable it with shinyJS
  # Disable also other buttons with shinyJS to standardize rendering
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
  
  # -----------------------------------------
  # UI output
  # -----------------------------------------
  
  # Rendering "Current selection" section
  output$InfoSelect <- renderUI({
    if (is.null(preprocessedData) || is.na(preprocessedData()[1,1])){
      HTML(paste("<b style='color:	#FF0000'>Data required for selection</b><br/><br/>", sep=""))
    }
    else {
      if (is.null(input$plot_brush)){
        HTML(paste("<b>Current selection</b><br/><i>p-value</i>: [0 ; ", input$pval,"]", "    <br/>    ", "<i>log2(FoldChange)</i>: ", "[-", input$Log2FC, " ; ", input$Log2FC,"]<br/><br/>", sep=""))
      }
      else {
        HTML(paste("<b>Current selection</b><br/><i>p-value</i>: [", format(exp(-input$plot_brush$ymax), scientific=TRUE, digits=3), " ; " , format(exp(-input$plot_brush$ymin), scientific=TRUE, digits=3),"]", "     <br/>    ", "<i>log2(FoldChange)</i>: ", "[", round(input$plot_brush$xmin, 3), " ; ", round(input$plot_brush$xmax, 3),"]<br/><br/>", sep=""))
      }}
  })
  
  # Rendering datatable
  output$table <- DT::renderDT({
    message("Rendering datatable")
    if (!is.na(preprocessedData()[1,1])){
      processedData()[processedData()$selected==TRUE,-c("selected","minuslog10")]}
    else {
      processedData()
    }})
  
  # Rendering volcano plot
  output$volcanoPlot <- renderPlot ({
    message("Rendering volcano plot")
    if (!is.na(preprocessedData()[1,1])){
      plot()}
    else {
    }
  })
  # GO Analysis reactive values
  # GO Analysis reactive values
  geneGroups <- reactiveValues(
    up = NULL,
    down = NULL,
    both = NULL
  )
  topGOdataList <- reactiveValues(up=NULL, down=NULL, both=NULL)
  
  goResults <- reactiveVal(NULL)
  
  
  observeEvent(input$go_analysisButton, {
    req(preprocessedData())
    
    withProgress(message = 'Running GO analysis...', {
      
      df <- preprocessedData()
      
      # Liste des gènes pour l'univers
      original_gene_list <- df$Log2FC
      names(original_gene_list) <- df$GeneID
      gene_list <- sort(na.omit(original_gene_list), decreasing = TRUE)
      
      # Sélection des gènes selon le seuil de p-value ajustée
      sig_genes_df <- subset(df, padj < input$go_pval)
      genes <- sig_genes_df$Log2FC
      names(genes) <- sig_genes_df$GeneID
      genes <- na.omit(genes)
      
      up_genes <- names(genes)[genes > 0]
      down_genes <- names(genes)[genes < 0]
      
      both_genes <- union(up_genes, down_genes)
      geneGroups$up <- up_genes
      geneGroups$down <- down_genes
      geneGroups$both <- both_genes
      # Sélection de l'ontologie GO et de l'organisme
      ontology <- input$inputGO
      organism <- input$organism
      universe <- names(gene_list)
      library(organism, character.only = TRUE)
      
      # Méthode d’analyse : ORA ou GSEA
      if ("Gene Set Enrichment Analysis (GSEA)" %in% input$go_analysisMethodChoice) {
        # GSEA à implémenter ici
        
      } else {
        # ORA selon les sélections
        analyze_up <- "Over expressed DEG" %in% input$go_oraChoice
        analyze_down <- "Under expressed DEG" %in% input$go_oraChoice
        
        if (!analyze_up && !analyze_down) {
          shinyalert("Erreur", "Veuillez sélectionner au moins un type de gène (Over ou Under expressed DEG).", type = "error")
          return(NULL)
        }
        
        result_up <- NULL
        result_down <- NULL
        result_both <- NULL  # Ajout des résultats pour Both-regulated
        
        if (analyze_up && length(up_genes) > 0) {
          result_up <- enrichGO(
            gene = up_genes,
            universe = universe,
            OrgDb = organism,
            keyType = "ENSEMBL",
            readable = TRUE,
            ont = ontology,
            pvalueCutoff = input$go_pval,
            qvalueCutoff = input$qvalueCutoff
          )
        }
        
        if (analyze_down && length(down_genes) > 0) {
          result_down <- enrichGO(
            gene = down_genes,
            universe = universe,
            OrgDb = organism,
            keyType = "ENSEMBL",
            readable = TRUE,
            ont = ontology,
            pvalueCutoff = input$go_pval,
            qvalueCutoff = input$qvalueCutoff 
          )
        }
        
        # Pour les résultats Both-regulated, tu pourrais combiner les gènes up et down (généralement ce n'est pas standard en GO, mais ça peut être utile selon ton cas d'usage)
        both_genes <- union(up_genes, down_genes)  # Union des gènes up et down
        if (length(both_genes) > 0) {
          result_both <- enrichGO(
            gene = both_genes,
            universe = universe,
            OrgDb = organism,
            keyType = "ENSEMBL",
            readable = TRUE,
            ont = ontology,
            pvalueCutoff = input$go_pval,
            qvalueCutoff = input$qvalueCutoff 
          )
        }
        
        # Storer les résultats dans goResults
        goResults(list(up = result_up, down = result_down, both = result_both))
      }
    })
    # Créer la liste binaire pour topGO
    allGenesBinary <- as.integer(names(gene_list) %in% both_genes)
    names(allGenesBinary) <- names(gene_list)
    
    # Même chose pour up et down
    upGenesBinary <- as.integer(names(gene_list) %in% up_genes)
    names(upGenesBinary) <- names(gene_list)
    
    downGenesBinary <- as.integer(names(gene_list) %in% down_genes)
    names(downGenesBinary) <- names(gene_list)
    
    # Créer les objets topGOdata
    GOdata_up <- new("topGOdata",
                     ontology = ontology,
                     allGenes = upGenesBinary,
                     geneSelectionFun = function(x)(x < 0.05),
                     annot = annFUN.org,
                     mapping = organism,
                     ID = "ensembl")
    
    GOdata_down <- new("topGOdata",
                       ontology = ontology,
                       allGenes = downGenesBinary,
                       geneSelectionFun = function(x)(x < 0.05),
                       annot = annFUN.org,
                       mapping = organism,
                       ID = "ensembl")
    
    GOdata_both <- new("topGOdata",
                       ontology = ontology,
                       allGenes = allGenesBinary,
                       geneSelectionFun = function(x)(x < 0.05),
                       annot = annFUN.org,
                       mapping = organism,
                       ID = "ensembl")
    
    # Stocker dans reactiveValues
    topGOdataList$up <- GOdata_up
    topGOdataList$down <- GOdata_down
    topGOdataList$both <- GOdata_both
    
  })
  
  # === PLOTS ORA ===
  # Barplots
  # output$goBarplotUp <- renderPlot({
  #   req(goResults()$up)
  #   barplot(goResults()$up,
  #           drop = TRUE,
  #           showCategory = input$topCategoriesUp,
  #           title = "GO Biological Pathways (Up-regulated)",
  #           font.size = 8)
  # })
  
  output$goBarplotDown <- renderPlot({
    if (!is.null(topGOdataList$down) && !is.null(geneGroups$down) && length(geneGroups$down) > 0) {
      barplot(topGOdataList$down, showCategory = 10)
    } else {
      plot.new()
      text(0.5, 0.5, "No data available for down-regulated genes")
    }
  })
  
  output$goBarplotBoth <- renderPlot({
    if (!is.null(topGOdataList$both) && !is.null(geneGroups$both) && length(geneGroups$both) > 0) {
      barplot(topGOdataList$both, showCategory = 10)
    } else {
      plot.new()
      text(0.5, 0.5, "No data available for both-regulated genes")
    }
  })
  
  
  # Dotplots
  output$goDotplotUp <- renderPlot({
    req(goResults()$up)
    dotplot(goResults()$up)
  })
  
  
  output$goDotplotDown <- renderPlot({
    req(goResults()$down)
    dotplot(goResults()$down)
  })
  
  output$goDotplotBoth <- renderPlot({
    req(goResults()$both)
    dotplot(goResults()$both)
  })
  
  # Network plots
  output$goNetplotUp <- renderPlot({
    req(goResults()$up)
    go_enrich <- pairwise_termsim(goResults()$up)
    # emapplot(go_enrich, layout = "kk", showCategory = 15)
    emapplot(go_enrich, layout.params = list(layout = "kk"), showCategory = 15)
    
  })
  
  output$goNetplotDown <- renderPlot({
    req(goResults()$down)
    go_enrich <- pairwise_termsim(goResults()$down)
    emapplot(go_enrich, layout = "kk", showCategory = 15)
  })
  
  output$goNetplotBoth <- renderPlot({
    req(goResults()$both)
    go_enrich <- pairwise_termsim(goResults()$both)
    emapplot(go_enrich, layout = "kk", showCategory = 15)
  })
  
  # === TABLES ORA ===
  output$goTableUp <- DT::renderDT({
    req(goResults()$up)
    df <- as.data.frame(goResults()$up)
    for (i in 1:length(df$ID)) {
      df$ID[i] <- paste0('<a href="https://www.ebi.ac.uk/QuickGO/term/', df$ID[i],
                         '" target="_blank">', df$ID[i], '</a>')
    }
    
    DT::datatable(df, escape = FALSE, options = list(pageLength = 10))
  })
  
  
  
  output$goTableDown <- DT::renderDT({
    req(goResults()$down)
    df <- as.data.frame(goResults()$down)
    
    
    for (i in 1:length(df$ID)) {
      df$ID[i] <- paste0('<a href="https://www.ebi.ac.uk/QuickGO/term/', df$ID[i],
                         '" target="_blank">', df$ID[i], '</a>')
    }
    
    DT::datatable(df, escape = FALSE, options = list(pageLength = 10))
  })
  
  output$goTableBoth <- DT::renderDT({
    req(goResults()$both)
    df <- as.data.frame(goResults()$both)
    
    for (i in 1:length(df$ID)) {
      df$ID[i] <- paste0('<a href="https://www.ebi.ac.uk/QuickGO/term/', df$ID[i],
                         '" target="_blank">', df$ID[i], '</a>')
    }
    
    DT::datatable(df, escape = FALSE)
  })
  ################upset plot
  
  output$upsetplotUp <- renderPlot({
    req(goResults()$up)
    upsetplot(goResults()$up)
    
  })
  
  output$upsetplotDown <- renderPlot({
    req(goResults()$down)
    upsetplot(goResults()$down)
  })
  
  output$upsetplotBoth <- renderPlot({
    req(goResults()$both)
    upsetplot(goResults()$both)})
  # ########"
  # # Network plots
  # output$goNetplotUp <- renderPlot({
  #   req(goResults()$up)
  #   go_enrich <- pairwise_termsim(goResults()$up)
  #   # emapplot(go_enrich, layout = "kk", showCategory = 15)
  #   emapplot(go_enrich, layout.params = list(layout = "kk"), showCategory = 15)
  #   
  # })
  # 
  # output$goNetplotDown <- renderPlot({
  #   req(goResults()$down)
  #   go_enrich <- pairwise_termsim(goResults()$down)
  #   emapplot(go_enrich, layout = "kk", showCategory = 15)
  # })
  # 
  # output$goNetplotBoth <- renderPlot({
  #   req(goResults()$both)
  #   go_enrich <- pairwise_termsim(goResults()$both)
  #   emapplot(go_enrich, layout = "kk", showCategory = 15)
  # })
  
  # CNE plot 
  output$cneUp <- renderPlot({
    req(goResults()$up)
    cnetplot(goResults()$up, showCategory = 10, foldChange = geneGroups$up )
  })
  
  output$cneDown <- renderPlot({
    req(goResults()$up)
    cnetplot(goResults()$up, showCategory = 10, foldChange = geneGroups$down)
  })
  
  output$cneBoth <- renderPlot({
    req(goResults()$both)
    cnetplot(goResults()$both, showCategory = 10, foldChange = geneGroups$both)
  })
  ###########################################
  
  
  # output$goTreeplotUp <- renderPlot({
  #   req(goResults()$up)
  #   str(goResults()$up)
  #   go_up <- enrichplot::pairwise_termsim(goResults()$up)
  #   treeplot(go_up)
  # })
  # 
  # 
  # output$goTreemapUp <- renderPlot({
  #   req(goResults()$up)
  #   df <- as.data.frame(goResults()$up)
  #   if (nrow(df) == 0) {
  #     plot.new(); title("Aucun terme GO significatif")
  #     return()
  #   }
  #   df_clean =df[1:8,]
  #   df_clean$Description <- ifelse(nchar(df_clean$Description) > 10,
  #                                  paste0(substr(df_clean$Description, 1, 5), "..."),
  #                                  df_clean$Description)
  #   treemap::treemap(df_clean,
  #                    index = "Description",
  #                    vSize = "Count",
  #                    vColor = "p.adjust",
  #                    type = "value",
  #                    palette = "RdYlBu")
  # })
  # 
  #####graph : 
  # 
  # # Création des objets topGOdata pour chaque groupe
  # GOdata_up <- new("topGOdata",
  #                  ontology = "BP",
  #                  allGenes =  geneGroups$up ,
  #                  annot = annFUN.org,
  #                  mapping = organism,
  #                  ID = "ensembl")
  # 
  # GOdata_down <- new("topGOdata",
  #                    ontology = "BP",
  #                    allGenes = geneGroups$down,                
  #                    annot = annFUN.org,
  #                    mapping = organism,
  #                    ID = "ensembl")
  # 
  # GOdata_both <- new("topGOdata",
  #                    ontology = "BP",
  #                    allGenes = gene_list,
  #                    annot = annFUN.org,
  #                    mapping = organism,
  #                    ID = "ensembl")
  # 
  # # Stocker dans reactiveValues pour pouvoir les utiliser dans renderPlot
  # topGOdataList <- reactiveValues(
  #   up = GOdata_up,
  #   down = GOdata_down,
  #   both = GOdata_both
  # )
  # 
  
  # # # Fonction helper pour plot un graph topGO
  # output$plotTopGO <- renderPlot({
  #   req(topGOdataList[[input$groupGO]])
  #    
  #   result <- runTest(topGOdataList[[input$groupGO]], algorithm = "classic", statistic = "fisher")
  #   
  #   showSigOfNodes(topGOdataList[[input$groupGO]], score(result), firstSigNodes = input$nbTopGO)
  # })
  # 
  output$plotTopGOUp <- renderPlot({
    req(topGOdataList$up)
    
    result <- runTest(topGOdataList$up, 
                      algorithm = "classic", 
                      statistic = "fisher")
    
    showSigOfNodes(topGOdataList$up, 
                   score(result), 
                   firstSigNodes = input$nbTopGO)
  })
  
  output$plotTopGODown <- renderPlot({
    req(topGOdataList$down)
    
    result <- runTest(topGOdataList$down, 
                      algorithm = "classic", 
                      statistic = "fisher")
    
    showSigOfNodes(topGOdataList$down, 
                   score(result), 
                   firstSigNodes = input$nbTopGO)
  })
  
  output$plotTopGOBoth <- renderPlot({
    req(topGOdataList$both)
    
    result <- runTest(topGOdataList$both, 
                      algorithm = "classic", 
                      statistic = "fisher")
    
    showSigOfNodes(topGOdataList$both, 
                   score(result), 
                   firstSigNodes = input$nbTopGO)
  })
  
  
  
  # plotTopGO <- function(GOdata, firstSigNodes) {
  #   message("Running runTest...")
  #   scoreFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  #   message("Running showSigOfNodes...")
  #   #showSigOfNodes(GOdata, score(scoreFisher), firstSigNodes = firstSigNodes)
  #   showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 5)
  #   
  #   printGraph(GOdata)
  #   
  # }
  # 
  # # Rendu plot UP
  # output$topGOgraph_up <- renderPlot({
  #   req(topGOdataList$up)
  #   req(input$firstSigNodes)
  #   scoreFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  #   message("Running showSigOfNodes...")
  #   #showSigOfNodes(GOdata, score(scoreFisher), firstSigNodes = firstSigNodes)
  #   showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 5)
  #   
  #   message("Rendering topGOgraph_up...")
  #  # plotTopGO(topGOdataList$up, input$firstSigNodes)
  # })
  # 
  # Rendu plot DOWN
  # output$topGOgraph_down <- renderPlot({
  #   req(topGOdataList$down)
  #   req(input$firstSigNodes)
  #   message("Rendering topGOgraph_down...")
  #   plotTopGO(topGOdataList$down, input$firstSigNodes)
  #   
  # })
  # 
  # # Rendu plot BOTH
  # output$topGOgraph_both <- renderPlot({
  #   req(topGOdataList$both)
  #   req(input$firstSigNodes)
  #   message("Rendering topGOgraph_both...")
  #   plotTopGO(topGOdataList$both, input$firstSigNodes)
  # })
  # 
  
}