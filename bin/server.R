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

emptyTable <- data.frame(Gene=NA, Log2FC=NA, p_value=NA)
emptyTableGo <- data.frame(GO=NA, Description=NA, p_value=NA, q_value=NA)
emptyTable2 <- data.frame(Pathway=NA, p_value=NA, q_value=NA)
brushInfo <- reactiveVal(NULL)
pathwayEnrichment <- reactiveVal(NULL)
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
                    label = "<b style='font-family:Mulish,sans-serif;font-size:20pt;'>Please select at least one enriched pathways</b>"), 
                size = 6, color = "red3", fill = NA, label.color = NA) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = NA))

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
  
  observe ({
    message("Preprocessing data")
    dataToPreprocess <- importedData()
    if (F %in% c(requiredNames %in% colnames(dataToPreprocess))){
      shinyalert(html = T, "Please select the variables", "Wrong column name", type = "info", confirmButtonCol = "#7e3535",
                 text = tagList(
                   selectInput("GeneNameCol", "Gene Name", colnames(dataToPreprocess)),
                   selectInput("GeneIDCol", "Gene ID", colnames(dataToPreprocess)),
                   selectInput("BaseMeanCol", "Base mean", colnames(dataToPreprocess)),
                   selectInput("Log2FCCol", "Log2(FoldChange)", colnames(dataToPreprocess)),
                   selectInput("pvalCol", "p-value", colnames(dataToPreprocess)),
                   selectInput("padjCol", "Adjusted p-value", colnames(dataToPreprocess)),
                   checkboxInput("saveColNames", "use these column names in the future"),
                   HTML(paste("Next time, use these names to import direcly the table : <i>", paste(requiredNames, collapse=", "), "</i><br>"))
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
                   } else {return(F)}})
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
      if (!isT(all.equal(emptyTable, preprocessedData()))){
        df <- preprocessedData()
        updateSliderInput(session,'Log2FC',max=ceiling(max(abs(df$Log2FC))))
        # Adapt the plot limit when user zoom in it
        if (Zoom()$Zoomed==T){
          df <- df[df$Log2FC>=Zoom()$coords[1]&df$Log2FC<=Zoom()$coords[2]&(-log(df$padj))>=Zoom()$coords[3]&(-log(df$padj))<=Zoom()$coords[4],]
        }
        # Update the selected points according to the selection mode
        if (selectionMode() == "Brush") {
          df$selected <- ifelse(df$GeneID%in%brushedPoints(df, brushInfo()())$GeneID, "T", "F")
        } else if (selectionMode() == "Sliders"){
          df$selected <- ifelse((df$Log2FC>input$Log2FC&df$padj<input$pval)|(df$Log2FC<(-input$Log2FC)&df$padj<input$pval), "T", "F")
        }
        updateNumericInput(session, "export_GeneNumber", value=nrow(df), min=0, max=nrow(df))
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
          scale_color_manual(values=c("F" = "#384246", "T" = "#E69F00")) +
          theme(legend.position = "none") +
          xlab("log2(FoldChange)") +
          ylab("-log10(p-value ajusted)") +
          theme(text = element_text(size = 14))    
        return(plot)})
    
    pathwayDotPlot <- reactive({
      if (!is.null(pathwayEnrichment())) {
        data = pathwayEnrichment()$enrichment
        result_df <- data@result
        result_df[result_df$p.adjust <= input$qval, ]
        data@result <- result_df
        plot = dotplot(data, showCategory=30) + ggtitle("Dotplot" )  
        return(plot)
      } else {
        return(emptyPlot)
      }}
      )
    
    pathwayTreePlot <- reactive({
      if (!is.null(pathwayEnrichment())) {
        data = pathwayEnrichment()$enrichment
        dataR <- pairwise_termsim(setReadable(data, orgs[orgs$organism==input$species, "db"], 'ENTREZID'))
        if (pathwayEnrichment()$parameters$analysis=="GSEA") {
          plot = treeplot(dataR, foldChange=data@geneList)
        } else if (pathwayEnrichment()$parameters$analysis=="ORA") {
          plot = treeplot(dataR)
        }  
        return(plot)
      } else {
        return(emptyPlot)
      }
    })
    
    pathwayCnetPlot <- reactive({
      if (!is.null(pathwayEnrichment())) {
        data = pathwayEnrichment()$enrichment
        dataR = setReadable(data, orgs[orgs$organism==input$species, "db"], 'ENTREZID')
        if (pathwayEnrichment()$parameters$analysis=="GSEA") {
          plot = cnetplot(dataR, foldChange=data@geneList)
        } else {
          plot = cnetplot(dataR)
          }  
        return(plot)
      } else {
        return(emptyPlot)
      }
    })
    
    pathwayGseaPlot <- reactive({
      if ((!is.null(input$pathwayGSEA)) & (length(input$pathwayGSEA)>0)) {
        data = pathwayEnrichment()$enrichment
        plot = gseaplot2(data, as.numeric(input$pathwayGSEA)) 
        return(plot)
      } else {
        return(emptyGseaPlot)
      }
     })
      
    output$conditional_gsea_row <- renderUI({
      if (!is.null(pathwayEnrichment())) {
        if (pathwayEnrichment()$parameters$analysis=="GSEA") {
          
        }
    } else {}
      })
    
    observeEvent(input$pathwayGSEA, {
      if (!is.null(pathwayEnrichment()) & (input$pathwayGSEA[1]=="None")) {
        data = pathwayEnrichment()$enrichment
        if (pathwayEnrichment()$parameters$DB=="KEGG"){
          updateSelectInput(session ,"pathwayGSEA", choices = setNames(1:length(data$ID), data$ID))
        }
        else if (pathwayEnrichment()$parameters$DB=="Reactome"){
          updateSelectInput(session ,"pathwayGSEA", choices = setNames(1:length(data$Description), data$Description))
        }}
    })
    
    # -----------------------------------------
    # User event
    # -----------------------------------------
      
    # Reactive function controling the selection mode 
    # As shiny do not allow to disable downloadButton, disable it with shinyJS
    # Disable also other buttons with shinyJS to standardize rendering
    
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
      
      up_genes <- names(genes)[genes > 0]
      down_genes <- names(genes)[genes < 0]
      
      # Selecting the GO ontology and organism
      ontology <- input$inputGO
      go_organism <- orgs[orgs$organism==input$species, "db"]
      universe <- names(gene_list)
      library(go_organism, character.only = T)
      
      if ("Gene Set Enrichment Analysis (GSEA)" %in% input$go_analysisMethodChoice) {
        gene_list <- df$Log2FC
        names(gene_list) <- df$GeneID
        gene_list <- sort(na.omit(gene_list), decreasing = T)
        gsea_result <- gseGO(geneList = gene_list, OrgDb = get(go_organism), ont = input$inputGO, keyType = "ENSEMBL", pvalueCutoff = input$go_pval, verbose = F)
        goGseaResults(gsea_result) 
      } 
      
      if ("Over Representation Analysis (ORA)" %in% input$go_analysisMethodChoice) {
        analyze_up <- "Over expressed DEG" %in% input$go_oraChoice
        analyze_down <- "Under expressed DEG" %in% input$go_oraChoice
        
        if (!analyze_up && !analyze_down) {
          shinyalert("Error", "Select at least one type of gene expression (Over or Under expressed DEG).", type = "error")
          return(NULL)
        }
        
        if (analyze_up && length(up_genes) > 0) {
          result_up <- enrichGO(
            gene = up_genes,
            universe = universe,
            OrgDb = go_organism,
            keyType = "ENSEMBL",
            readable = T,
            ont = ontology,
            pvalueCutoff = input$go_pval,
            qvalueCutoff = input$qvalueCutoff
          )
        } else { result_up <- NULL }
        
        if (analyze_down && length(down_genes) > 0) {
          result_down <- enrichGO(
            gene = down_genes,
            universe = universe,
            OrgDb = go_organism,
            keyType = "ENSEMBL",
            readable = T,
            ont = ontology,
            pvalueCutoff = input$go_pval,
            qvalueCutoff = input$qvalueCutoff 
          )
        } else { result_down <- NULL }
        
        both_genes <- union(up_genes, down_genes)  
        if (length(both_genes) > 0) {
          result_both <- enrichGO(
            gene = both_genes,
            universe = universe,
            OrgDb = go_organism,
            keyType = "ENSEMBL",
            readable = T,
            ont = ontology,
            pvalueCutoff = input$go_pval,
            qvalueCutoff = input$qvalueCutoff 
          )
        } else { result_both <- NULL }
        
        # Storing results in goOraResults
        goOraResults(list(up = result_up, down = result_down, both = result_both))
        
      }
    })
    })
  
    #--- Pathway Enrichment ---#
    
    observe({
      if (is.na(preprocessedData()[1,1])){
        shinyalert("Bad input", "Analysis require data!", type = "error", confirmButtonCol = "#7e3535")
      }
      else if (input$analysisMethodChoice=="ORA" & is.null(input$oraChoice)){
        shinyalert("Incomplete selection!", "You should select an interest for ORA method", type = "error", confirmButtonCol = "#7e3535")
      }
      else {
      withProgress(message = "Pathway enrichment...",
        {
          message("Starting Pathway Enrichment...")
          df = as.data.frame(preprocessedData())
          pathwayEnrichment(pathway(df, organism = input$species, DB=input$dbPathwaychoice, analysis=input$analysisMethodChoice, oraInterest=input$oraChoice, pval=input$pvalPathway))},
        )}
    }) |> bindEvent(input$analysisPathwayButton)
    
    observeEvent(pathwayEnrichment(), {
      data = pathwayEnrichment()$enrichment
      if (pathwayEnrichment()$parameters$DB=="KEGG"){
        print(data)
        updateSelectInput(session ,"pathway", choices = data$ID)
        if (pathwayEnrichment()$parameters$analysis == "GSEA") {
          updateSelectInput(session ,"pathwayGSEA", choices = setNames(1:length(data$ID), data$ID))
        }
      }
      else if (pathwayEnrichment()$parameters$DB=="Reactome"){
        updateSelectInput(session ,"pathway", choices = data$Description)
        if (pathwayEnrichment()$parameters$analysis == "GSEA") {
          updateSelectInput(session ,"pathwayGSEA", choices = setNames(1:length(data$Description), data$Description))
        }
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
  
  output$goBarplotUp <- renderPlot({
    if (!is.null(goOraResults()$up)) {
      barplot(goOraResults()$up)
    } else {
      emptyPlot
    }
  })
  
  output$goBarplotDown <- renderPlot({
    if (!is.null(goOraResults()$down)) {
    #barplot(goOraResults()$down, drop = T, showCategory = input$topCategoriesDown, itle = "GO Biological Pathways (Down-regulated)", font.size = 8)
    barplot(goOraResults()$down)
    } else {
        emptyPlot
    }
  })
  
  output$goBarplotBoth <- renderPlot({
    if (!is.null(goOraResults()$both)) {
      barplot(goOraResults()$both)
      #barplot(goOraResults()$both, drop = T, showCategory = input$topCategoriesUp, title = "GO Biological Pathways (Both-regulated)", font.size = 8)
    } else {
      emptyPlot
    }
  })
  
  output$goDotplotUp <- renderPlot({
    if (!is.null(goOraResults()$up)) {
      dotplot(goOraResults()$up)
    } else {
      emptyPlot
    }
  })
  
  output$goDotplotDown <- renderPlot({
    if (!is.null(goOraResults()$down)) {
      dotplot(goOraResults()$down)
    } else {
      emptyPlot
    }
  })
  
  output$goDotplotBoth <- renderPlot({
    if (!is.null(goOraResults()$both)) {
      dotplot(goOraResults()$both)
    } else {
      emptyPlot
    }
  })
  
  # Network plots
  output$goNetplotUp <- renderPlot({
    if (!is.null(goOraResults()$up)) {
      go_enrich <- pairwise_termsim(goOraResults()$up)
      emapplot(go_enrich, layout = "kk", showCategory = 15)
      #emapplot(go_enrich, layout.params = list(layout = "kk"), showCategory = 15)
    } else {
      emptyPlot
    }
    
  })
  
  output$goNetplotDown <- renderPlot({
    if (!is.null(goOraResults()$down)) {
      go_enrich <- pairwise_termsim(goOraResults()$down)
      emapplot(go_enrich, layout = "kk", showCategory = 15)
    } else {
      emptyPlot
    }
  })
  
  output$goNetplotBoth <- renderPlot({
    if (!is.null(goOraResults()$both)) {
    go_enrich <- pairwise_termsim(goOraResults()$both)
    emapplot(go_enrich, layout = "kk", showCategory = 15)
    } else {
      emptyPlot
    }
  })
  
  # === TABLES ORA ===
  output$goTableUp <- DT::renderDT({
    if (!(is.null(goOraResults()$up))) {
      DT::datatable(as.data.frame(goOraResults()$up), options = list(pageLength = 10)) 
    } else { emptyTableGo }
  })
  
  output$goTableDown <- DT::renderDT({
    if (!(is.null(goOraResults()$down))) {
      DT::datatable(as.data.frame(goOraResults()$down), options = list(pageLength = 10)) 
    } else { emptyTableGo }
  })
  
  output$goTableBoth <- DT::renderDT({
    if (!(is.null(goOraResults()$both))) {
      DT::datatable(as.data.frame(goOraResults()$both), options = list(pageLength = 10)) 
    } else { emptyTableGo }
  })
  
  output$goGseaDotplot <- renderPlot({
    if (!(is.null(goGseaResults()))) {
      dotplot(goGseaResults())
    } else { emptyPlot }
  })
  
  output$goGseaEnrichmentPlot <- renderPlot({
    if (!(is.null(goGseaResults()))) {
      top_term <- goGseaResults()@result$ID[1]
      gseaplot(goGseaResults(), by = "all", geneSetID = top_term)
      #    pathwayGseaPlot2(goGseaResults(), geneSetID = top_term)
      } else { emptyPlot }
    })
  
  output$goGseaRidgeplot <- renderPlot({
    if (!(is.null(goGseaResults()))) {
      ridgeplot(goGseaResults(), showCategory = 13)
    } else { emptyPlot }
  })
  
  output$goGseaTable <- DT::renderDT({
    if (!(is.null(goGseaResults()))) {
      DT::datatable(as.data.frame(goGseaResults()), options = list(pageLength = 10))
    } else {
      emptyTableGo
    }
  })
  
  #--- Pathway Enrichment ---#

    output$analysisDesc <- renderUI({
      if (!is.null(pathwayEnrichment())) {
        HTML(paste(
          "<div style='font-family: Arial, sans-serif; line-height: 1.6;'>",
          "<p><strong>Organism:</strong><em>", pathwayEnrichment()$parameters$organism, "</em></p>",
          "<p><strong>Analysis:</strong> ", pathwayEnrichment()$parameters$analysis, "</p>",
          "<p><strong>Requested database:</strong> ", pathwayEnrichment()$parameters$DB, "</p>",
          "<p><strong>Log2FC threshold:</strong> ", pathwayEnrichment()$parameters$threshold, "</p>",
          "<p><strong>Ajusted p-value threshold:</strong> ", pathwayEnrichment()$parameters$pval, "</p>",
          "<p><strong>Number of enriched pathways:</strong> ", length(unique(pathwayEnrichment()$enrichment$ID)), "</p>",
          "</div>"
        ))
      } else {
        HTML("<i>You must launch a pathway enrichment to explore its results.</i>")
      }
    })
    
    output$pathwaytable <- DT::renderDT ({
      message("Rendering Pathway DataTable")
      if (!is.null(pathwayEnrichment())){
        df = as.data.frame(pathwayEnrichment()$enrichment)
        df = df[df$p.adjust <= input$qval,] 
        for (i in 1:length(df$ID)) {
          if (pathwayEnrichment()$parameters$DB=="KEGG"){
            df$ID[i] = paste('<a href="https://www.genome.jp/pathway/', df$ID[i], '" target="_blank">', df$ID[i], '</a>', sep="")
          } else if (pathwayEnrichment()$parameters$DB=="Reactome"){
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
    
    output$pathwayplotout <- renderPlot ({
      message("Rendering Pathway Dot plot")
      if (!is.null(pathwayDotPlot())){
        pathwayDotPlot()}
      else {
      }
    })
      
    output$pathwayplotout2 <- renderPlot ({
      message("Rendering Pathway Tree plot")
      if (!is.null(pathwayTreePlot())){
        pathwayTreePlot()}
      else {
      }
    })  
    
    output$pathwayplotout3 <- renderPlot ({
      message("Rendering Pathway Gene-Concept Network plot")
      if (!is.null(pathwayCnetPlot())){
        pathwayCnetPlot()}
      else {
      }
    })  
    
    output$pathwayGseaPlot <- renderPlot ({
      message("Rendering Pathway GSEA plot")
      if (!is.null(pathwayGseaPlot()) & pathwayEnrichment()$parameters$analysis=="GSEA"){
        pathwayGseaPlot()}
      else {
      }
  })
    
    output$pathway <- renderUI({
      if (!is.null(pathwayEnrichment())) {
        if (pathwayEnrichment()$parameters$DB == "KEGG"){
          imageOutput("pathwayKegg", height = "750px")
        } else {
          imageOutput("pathwayReactome", height = "750px")
        }
      } else {
          plotOutput("pathwayEmpty")
      }
    }) |> bindEvent(input$pathway)
    
    output$pathwayEmpty <- renderPlot ({emptyPlot})
    
    output$pathwayKegg <- renderImage(deleteFile=F, {
      message("Rendering Pathway plot")
      if (input$pathway!="None" & pathwayEnrichment()$parameters$DB == "KEGG"){
        list(
          src = paste(getwd(), "/out/", input$pathway, ".png", sep=""),
          contentType = 'image/png',
          alt = "Image PNG",
          height="750Px"
        )} else {
        }
  })
    
    output$pathwayReactome <- renderImage(deleteFile=F, {
     # message("Warning : plotting fold change only for non-duplicate gene")
     if (input$pathway!="None" & pathwayEnrichment()$parameters$DB == "Reactome"){
       list(
         src = paste(getwd(), "/out/", gsub(" ", "_", input$pathway), ".png", sep=""),
         contentType = 'image/png',
         alt = "Image PNG",
         height="750Px"
       )} else {
       }}
      )
    
  # -----------------------------------------
  # Save and export functions
  # -----------------------------------------

   output$export_options <- renderUI({
      choices <- list()
      if (!is.null(preprocessedData()) && !is.na(preprocessedData()[1,1])) {
        choices[["Whole Data Inspection (current VolcanoPlot & Table)"]] <- "WDI"
      }
      if (!is.null(goOraResults())) {
        res <- goOraResults()
        if (!is.null(res$up) || !is.null(res$down) || !is.null(res$both)) {
          choices[["GO enrichment"]] <- "GO"
        }
      }
      if (!is.null(pathwayEnrichment())) { choices[["Pathway enrichment"]] <- "PATH" }
      checkboxGroupInput("export_choices",
                         "Select results to export:",
                         choices = choices)
    })
    
  output$exportOptions <- downloadHandler(
    filename = function() {
      paste0("HEATraN_results_", Sys.Date(), ".html")
    },
    content = function(file) {
      tempReport <- file.path(tempdir(), "template.Rmd")
      file.copy("www/template.Rmd", tempReport, overwrite = T)
  
      pathwayViewPlots <- list()
      
      if (!is.null(pathwayEnrichment())) {
        data <- pathwayEnrichment()$enrichment
        db <- pathwayEnrichment()$parameters$DB
        
        if (db == "KEGG") {
          for (i in 1:nrow(data)) {
            tryCatch({
              pathway_id <- data$ID[i]
              img_path <- paste0("./out/", pathway_id, ".pathview.png")
              if (file.exists(img_path)) {
                pathwayViewPlots[[i]] <- png::readPNG(img_path)
              } else {
                pathwayViewPlots[[i]] <- NULL
              }
            }, error = function(e) {
              message(paste("Error loading KEGG pathway image:", e$message))
              pathwayViewPlots[[i]] <- NULL
            })
          }
        } else if (db == "Reactome") {
          for (i in 1:nrow(data)) {
            tryCatch({
              pathway_name <- data$Description[i]
              safe_name <- gsub("[^a-zA-Z0-9]", "_", pathway_name)
              img_path <- paste0("./out/", safe_name, ".png")
              
              if (file.exists(img_path)) {
                pathwayViewPlots[[i]] <- png::readPNG(img_path)
              } else if (!is.null(pathwayEnrichment()$pathway_images) && !is.null(pathwayEnrichment()$pathway_images[[i]])) {
                img_path <- pathwayEnrichment()$pathway_images[[i]]
                if (file.exists(img_path)) {
                  pathwayViewPlots[[i]] <- png::readPNG(img_path)
                }
              } else {
                pathway_id <- data$ID[i]
                api_url <- paste0("https://reactome.org/ContentService/exporter/diagram/", 
                                  pathway_id, ".png?quality=7")
                
                tmp_file <- tempfile(fileext = ".png")
                download.file(api_url, tmp_file, mode = "wb", quiet = T)
                if (file.exists(tmp_file)) {
                  pathwayViewPlots[[i]] <- png::readPNG(tmp_file)
                  file.remove(tmp_file)
                } else {
                  pathwayViewPlots[[i]] <- NULL
                }
              }
            }, error = function(e) {
              message(paste("Error loading Reactome pathway image:", e$message))
              pathwayViewPlots[[i]] <- NULL
            })
          }
        }
      }
      
      includeWDI = ifelse("WDI" %in% input$exportSections, T, F)
      includeGO = ifelse("GO" %in% input$exportSections, T, F)
      includePATHWAY = ifelse("PATHWAY" %in% input$exportSections, T, F)
      
      # Rendre le document R Markdown
      rmarkdown::render(
        input = tempReport,
        output_file = file,
        params = list(
          nGenes = input$nGenesExport,
          volcanoPlot = plot(),
          tableData = processedData()[processedData()$selected == T, -c("selected", "minuslog10")],
          enrichmentData = pathwayEnrichment(),
          pathwayGseaPlot = if (!is.null(pathwayEnrichment()) && pathwayEnrichment()$parameters$analysis == "GSEA") {
            lapply(1:nrow(pathwayEnrichment()$enrichment), function(i) {
              pathwayGseaPlot2(pathwayEnrichment()$enrichment, i)
            })
          } else NULL,
          pathwayTreePlot = if (!is.null(pathwayEnrichment())) pathwayTreePlot() else NULL,
          pathwayDotPlot = if (!is.null(pathwayEnrichment())) pathwayDotPlot() else NULL,
          pathwayCnetPlot = if (!is.null(pathwayEnrichment())) pathwayCnetPlot() else NULL,
          pathwayViewPlots = pathwayViewPlots,
          organismInfo = input$species,
          includeWDI = includeWDI,
          includeGO = includeGO,
          includePATHWAY = includePATHWAY,
          includepathwayGseaPlots = input$includeGSEA,
          includePathwayViews = input$includePathwayViews
        ),
        envir = new.env()
        )
      }
    )

  
  
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

}

