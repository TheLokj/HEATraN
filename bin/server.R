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
emptyTable2 <- data.frame(Pathway=NA, p_value=NA, q_value=NA)
brushInfo <- reactiveVal(NULL)
enrichment <- reactiveVal(NULL)
selectionMode <- reactiveVal("None")
preprocessedData <- reactiveVal(emptyTable)
requiredNames <- c(config$DATA$gene_name, config$DATA$gene_id, config$DATA$basemean, config$DATA$log2FC,  config$DATA$pvalue, config$DATA$adjusted_pvalue)
emptyPlot <- ggplot() + 
  geom_richtext(aes(x = 0.5, y = 0.5, 
                    label = "<b style='font-family:Mulish,sans-serif;font-size:20pt;'>Data required for visual representation</b>"), 
                size = 6, color = "red3", fill = NA, label.color = NA) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = NA))

emptyGSEAPlot <- ggplot() + 
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
    
    dotPlot <- reactive({
      if (!is.null(enrichment())) {
        data = enrichment()$enrichment
        result_df <- data@result
        result_df[result_df$p.adjust <= input$qval, ]
        data@result <- result_df
        plot = dotplot(data, showCategory=30) + ggtitle("DotPlot" )  
        return(plot)
      } else {
        return(emptyPlot)
      }}
      )
    
    treePlot <- reactive({
      if (!is.null(enrichment())) {
        data = enrichment()$enrichment
        dataR <- pairwise_termsim(setReadable(data, orgs[orgs$organism==input$species, "db"], 'ENTREZID'))
        if (enrichment()$parameters$analysis=="GSEA") {
          plot = treeplot(dataR, foldChange=data@geneList)
        } else if (enrichment()$parameters$analysis=="ORA") {
          plot = treeplot(dataR)
        }  
        return(plot)
      } else {
        return(emptyPlot)
      }
      })
    
    cnetPlot <- reactive({
      if (!is.null(enrichment())) {
        data = enrichment()$enrichment
        dataR = setReadable(data, orgs[orgs$organism==input$species, "db"], 'ENTREZID')
        if (enrichment()$parameters$analysis=="GSEA") {
          plot = cnetplot(dataR, foldChange=data@geneList)
        } else {
          plot = cnetplot(dataR)
          }  
        return(plot)
      } else {
        return(emptyPlot)
      }
      })
    
    #### Change here to select the GSEA Pathway according to ID 
    gseaPlot <- reactive({
      if ((!is.null(input$pathwayGSEA)) & (length(input$pathwayGSEA)>0)) {
        data = enrichment()$enrichment
        plot = gseaplot2(data, as.numeric(input$pathwayGSEA)) 
        return(plot)
      } else {
        return(emptyGSEAPlot)
      }
      })
      
    output$conditional_gsea_row <- renderUI({
      if (!is.null(enrichment())) {
        if (enrichment()$parameters$analysis=="GSEA") {
          fluidRow(
            box(
              title = HTML("<b>GSEA plot</b>"),
              id = "gseaplot", width = 12,
              
              selectInput("pathwayGSEA", "Select one or more pathways:", choices = c("None"), selected = "None", multiple=TRUE),
              plotOutput("gseaplot",
                         height = "425px",
                         click = "plot_click",
                         dblclick = "plot_dblclick",
                         hover = "plot_hover",
                         brush = brushOpts(id = "plot_brush", delay = 3000, delayType = "debounce", fill="#7e3535", stroke="#7e3535"))
            )
          )
        }
    } else {}
      })
    
    observeEvent(input$pathwayGSEA, {
      if (!is.null(enrichment()) & (input$pathwayGSEA[1]=="None")) {
        data = enrichment()$enrichment
        if (enrichment()$parameters$DB=="KEGG"){
          updateSelectInput(session ,"pathwayGSEA", choices = setNames(1:length(data$ID), data$ID))
        }
        else if (enrichment()$parameters$DB=="Reactome"){
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
          enrichment(pathway(df, organism = input$species, DB=input$dbPathwaychoice, analysis=input$analysisMethodChoice, oraInterest=input$oraChoice, pval=input$pvalPathway))},
        )}
    }) |> bindEvent(input$analysisPathwayButton)
  
    # -----------------------------------------
    # UI output 
    # -----------------------------------------

    # Rendering "Current selection" section 
    output$InfoSelect <- renderUI({
      if (is.null(preprocessedData) || is.na(preprocessedData()[1,1])){
        HTML(paste("<b style='color:	#FF0000'>Data required for selection</b><br/><br/>", sep=""))
      }
      else {
        HTML(paste("<b>Current selection</b><br/><i>p-value</i>: [", format(exp(-input$plot_brush$ymax), scientific=TRUE, digits=3), " ; " , format(exp(-input$plot_brush$ymin), scientific=TRUE, digits=3),"]", "     <br/>    ", "<i>log2(FoldChange)</i>: ", "[", round(input$plot_brush$xmin, 3), " ; ", round(input$plot_brush$xmax, 3),"]<br/><br/>", sep=""))
      }}
  )
  
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
  goResults <- reactiveVal(NULL)
  # Reactive value for GSEA results
  gseaResults <- reactiveVal(NULL)
  
  observeEvent(input$go_analysisButton, {
    req(preprocessedData())
    
    withProgress(message = 'Running GO analysis...', {
      
      df <- preprocessedData()
      
      # Liste des gènes pour l'univers
      original_gene_list <- df$Log2FC
      names(original_gene_list) <- df$GeneID
      gene_list <- sort(na.omit(original_gene_list), decreasing = TRUE)
      
      # Sélection des gènes selon le seuil de p-value ajustée
      sig_genes_df <- df[df$pval < input$go_pval,]
      genes <- sig_genes_df$Log2FC
      names(genes) <- sig_genes_df$GeneID
      genes <- na.omit(genes)
      
      up_genes <- names(genes)[genes > 0]
      down_genes <- names(genes)[genes < 0]
      
      # Sélection de l'ontologie GO et de l'organisme
      ontology <- input$inputGO
      go_organism <- orgs[orgs$organism==input$species, "db"]
      universe <- names(gene_list)
      library(go_organism, character.only = TRUE)
      
      # Méthode d’analyse : ORA ou GSEA
      if ("Gene Set Enrichment Analysis (GSEA)" %in% input$go_analysisMethodChoice) {
        df <- preprocessedData()
        gene_list <- df$Log2FC
        names(gene_list) <- df$GeneID
        gene_list <- sort(na.omit(gene_list), decreasing = TRUE)
        
        #analyse GSEA
        gsea_result <- gseGO(
          geneList      = gene_list,
          OrgDb         = get(input$organism),  
          ont           = input$inputGO,        
          keyType       = "ENSEMBL",
          pvalueCutoff  = input$go_pval,
          verbose       = FALSE
        )
        
     
        
        gseaResults(gsea_result) 
       
 
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
            OrgDb = go_organism,
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
            OrgDb = go_organism,
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
            OrgDb = go_organism,
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
  
  output$goBarplotUp <- renderPlot({
    req(goResults()$up)
    barplot(goResults()$up)
  })
  
  output$goBarplotDown <- renderPlot({
    req(goResults()$down)
    barplot(goResults()$down,
            drop = TRUE,
            showCategory = input$topCategoriesDown,
            title = "GO Biological Pathways (Down-regulated)",
            font.size = 8)
  })
  
  output$goBarplotBoth <- renderPlot({
    req(goResults()$both)
    barplot(goResults()$both,
            drop = TRUE,
            showCategory = input$topCategoriesUp,  # tu peux ajuster si nécessaire
            title = "GO Biological Pathways (Both-regulated)",
            font.size = 8)
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
    DT::datatable(as.data.frame(goResults()$up), options = list(pageLength = 10))
  })
  
  output$goTableDown <- DT::renderDT({
    req(goResults()$down)
    DT::datatable(as.data.frame(goResults()$down), options = list(pageLength = 10))
  })
  
  output$goTableBoth <- DT::renderDT({
    req(goResults()$both)
    DT::datatable(as.data.frame(goResults()$both), options = list(pageLength = 10))
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
        return(emptyPlot)
      }
    })
    
    ### Pathway
    observeEvent(enrichment(), {
      data = enrichment()$enrichment
      if (enrichment()$parameters$DB=="KEGG"){
        print(data)
        updateSelectInput(session ,"pathway", choices = data$ID)
        if (enrichment()$parameters$analysis == "GSEA") {
        updateSelectInput(session ,"pathwayGSEA", choices = setNames(1:length(data$ID), data$ID))
          }
      }
      else if (enrichment()$parameters$DB=="Reactome"){
        updateSelectInput(session ,"pathway", choices = data$Description)
        if (enrichment()$parameters$analysis == "GSEA") {
        updateSelectInput(session ,"pathwayGSEA", choices = setNames(1:length(data$Description), data$Description))
        }
      }
    }) 
    
    output$analysisDesc <- renderUI({
      if (!is.null(enrichment())) {
        HTML(paste(
          "<div style='font-family: Arial, sans-serif; line-height: 1.6;'>",
          "<p><strong>Organism:</strong><em>", enrichment()$parameters$organism, "</em></p>",
          "<p><strong>Analysis:</strong> ", enrichment()$parameters$analysis, "</p>",
          "<p><strong>Requested database:</strong> ", enrichment()$parameters$DB, "</p>",
          "<p><strong>Log2FC threshold:</strong> ", enrichment()$parameters$threshold, "</p>",
          "<p><strong>Ajusted p-value threshold:</strong> ", enrichment()$parameters$pval, "</p>",
          "<p><strong>Number of enriched pathways:</strong> ", length(unique(enrichment()$enrichment$ID)), "</p>",
          "</div>"
        ))
      } else {
        HTML("<i>You must launch a pathway enrichment to explore its results.</i>")
      }
    })
    
    output$pathwaytable <- DT::renderDT ({
      message("Rendering Pathway DataTable")
      if (!is.null(enrichment())){
        df = as.data.frame(enrichment()$enrichment)
        df = df[df$p.adjust <= input$qval,] 
        for (i in 1:length(df$ID)) {
          if (enrichment()$parameters$DB=="KEGG"){
            df$ID[i] = paste('<a href="https://www.genome.jp/pathway/', df$ID[i], '" target="_blank">', df$ID[i], '</a>', sep="")
          } else if (enrichment()$parameters$DB=="Reactome"){
            df$ID[i] = paste('<a href="https://reactome.org/PathwayBrowser/#/', df$ID[i], '" target="_blank">', df$ID[i], '</a>', sep="")
          }
        }
        if (nrow(df) != 0) {
          datatable(df, escape = FALSE)
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
      if (!is.null(dotPlot())){
        dotPlot()}
      else {
      }
    })
      
    output$pathwayplotout2 <- renderPlot ({
      message("Rendering Pathway Tree plot")
      if (!is.null(treePlot())){
        treePlot()}
      else {
      }
    })  
    
    output$pathwayplotout3 <- renderPlot ({
      message("Rendering Pathway Gene-Concept Network plot")
      if (!is.null(cnetPlot())){
        cnetPlot()}
      else {
      }
    })  
    
    output$gseaplot <- renderPlot ({
      message("Rendering Pathway GSEA plot")
      if (!is.null(gseaPlot()) & enrichment()$parameters$analysis=="GSEA"){
        gseaPlot()}
      else {
      }
  })
    
    output$pathway <- renderUI({
      if (!is.null(enrichment())) {
        if (enrichment()$parameters$DB == "KEGG"){
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
      if (input$pathway!="None" & enrichment()$parameters$DB == "KEGG"){
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
     if (input$pathway!="None" & enrichment()$parameters$DB == "Reactome"){
       list(
         src = paste(getwd(), "/out/", gsub(" ", "_", input$pathway), ".png", sep=""),
         contentType = 'image/png',
         alt = "Image PNG",
         height="750Px"
       )} else {
       }}
      )
    
    onStop(function() {
      if (config$FILE$clear_cache == "True") {
        if ((dir_size("./out") / (1024^2)) > config$FILE$max_cache_mb) {
          message("Cache exceeding the authorised limit: cleaning cache...")
          files_to_delete <- list.files(path = "./out", full.names = TRUE)
          unlink(files_to_delete)
        }
      }
      message("Thanks for using HEATraN!")
    })
      
# Téléchargement du rapport complet
output$downloadFullReport <- downloadHandler(
  {print(input$export_ListOptions)
  if(!is.null(enrichment()) && enrichment()$parameters$analysis=="GSEA" && "pathways"%in%input$export_ListOptions) {
    gseaPlotList = list()
    for (i in 1:length(enrichment()$enrichment$ID)) {
      gseaPlotList[[i]] = gseaplot2(data, i)
    } 
  } else {
    gseaPlotList = NULL
  }},
  
  filename = function() {
    paste0("HEATraN_Rapport_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".HTML")
  },
  content = function(file) {
    tempReport <- file.path(tempdir(), "rapport.Rmd")
    file.copy("www/template.Rmd", tempReport, overwrite = TRUE)
    params <- list(
      nGenes = input$export_GeneNumber,
      volcanoPlot = if(!is.na(preprocessedData()[1,1])) { plot() } else { NULL },
      tableData = if(!is.na(preprocessedData()[1,1])) { 
        processedData()[processedData()$selected==TRUE,-c("selected","minuslog10")] 
      } else { NULL },
      enrichmentData = enrichment(),
      dotPlot = if(!identical(dotPlot, emptyPlot)){ dotPlot() } else { NULL },
      treePlot = if(!identical(treePlot, emptyPlot)){ treePlot() } else { NULL },
      cnetPlot = if(!identical(cnetPlot, emptyPlot)){ cnetPlot() } else { NULL },
      gseaPlot = gseaPlotList,
      organismInfo = input$species
    )
    
    # Générer le rapport
    rmarkdown::render(tempReport, 
                      output_file = file,
                      output_format = paste0("html_document"),
                      params = params,
                      envir = new.env(parent = globalenv()))
  }
)


# Fonction pour exporter les résultats d'analyse
output$exportReport <- downloadHandler(
  filename = function() {
    paste0("HEATraN_results_", Sys.Date(), ".html")
  },
  # GSEA Dotplot
  output$gseaDotplot <- renderPlot({
    req(gseaResults())
    dotplot(gseaResults())
  })
  
  # GSEA Enrichment plot (top term)
  output$gseaEnrichmentPlot <- renderPlot({
   req(gseaResults())
   top_term <- gseaResults()@result$ID[1]
   gseaplot(gseaResults(), by = "all", geneSetID = top_term)
#    gseaplot2(gseaResults(), geneSetID = top_term)
  }) 
  #GSEA RIDGE PLOT
  # GSEA Ridgeplot
  output$gseaRidgeplot <- renderPlot({
    req(gseaResults())
    ridgeplot(gseaResults(), showCategory = 13)
  })
  
  # GSEA TABLE
  output$gseaTable <- DT::renderDT({
    req(gseaResults())
    DT::datatable(as.data.frame(gseaResults()), options = list(pageLength = 10))
  })
  
  
  content = function(file) {
    tempReport <- file.path(tempdir(), "template.Rmd")
    file.copy("www/template.Rmd", tempReport, overwrite = TRUE)
    
    # Création d'une liste pour stocker les graphiques de visualisation des pathways
    pathwayViewPlots <- list()
    
    if (!is.null(enrichment())) {
      data <- enrichment()$enrichment
      db <- enrichment()$parameters$DB
      
      if (db == "KEGG") {
        # Pour KEGG, utiliser les images générées par pathview
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
        # Pour Reactome, utiliser les images sauvegardées ou les récupérer si nécessaire
        for (i in 1:nrow(data)) {
          tryCatch({
            pathway_name <- data$Description[i]
            safe_name <- gsub("[^a-zA-Z0-9]", "_", pathway_name)
            img_path <- paste0("./out/", safe_name, ".png")
            
            if (file.exists(img_path)) {
              pathwayViewPlots[[i]] <- png::readPNG(img_path)
            } else if (!is.null(enrichment()$pathway_images) && !is.null(enrichment()$pathway_images[[i]])) {
              # Utiliser le chemin stocké dans l'objet enrichment s'il existe
              img_path <- enrichment()$pathway_images[[i]]
              if (file.exists(img_path)) {
                pathwayViewPlots[[i]] <- png::readPNG(img_path)
              }
            } else {
              # Si l'image n'existe pas encore, générer à la volée
              pathway_id <- data$ID[i]
              api_url <- paste0("https://reactome.org/ContentService/exporter/diagram/", 
                                pathway_id, ".png?quality=7")
              
              # Télécharger l'image
              tmp_file <- tempfile(fileext = ".png")
              download.file(api_url, tmp_file, mode = "wb", quiet = TRUE)
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
    
    includeWDI = ifelse("WDI" %in% input$exportSections, TRUE, FALSE)
    includeGO = ifelse("GO" %in% input$exportSections, TRUE, FALSE)
    includePATHWAY = ifelse("PATHWAY" %in% input$exportSections, TRUE, FALSE)
    
    # Rendre le document R Markdown
    rmarkdown::render(
      input = tempReport,
      output_file = file,
      params = list(
        nGenes = input$nGenesExport,
        volcanoPlot = plot(),
        tableData = processedData()[processedData()$selected == TRUE, -c("selected", "minuslog10")],
        enrichmentData = enrichment(),
        gseaPlot = if (!is.null(enrichment()) && enrichment()$parameters$analysis == "GSEA") {
          lapply(1:nrow(enrichment()$enrichment), function(i) {
            gseaplot2(enrichment()$enrichment, i)
          })
        } else NULL,
        treePlot = if (!is.null(enrichment())) treePlot() else NULL,
        dotPlot = if (!is.null(enrichment())) dotPlot() else NULL,
        cnetPlot = if (!is.null(enrichment())) cnetPlot() else NULL,
        pathwayViewPlots = pathwayViewPlots,
        organismInfo = input$species,
        includeWDI = includeWDI,
        includeGO = includeGO,
        includePATHWAY = includePATHWAY,
        includeGSEAplots = input$includeGSEA,
        includePathwayViews = input$includePathwayViews
      ),
      envir = new.env()
    )
  }
)
}
