# Developed by LESAGE Louison (@thelokj).
# louison.lesage@univ-rouen.fr
# Student at Rouen Normandy University
# University project 2024-2025
# Last updated : 18/11/2024
# HEATraN version 0.2.0-a.5

setwd("../")
source("./bin/fun/pathway.R")

emptyTable <- data.frame(Gene=NA, Log2FC=NA, p_value=NA)
emptyTable2 <- data.frame(Pathway=NA, p_value=NA, q_value=NA)
brushInfo <- reactiveVal(NULL)
enrichment <- reactiveVal(NULL)
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
    
    dotPlot <- reactive({
      data = enrichment()$enrichment
      plot = dotplot(data, showCategory=30) + ggtitle("DotPlot" )   
      return(plot)})
    
    barPlot <- reactive({
      data = enrichment()$enrichment
      plot = barplot(data, showCategory=30) + ggtitle("BarPlot" )   
      return(plot)})
    
    gseaPlot <- reactive({
      data = enrichment()$enrichment
      plot = gseaplot(data, 1) 
      return(plot)})
    
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
    
    ### Pathway
    observeEvent(enrichment(), {
      data = enrichment()$enrichment
      if (enrichment()$parameters$DB=="KEGG"){
      updateSelectInput(session ,"pathway",
                      choices = data$ID)}
      else if (enrichment()$parameters$DB=="Reactome"){
        updateSelectInput(session ,"pathway",
                          choices = data$Description)}
    }) 
    
    output$analysisDesc <- renderUI({
      if (!is.null(enrichment())){
        HTML(paste("<b>Organism :</b>", enrichment()$parameters$organism,
              "| <b>Requested database :</b>", enrichment()$parameters$DB,
              "| <b>Analysis :</b>", enrichment()$parameters$analysis,
              "| <b>Log2FC threshold :</b>", enrichment()$parameters$threshold,
              "| <b>Pval :</b>", enrichment()$parameters$pval,
              "<br><b> No. of enriched pathways :</b>", length(unique(enrichment()$enrichment$ID)))
             )
      } else {
        HTML("<i>You must launch a pathway enrichment to explore its results.</i>")}
    })
    
    output$pathwaytable <- DT::renderDT ({
      message("Rendering Pathway DataTable")
      if (!is.null(enrichment())){
        df = as.data.frame(enrichment()$enrichment)
        for (i in 1:length(df$ID)) {
          if (enrichment()$parameters$DB=="KEGG"){
            df$ID[i] = paste('<a href="https://www.genome.jp/pathway/', df$ID[i], '" target="_blank">', df$ID[i], '</a>', sep="")
          } else if (enrichment()$parameters$DB=="Reactome"){
            df$ID[i] = paste('<a href="https://reactome.org/PathwayBrowser/#/', df$ID[i], '" target="_blank">', df$ID[i], '</a>', sep="")
          }
        }
        datatable(df, escape = FALSE)}
      else {
        emptyTable2
      }
    })
    
    output$pathwayplotout <- renderPlot ({
      message("Rendering Pathway dotPlot")
      if (!is.null(dotPlot())){
        dotPlot()}
      else {
      }
    })
      
    output$pathwayplotout2 <- renderPlot ({
      message("Rendering Pathway batPlot")
      if (!is.null(barPlot())){
        barPlot()}
      else {
      }
    })  
    
    output$gseaplot <- renderPlot ({
      message("Rendering Pathway GSEA plot")
      if (!is.null(gseaPlot())){
        gseaPlot()}
      else {
      }
  })
    
    output$pathway <- renderUI({
      if (enrichment()$parameters$DB == "KEGG"){
        imageOutput("pathwayKegg", width = "100%", height = "auto")
      } else {
        imageOutput("pathwayReactome", width = "100%", height = "auto")
      }
    }) |> bindEvent(input$pathway)
    
    output$pathwayKegg <- renderImage(deleteFile=F, {
      message("Rendering Pathway plot")
      if (input$pathway!="None" & enrichment()$parameters$DB == "KEGG"){
        list(
          src = paste(getwd(), "/out/", input$pathway, ".png", sep=""),
          contentType = 'image/png',
          alt = "Image PNG"
        )} else {
        }
  })
    
    output$pathwayReactome <- renderImage(deleteFile=F, {
     # message("Warning : plotting fold change only for non-duplicate gene")
     if (input$pathway!="None" & enrichment()$parameters$DB == "Reactome"){
       list(
         src = paste(getwd(), "/out/", gsub(" ", "_", input$pathway), ".png", sep=""),
         contentType = 'image/png',
         alt = "Image PNG"
       )} else {
       }}
      )
}
