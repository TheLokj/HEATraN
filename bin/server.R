# Developed by DAHER Rayan, NAIT EL DJOUDI Lamia & LESAGE Louison (@thelokj).
# rayan.daher@univ-rouen.fr
# lamia.nait-el-djoudi@univ-rouen.fr
# louison.lesage@univ-rouen.fr
# Students at Rouen Normandy University
# Master of Bioinformatics, class M2.2 BIMS 2026 
# Last updated : 04/06/2025
# HEATraN version 1.0.0

config <- read.ini("./conf.ini")

app_palette = colorRampPalette(c("#ef940b", "#7e3535"))
options(enrichplot.colours = c("#ef940b", "#7e3535"))

empty_table <- data.frame(Gene=NA, Log2FC=NA, p_value=NA)
empty_go_table <- data.frame(GO=NA, Description=NA, p_value=NA, q_value=NA)
empty_pathway_table <- data.frame(Pathway=NA, p_value=NA, q_value=NA)
brush_info <- reactiveVal(NULL)
sig_go <- reactiveVal(NULL)
sig_go_result <- reactiveVal(NULL)
pathway_ora_enrichment <- reactiveVal(NULL)
pathway_gsea_enrichment <- reactiveVal(NULL)
go_ora_results <- reactiveVal(NULL)
go_gsea_results <- reactiveVal(NULL)
selection_mode <- reactiveVal("None")
preprocessed_data <- reactiveVal(empty_table)
required_names <- c(config$DATA$gene_name, config$DATA$gene_id, config$DATA$basemean, config$DATA$log2FC,  config$DATA$pvalue, config$DATA$adjusted_pvalue)

empty_plot <- ggplot() + 
  geom_richtext(aes(x = 0.5, y = 0.5, 
                    label = "<b style='font-family:Mulish,sans-serif;font-size:20pt;'>Data required for visual representation</b>"), 
                size = 6, color = "red3", fill = NA, label.color = NA) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = NA))

empty_gsea_plot <- ggplot() + 
  geom_richtext(aes(x = 0.5, y = 0.5, 
                    label = "<b style='font-family:Mulish,sans-serif;font-size:20pt;'>Please select a pathway</b>"), 
                size = 3, color = "red3", fill = NA, label.color = NA) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = NA))

empty_pathway_plot <- ggplot() + 
  geom_richtext(aes(x = 0.5, y = 0.5, 
                    label = "<b style='font-family:Mulish,sans-serif;font-size:20pt;'>Please select an enriched data</b>"), 
                size = 6, color = "red3", fill = NA, label.color = NA) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = NA))

not_enough_data_plot <- ggplot() + 
  geom_richtext(aes(x = 0.5, y = 0.5, 
                    label = "<b style='font-family:Mulish,sans-serif;font-size:20pt;'>Insufficient enriched data for clustering</b>"), 
                size = 6, color = "red3", fill = NA, label.color = NA) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = NA))

server <- function(input, output, session) {
  
  hideTab(inputId = "go_tabs", target = "Results (ORA)", session = session)
  hideTab(inputId = "go_tabs", target = "Results (GSEA)", session = session)
  
  hideTab(inputId = "pathway_tabs", target = "Results (ORA)", session = session)
  hideTab(inputId = "pathway_tabs", target = "Results (GSEA)", session = session)
  hideTab(inputId = "pathway_tabs", target = "View Pathway", session = session)
  
  # -----------------------------------------
  # Data processing
  # -----------------------------------------
  
  # Reactive function controlling the imported data
  imported_data <- reactive ({
    message("Importing data")
    # While there is no input, return NULL
    if(is.null(input$table_input)){
      return(NULL)}
    # After the importation, print the loaded table
    else{
      input_path <- (input$table_input)[1,4]
      format <- tail(strsplit(input_path, ".", fixed=T)[[1]], 1)
      # Check the file format
      if(format=="csv" || format=="tsv"){
        table = fread(input_path, sep="auto", h=T)
      } else if(format=="xlsx" | format=="xls"){
        # Inform the user for the excel files
        shinyalert("Excel file imported", "The first datasheet of this file will be imported.", confirmButtonCol = "#7e3535")
        table <- read_excel(input_path, sheet=1)
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
  }) |> bindEvent(input$table_input, ignoreNULL=F, ignoreInit=T)
  
  # Observe importation in order to preprocess data
  observe ({
    message("Preprocessing data")
    data_to_preprocess <- imported_data()
    # If the required columns aren't found, 
    if (FALSE %in% c(required_names %in% colnames(data_to_preprocess))){
      shinyalert(html = TRUE, "Please select the variables", "Wrong column name", type = "info", confirmButtonCol = "#7e3535",
                 text = tagList(
                   selectInput("gene_name_col", "Gene Name", colnames(data_to_preprocess)),
                   selectInput("gene_id_col", "Gene ID", colnames(data_to_preprocess)),
                   selectInput("basemean_col", "Base mean", colnames(data_to_preprocess)),
                   selectInput("log2fc_col", "Log2(FoldChange)", colnames(data_to_preprocess)),
                   selectInput("pval_col", "p-value", colnames(data_to_preprocess)),
                   selectInput("padj_col", "Adjusted p-value", colnames(data_to_preprocess)),
                   checkboxInput("save_colnames", "use these column names in the future"),
                   HTML(paste("Next time, use these names to import direcly the table : <i>", paste(required_names, collapse=", "), "</i><br>")),
                 ),
                 callbackR = function(finished) { 
                   if(finished) {
                     if (length(unique(c(input$gene_name_col, input$gene_id_col, input$basemean_col,input$log2fc_col,input$pval_col, input$padj_col)))==6){
                       if (input$save_colnames==T) { 
                         config$DATA$gene_name = input$gene_name_col
                         config$DATA$gene_id = input$gene_id_col
                         config$DATA$basemean = input$basemean_col
                         config$DATA$log2FC = input$log2fc_col
                         config$DATA$pvalue = input$pval_col
                         config$DATA$adjusted_pvalue = input$padj_col
                         write.ini(config, "./conf.ini")
                       }
                       data_to_preprocess <- imported_data()
                       colnames(data_to_preprocess)[colnames(data_to_preprocess)==input$gene_name_col] = "GeneName"
                       colnames(data_to_preprocess)[colnames(data_to_preprocess)==input$gene_id_col] = "GeneID"
                       colnames(data_to_preprocess)[colnames(data_to_preprocess)==input$basemean_col] = "baseMean"
                       colnames(data_to_preprocess)[colnames(data_to_preprocess)==input$log2fc_col] = "Log2FC"
                       colnames(data_to_preprocess)[colnames(data_to_preprocess)==input$pval_col] = "pval"
                       colnames(data_to_preprocess)[colnames(data_to_preprocess)==input$padj_col] = "padj"
                       data_to_preprocess$minuslog10 <- -log(data_to_preprocess$pval)
                       preprocessed_data(data_to_preprocess)
                       selection_mode("Sliders")
                     } else {
                       shinyalert("Incorrect choice!", "Each column must be unique.", type = "error", confirmButtonCol = "#7e3535")
                     }
                   } else {return(FALSE)}})
    } else {
      colnames(data_to_preprocess)[colnames(data_to_preprocess)==config$DATA$gene_name] = "GeneName"
      colnames(data_to_preprocess)[colnames(data_to_preprocess)==config$DATA$gene_id] = "GeneID"
      colnames(data_to_preprocess)[colnames(data_to_preprocess)==config$DATA$base_mean] = "baseMean"
      colnames(data_to_preprocess)[colnames(data_to_preprocess)==config$DATA$log2FC] = "Log2FC"
      colnames(data_to_preprocess)[colnames(data_to_preprocess)==config$DATA$pvalue] = "pval"
      colnames(data_to_preprocess)[colnames(data_to_preprocess)==config$DATA$adjusted_pvalue] = "padj"
      data_to_preprocess$minuslog10 <- -log(data_to_preprocess$pval)
      preprocessed_data(data_to_preprocess)
      selection_mode("Sliders")}}) |> bindEvent(imported_data())
  
  # Reactive function containing the selected points
  processed_data <- reactive({
    message("Processing data")
    if (!isTRUE(all.equal(empty_table, preprocessed_data()))){
      df <- preprocessed_data()
      updateSliderInput(session,'Log2FC',max=ceiling(max(abs(df$Log2FC))))
      # Adapt the plot limit when user zoom in it
      if (Zoom()$Zoomed==T){
        df <- df[df$Log2FC>=Zoom()$coords[1]&df$Log2FC<=Zoom()$coords[2]&(-log(df$padj))>=Zoom()$coords[3]&(-log(df$padj))<=Zoom()$coords[4],]
      }
      # Update the selected points according to the selection mode
      if (selection_mode() == "Brush") {
        df$selected <- ifelse(df$GeneID%in%brushedPoints(df, brush_info()())$GeneID, "TRUE", "FALSE")
      } else if (selection_mode() == "Sliders"){
        df$selected <- ifelse((df$Log2FC>input$Log2FC&df$padj<input$pval)|(df$Log2FC<(-input$Log2FC)&df$padj<input$pval), "TRUE", "FALSE")
      }
      updateNumericInput(session, "export_GeneNumber", value=nrow(df), min=0, max=nrow(df))
      return(df)
    } else {
      return(empty_table)
    }
  })
  
  output$go_gsea_enrichplot <- renderPlot({
    if (!is.null(go_gsea_results()) && length(input$go_selected_gsea) > 0 && input$go_selected_gsea[1] != "None") {
      data = go_gsea_results()
      plot = gseaplot2(data, as.numeric(input$go_selected_gsea), 
                       pvalue_table = TRUE,
                       color = c("#E495A5", "#86B875", "#7DB0DD"))
      return(plot)
    } else {
      return(empty_gsea_plot)
    }
  })
  
    # Reactive function building the plot
    plot <- reactive({
      message("WDI: loading volcano plot")
      plot = ggplot(processed_data(), aes(x=Log2FC, y=minuslog10, col=selected)) +
        geom_point() +
        scale_color_manual(values=c("FALSE" = "#384246", "TRUE" = "#E69F00")) +
        theme(legend.position = "none") +
        xlab("log2(FoldChange)") +
        ylab("-log10(p-value ajusted)") +
        theme(text = element_text(size = 14))    
      return(plot)})
    
    pathway_gsea_dotplot <- reactive({
      if (!is.null(pathway_gsea_enrichment())) {
        data = pathway_gsea_enrichment()$enrichment
        plot = dotplot(data, showCategory=30) + ggtitle("Dotplot" )  
        return(plot)
      } else {
        return(empty_plot)
      }}
      )
    
    pathway_gsea_enrichplot <- reactive({
      if ((!is.null(input$pathway_selected_gsea)) & (length(input$pathway_selected_gsea)>0)) {
        data = pathway_gsea_enrichment()$enrichment
        plot = gseaplot2(data, as.numeric(input$pathway_selected_gsea)) 
        return(plot)
      } else {
        return(empty_gsea_plot)
      }
     })
    
    pathway_gsea_upsetplot <- reactive({
      if (!is.null(pathway_gsea_enrichment())) {
        plot = upsetplot(pathway_gsea_enrichment()$enrichment, col="#7e3535") 
        return(plot)
      } else {
        return(empty_plot)
      }
    })
    
    pathway_ora_upsetplot <- reactive({
      if (!is.null(pathway_ora_enrichment())) {
        plot = upsetplot(pathway_ora_enrichment()$enrichment) 
        return(plot)
      } else {
        return(empty_plot)
      }
    })
    
    pathway_ora_dotplot <- reactive({
      if (!is.null(pathway_ora_enrichment())) {
        data = pathway_ora_enrichment()$enrichment
        plot = dotplot(data, showCategory=30) + ggtitle("Dotplot" )  
        return(plot)
      } else {
        return(empty_plot)
      }}
    )
    
    pathway_ora_treeplot <- reactive({
      data <- pathway_ora_enrichment()
      if (is.null(data)) {
        return(empty_plot)
      }
      dataR <- pairwise_termsim(setReadable(data$enrichment,orgs[orgs$organism == input$species, "db"],keyType = "ENTREZID"))
      tryCatch({treeplot(dataR)}, error = function(e) { not_enough_data_plot })
    })
    
    
    pathway_ora_cnetplot <- reactive({
      if (!is.null(pathway_ora_enrichment())) {
        data = pathway_ora_enrichment()$enrichment
        dataR = setReadable(data, orgs[orgs$organism==input$species, "db"], 'ENTREZID')
        plot = cnetplot(dataR, layout.params = list(layout = "kk"))
        return(plot)
      } else {
        return(empty_plot)
      }
    })
    
    pathway_ora_emapplot <- reactive({
      if (!is.null(pathway_ora_enrichment())) {
        data = pathway_ora_enrichment()$enrichment
        dataR = setReadable(data, orgs[orgs$organism==input$species, "db"], 'ENTREZID')
        plot = tryCatch({emapplot(pairwise_termsim(dataR), layout.params = list(layout = "kk"), showCategory = 15)}, error = function(e) { not_enough_data_plot })
        return(plot)
      } else {
        return(empty_plot)
      }
    })
    
    pathway_gsea_emapplot <- reactive({
      if (!is.null(pathway_gsea_enrichment())) {
        data = pathway_gsea_enrichment()$enrichment
        dataR = setReadable(data, orgs[orgs$organism==input$species, "db"], 'ENTREZID')
        plot = tryCatch({emapplot(pairwise_termsim(dataR), layout.params = list(layout = "kk"), showCategory = 15)}, error = function(e) { not_enough_data_plot })
        return(plot)
      } else {
        return(empty_plot)
      }
    })
    
    pathway_gsea_treeplot <- reactive({
      data <- pathway_gsea_enrichment()
      if (is.null(data)) {
        return(empty_plot)
      }
      dataR <- pairwise_termsim(setReadable(data$enrichment,orgs[orgs$organism == input$species, "db"],keyType = "ENTREZID"))
      tryCatch({treeplot(dataR, foldChange=data$geneList)}, error = function(e) { not_enough_data_plot })
    })
    
    pathway_gsea_cnetplot <- reactive({
      if (!is.null(pathway_gsea_enrichment())) {
        data = pathway_gsea_enrichment()$enrichment
        dataR = setReadable(data, orgs[orgs$organism==input$species, "db"], 'ENTREZID')
        plot = cnetplot(dataR, foldChange=data@geneList, layout.params = list(layout = "kk"))
        return(plot)
      } else {
        return(empty_plot)
      }
    })
    
    observeEvent(input$pathway_selected_gsea, {
      if (!is.null(pathway_gsea_enrichment()) & (input$pathway_selected_gsea[1]=="None")) {
        data = pathway_gsea_enrichment()$enrichment
        if (pathwayEnrichment()$parameters$DB=="KEGG"){
          updateSelectInput(session ,"pathway_selected_gsea", choices = setNames(1:length(data$ID), data$Description))
        }
        else if (pathway_gsea_enrichment()$parameters$DB=="Reactome"){
          updateSelectInput(session ,"pathway_selected_gsea", choices = setNames(1:length(data$Description), data$Description))
        }}
    })
    
    observe({
      if (!is.null(go_ora_results())) { 
        showTab(inputId = "go_tabs", target = "Results (ORA)", session = session)
      } else {
        hideTab(inputId = "go_tabs", target = "Results (ORA)", session = session)
      }
      if (!is.null(go_gsea_results())) { 
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
    
    observeEvent(selection_mode(), {
      message(paste("WDI selection mode:", selection_mode()))
      if (selection_mode() != "None") {
        enable('Download')
        enable('download_table')
        enable('SelectAll')
        enable('ResetButton')
        removeClass("Download", "disabled-button")
        removeClass("download_table", "disabled-button")
        removeClass("SelectAll", "disabled-button")
        removeClass("ResetButton", "disabled-button")
        if (selection_mode() == "Sliders"){
          enable('pval')
          enable('Log2FC')
          if (Zoom()$Zoomed==F){
            disable('ZoomButton')
            addClass("ZoomButton", "disabled-button")
          } else {
            updateActionButton(session, "ZoomButton", label = "Unzoom selection", icon=icon('zoom-out', lib='glyphicon'))
          }
        } else if (selection_mode() == "Brush"){
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
        addClass("download_table", "disabled-button")
      }})
  
  # Download button event
  output$Download <- downloadHandler(
    filename <- function() { paste(unlist(strsplit(input$table_input[,1], ".", fixed=T))[1], "_HEATraNplot.pdf") },
    content <- function(file) {
      ggsave(file, plot(), units = "mm", height=210, width=297)
    }
  )
  # Export table event
  output$download_table <- downloadHandler(
    filename <- function() { paste(unlist(strsplit(input$table_input[,1], ".", fixed=T))[1], "_HEATraNtable.csv") },
    content <- function(file) {
      write.csv(processed_data()[,-c("selected","minuslog10")], file)
    }
  )
  
  # Reset button event
  observe({
    updateSliderInput(session,'Log2FC',value = 1)
    updateSliderInput(session,'pval',value = 0.05)
    if (selection_mode() == "Brush"){
      selection_mode("Sliders")
      session$resetBrush("plot_brush")
      brush_info(NULL)
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
    if (selection_mode() == "Brush"){
      selection_mode("Sliders")
      session$resetBrush("plot_brush")
      brush_info(NULL)
      updateActionButton(session, "ZoomButton", label = "Unzoom selection", icon=icon('zoom-out', lib='glyphicon'))
      if (Zoom()$Zoomed==F){
        disable('ZoomButton')
        addClass("ZoomButton", "disabled-button")
      }
    }
  }) |> bindEvent(input$SelectAll)
  
  # Brush use
  observe({
    message("WDI: brush use")
    selection_mode("Brush")
    brush_info(reactiveVal(input$plot_brush))
  }) |> bindEvent(input$plot_brush)
  
  # Zoom button event
  Zoom = reactive({
    if (selection_mode() == "Brush"){
      message("WDI: zoom")
      coords <- c(brush_info()()$xmin, brush_info()()$xmax, brush_info()()$ymin, brush_info()()$ymax)
      session$resetBrush("plot_brush")
      brush_info(NULL)
      selection_mode("Sliders")
      return(list(Zoomed=T, coords=coords))}
    else {
      message("WDI: unzoom")
      disable('ZoomButton')
      addClass("ZoomButton", "disabled-button")
      return(list(Zoomed=F, coords=NULL))
    }
  }) |> bindEvent(input$ZoomButton, ignoreNULL=F)
  
  observe({
    if (!is.null(pathway_ora_enrichment()) || !is.null(pathway_gsea_enrichment())) {
      if (!is.null(pathway_ora_enrichment())) {
        showTab(inputId = "pathway_tabs", target = "Results (ORA)", session = session)
        showTab(inputId = "pathway_tabs", target = "View Pathway", session = session)
      } 
      if (!is.null(pathway_gsea_enrichment())) {
        showTab(inputId = "pathway_tabs", target = "Results (GSEA)", session = session)
        showTab(inputId = "pathway_tabs", target = "View Pathway", session = session)
      }
    }
  })

 #--- GO Enrichment ---#
  
  observeEvent(input$go_analysisButton , {
    req(preprocessed_data())
    
    withProgress(message = 'Computing GO enrichment...', {
      df <- preprocessed_data()
      
      # List of genes for the universe
      original_gene_list <- df$Log2FC
      names(original_gene_list) <- df$GeneID
      gene_list <- sort(na.omit(original_gene_list), decreasing = T)
      
      # Gene selection by adjusted p-value threshold
      sig_genes_df <- df[df$padj < input$go_pval,]
      genes <- sig_genes_df$Log2FC
      names(genes) <- sig_genes_df$GeneID
      genes <- na.omit(genes)
      
      up_genes <- names(genes)[genes > log2(input$go_ora_fc)]
      down_genes <- names(genes)[genes < -log2(input$go_ora_fc)]
      
      # Selecting GO ontology and organism
      ontology <- input$inputGO
      go_organism <- orgs[orgs$organism==input$species, "db"]
      universe <- names(gene_list)
      library(go_organism, character.only = T)
      
      if ("Gene Set Enrichment Analysis (GSEA)" %in% input$go_analysisMethodChoice) {
        setProgress(message("Computing GSEA enrichment for GO term..."), value=0.5)
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
          go_gsea_results(gsea_result) 
        } else {
          go_gsea_results(NULL)
        }
        
      } 
      
      if ("Over Representation Analysis (ORA)" %in% input$go_analysisMethodChoice) {
        setProgress(message("Computing ORA enrichment for GO term..."), value=0.75)
        analyze_up <- "Over expressed DEG" %in% input$go_ora_choice
        analyze_down <- "Under expressed DEG" %in% input$go_ora_choice
        
        if (!analyze_up && !analyze_down) {
          shinyalert("Error", "Select at least one type of gene expression (Over or Under expressed DEG).", type = "error")
          return(NULL)
        }
        
        if (analyze_up && analyze_down) {
          # Both
          genes_to_analyze <- union(up_genes, down_genes)
          analysis_type <- "both"
        } else if (analyze_up) {
          # Only up-regulated
          genes_to_analyze <- up_genes
          analysis_type <- "up"
        } else {
          # Only down-regulated
          genes_to_analyze <- down_genes
          analysis_type <- "down"
        }
        
        ratio = (length(genes_to_analyze) / length(universe))
        message(paste("ORA: selected genes represent ", round(ratio*100,2), "% of Universe", sep=""))
        if (ratio > 0.1) {
          shinyalert("Statistical issue", text=paste("The genes selected for the ORA analysis represent more than 10% of the initial universe (", round(ratio*100, 2), "%), which can greatly distort and negatively impact the statistical analysis. Please increase the minimum fold change threshold or decrease the maximum p-value threshold.", sep=""), type = "error")
          go_ora_results(NULL)
          return()
        }
        
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
            go_ora_results(list(result = result_go, type = analysis_type))
            
            all_genes_binary <- as.integer(names(gene_list) %in% genes_to_analyze)
            names(all_genes_binary) <- names(gene_list)
            
            go_data <- new("topGOdata",
                          ontology = ontology,
                          allGenes = all_genes_binary,
                          geneSelectionFun = function(x)(x < 0.05),
                          annot = annFUN.org,
                          mapping = go_organism,
                          ID = "ensembl")
            
            sig_go(go_data)
            sig_go_result(runTest(sig_go(), algorithm = "classic", statistic = "fisher"))
          }
          else {
            go_ora_results(NULL)
          }
        } else {
          go_ora_results(NULL)
        }

        }
      })
    })
  
  observeEvent(go_gsea_results(), {
    if (!is.null(go_gsea_results())) {
      data = go_gsea_results()
      updateSelectInput(session, "go_selected_gsea", 
                        choices = setNames(1:length(data$Description), 
                                           data$Description))
    }
  })
  
    #--- Pathway Enrichment ---#
    
    observe({
      if (is.na(preprocessed_data()[1,1])){
        shinyalert("Bad input", "Analysis require data!", type = "error", confirmButtonCol = "#7e3535")
      }
      else if ("ORA" %in% input$analysisMethodChoice & is.null(input$pathway_ora_choice)){
        shinyalert("Incomplete selection!", "You should select an interest for ORA method", type = "error", confirmButtonCol = "#7e3535")
      }
      else {
        withProgress(message = "Computing Pathway enrichment...", {
            message("Computing Pathway enrichment...")
            df = as.data.frame(preprocessed_data())
            if ("ORA" %in% input$analysisMethodChoice) {
              setProgress(message = "Computing ORA enrichment for Pathway...", value=0.33)
              
              # 1. Run enrichment and store in a temp var
              res <- pathway(
                df,
                organism   = input$species,
                DB         = input$dbPathwaychoice,
                analysis   = "ORA",
                ora_interest= input$pathway_ora_choice,
                pval       = input$pvalPathway,
                threshold = input$pathway_ora_fc
              )
              # 2. Test whether any enriched terms were returned
              if (is.null(res) || (is.null(res$enrichment)) || nrow(res$enrichment@result) == 0) {
                # 3a. No enriched terms: set to NULL
                shinyalert("No results for ORA", text="Check your significance and fold-change threshold, which may be too restrictive.", type = "error")
                pathway_ora_enrichment(NULL)
              } else {
                # 3b. Enriched terms found: pass result and show tab
                pathway_ora_enrichment(res)
                showTab(inputId = "pathway_tabs", target = "View Pathway", session = session)
              }
            } else {
              pathway_ora_enrichment(NULL)
            }
            
            if ("GSEA" %in% input$analysisMethodChoice) {
              setProgress(message = "Computing GSEA enrichment for Pathway...", value=0.33)
              
              res <- pathway(
                df,
                organism   = input$species,
                DB         = input$dbPathwaychoice,
                analysis   = "GSEA",
                ora_interest= input$pathway_ora_choice,
                pAdjustMethod = config$STAT$adjust_method,
                pval       = input$pvalPathway
              )
              
              if (is.null(res) || (is.null(res$enrichment)) ||  nrow(res$enrichment@result) == 0) {
                shinyalert("No results for GSEA", text="Check your significance threshold, which may be too restrictive.", type = "error")
                pathway_gsea_enrichment(NULL)
              } else {
                pathway_gsea_enrichment(res)
                showTab(inputId = "pathway_tabs", target = "View Pathway", session = session)
              }
            } else {
              pathway_gsea_enrichment(NULL)
            }
        })
    }}) |> bindEvent(input$pathway_analysis_button)
    
    observeEvent(pathway_gsea_enrichment(), {
      data = pathway_gsea_enrichment()$enrichment
      if (pathway_gsea_enrichment()$parameters$DB=="KEGG"){
        updateSelectInput(session ,"pathway_selected_gsea", choices = setNames(1:length(data$ID), data$Description))
        } else if (pathway_gsea_enrichment()$parameters$DB=="Reactome"){
        updateSelectInput(session ,"pathway_selected_gsea", choices = setNames(1:length(data$Description), data$Description))
        }
    }) 
    
    observeEvent(pathway_combined_results(), {
      combined_results <- pathway_combined_results()
      
      if (!is.null(combined_results)) {
        showTab(inputId = "pathway_tabs", target = "View Pathway", session = session)
        
        updateSelectInput(session, "pathway", 
                          choices = combined_results$choices,
                          selected = combined_results$choices[1])
      }
    })
  
  # -----------------------------------------
  # UI output 
  # -----------------------------------------

  #--- Whole Data Inspection ---#
    
  output$info_select <- renderUI({
    if (is.null(preprocessed_data) || is.na(preprocessed_data()[1,1])){
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
    message("WDI: rendering DataTable")
    if (!is.na(preprocessed_data()[1,1])){
      processed_data()[processed_data()$selected==T,-c("selected","minuslog10")]}
    else {
      processed_data()
    }})
  
  output$volcano_plot <- renderPlot ({
    message("WDI: rendering volcano plot")
    if (!is.na(preprocessed_data()[1,1])){
      plot()}
    else {
      return(empty_plot)
    }
  })
  
  #--- GO Enrichment ---#
  
  output$go_ora_barplot <- renderPlot({
    message("GO: rendering ORA Barplot")
    if (!is.null(go_ora_results()$result)) {
      barplot(go_ora_results()$result)
    } else {
      empty_plot
    }
  })
  
  output$go_ora_dotplot <- renderPlot({
    message("GO: rendering ORA Dotplot")
    if (!is.null(go_ora_results()$result)) {
      dotplot(go_ora_results()$result)
    } else {
      empty_plot
    }
  })
  
  output$go_ora_netplot <- renderPlot({
    message("GO: rendering Gene-Concept Network")
    if (!is.null(go_ora_results()$result)) {
      go_enrich <- pairwise_termsim(go_ora_results()$result)
      cnetplot(go_enrich, layout.params = list(layout = "kk"), showCategory = 15)
    } else {
      empty_plot
    }
  })
  
  output$go_ora_emapplot <- renderPlot({
    message("GO: rendering ORA Enrichment Map")
    if (!(is.null(go_ora_results()))) {
      go_enrich <- pairwise_termsim(go_ora_results()$result)
      tryCatch({emapplot(go_enrich, layout.params = list(layout = "kk"), showCategory = 15)}, error = function(e) { not_enough_data_plot })
    } else { empty_plot }
  })
  
  output$go_ora_table <- DT::renderDT({
    message("GO: rendering ORA DataTable")
    if (!(is.null(go_ora_results()$result))) {
      df = as.data.frame(go_ora_results()$result)
      for (i in 1:length(df$ID)) {
        df$ID[i] = paste('<a href="https://www.ebi.ac.uk/QuickGO/term/', df$ID[i], '" target="_blank">', df$ID[i], '</a>', sep="")
      }
      if (nrow(df) != 0) {
        datatable(df, escape = F, options = list(scrollX = TRUE, scrollY = TRUE, autoWidth = TRUE, autoHeight = TRUE))
      }
      else {
        empty_go_table}
    } else {
      empty_go_table
    }
  })
  
  output$go_ora_upsetplot <- renderPlot({
    message("GO: rendering ORA Upsetplot")
    if (!(is.null(go_ora_results()$result))) {
      upsetplot(go_ora_results()$result)
    } else {
      empty_plot
    }
    
  })
  
  output$go_gsea_upsetplot <- renderPlot({
    message("GO: rendering GSEA Upsetplot")
    if (!(is.null(go_gsea_results()))) {
      upsetplot(go_gsea_results())
    } else { empty_plot }
  })
  
  output$go_gsea_netplot <- renderPlot({
    message("GO: rendering GSEA Gene-Concept Network")
    if (!(is.null(go_gsea_results()))) {
      go_enrich <- pairwise_termsim(go_gsea_results())
      cnetplot(go_enrich, layout.params = list(layout = "kk"), showCategory = 15)
    } else { empty_plot }
  })
  
  output$go_gsea_emapPlot <- renderPlot({
    message("GO: rendering GSEA Enrichment map")
    if (!(is.null(go_gsea_results()))) {
      go_enrich <- pairwise_termsim(go_gsea_results())
      tryCatch({emapplot(go_enrich, layout.params = list(layout = "kk"), showCategory = 15)}, error = function(e) { not_enough_data_plot })
    } else { empty_plot }
  })
  
  output$go_gsea_dotplot <- renderPlot({
    message("GO: rendering GSEA Dotplot")
    if (!(is.null(go_gsea_results()))) {
      dotplot(go_gsea_results())
    } else { empty_plot }
  })
  
  output$go_gsea_ridgeplot <- renderPlot({
    message("GO: rendering GSEA Ridgeplot")
    if (!(is.null(go_gsea_results()))) {
      ridgeplot(go_gsea_results(), showCategory = 13)
    } else { empty_plot }
  })
  
  output$go_gsea_table <- DT::renderDT({
    message("GO: rendering GSEA DataTable")
    if (!(is.null(go_gsea_results()))) {
      df = as.data.frame(go_gsea_results())
      for (i in 1:length(df$ID)) {
        df$ID[i] = paste('<a href="https://www.ebi.ac.uk/QuickGO/term/', df$ID[i], '" target="_blank">', df$ID[i], '</a>', sep="")
      }
      if (nrow(df) != 0) {
        datatable(df, escape = F, options = list(scrollX = TRUE, scrollY = TRUE, autoWidth = TRUE, autoHeight = TRUE))
      } else {
        empty_go_table}
    } else {
      empty_go_table
    }
  })
  
  output$go_ora_goplot <- renderPlot({
    message("GO: rendering ORA GOplot")
    if (!is.null(sig_go())) {
      showSigOfNodes(sig_go(), score(sig_go_result()), firstSigNodes = input$goLevel, useInfo = "def")
    } else {
      empty_plot
    }
    
  })
  
  #--- Pathway Enrichment ---#
    
    output$pathway_gsea_table <- DT::renderDT ({
      message("Pathway: Rendering GSEA DataTable")
      if (!is.null(pathway_gsea_enrichment())){
        df = as.data.frame(pathway_gsea_enrichment()$enrichment)
        for (i in 1:length(df$ID)) {
          if (pathway_gsea_enrichment()$parameters$DB=="KEGG"){
            df$ID[i] = paste('<a href="https://www.kegg.jp/entry/', df$ID[i], '" target="_blank">', df$ID[i], '</a>', sep="")
          } else if (pathway_gsea_enrichment()$parameters$DB=="Reactome"){
            df$ID[i] = paste('<a href="https://reactome.org/PathwayBrowser/#/', df$ID[i], '" target="_blank">', df$ID[i], '</a>', sep="")
          }
        }
        if (nrow(df) != 0) {
          datatable(df, escape = F, options = list(scrollX = TRUE, scrollY = TRUE, autoWidth = TRUE, autoHeight = TRUE))
        }
        else {
          empty_pathway_table}
       }
      else {
        empty_pathway_table
      }
    })
    
    output$pathway_ora_table <- DT::renderDT ({
      message("Pathway: Rendering ORA DataTable")
      if (!is.null(pathway_ora_enrichment())){
        df = as.data.frame(pathway_ora_enrichment()$enrichment)
        for (i in 1:length(df$ID)) {
          if (pathway_ora_enrichment()$parameters$DB=="KEGG"){
            df$ID[i] = paste('<a href="https://www.genome.jp/pathway/', df$ID[i], '" target="_blank">', df$ID[i], '</a>', sep="")
          } else if (pathway_ora_enrichment()$parameters$DB=="Reactome"){
            df$ID[i] = paste('<a href="https://reactome.org/PathwayBrowser/#/', df$ID[i], '" target="_blank">', df$ID[i], '</a>', sep="")
          }
        }
        if (nrow(df) != 0) {
          datatable(df, escape = F, options = list(scrollX = TRUE, scrollY = TRUE, autoWidth = TRUE, autoHeight = TRUE))
        }
        else {
          empty_pathway_table}
      }
      else {
        empty_pathway_table
      }
    })
    
    output$pathway_gsea_upsetplot <- renderPlot ({
      message("Pathway: rendering GSEA Upsetplot")
      if (!is.null(pathway_gsea_upsetplot())){
        pathway_gsea_upsetplot()}
      else {
      }
    })
    
    output$pathway_ora_upsetplot <- renderPlot ({
      message("Pathway: rendering ORA Upsetplot")
      if (!is.null(pathway_ora_upsetplot())){
        pathway_ora_upsetplot()}
      else {
      }
    })
    
    output$pathway_gsea_dotplot <- renderPlot ({
      message("Pathway: rendering GSEA Dotplot")
      if (!is.null(pathway_gsea_dotplot())){
        pathway_gsea_dotplot()}
      else {
      }
    })
    
    output$pathway_ora_dotplot <- renderPlot ({
      message("Pathway: rendering ORA Dotplot")
      if (!is.null(pathway_ora_dotplot())){
        pathway_ora_dotplot()}
      else {
      }
    })
      
    output$pathway_ora_treeplot <- renderPlot ({
      message("Pathway: rendering ORA Treeplot")
      if (!is.null(pathway_ora_dotplot())){
        pathway_ora_treeplot()}
      else {
      }
    })
    
    output$pathway_ora_cnetplot <- renderPlot ({
      message("Pathway: rendering ORA Gene-Concept Network")
      if (!is.null(pathway_ora_cnetplot())){
        pathway_ora_cnetplot()}
      else {
      }
    })  
    
    output$pathway_ora_emapplot <- renderPlot ({
      message("Pathway: rendering ORA Enrichment map")
      if (!is.null(pathway_ora_emapplot())){
        pathway_ora_emapplot()}
      else {
      }
    })  
    
    output$pathway_gsea_emapplot <- renderPlot ({
      message("Pathway: rendering GSEA Enrichment map")
      if (!is.null(pathway_gsea_emapplot())){
        pathway_gsea_emapplot()}
      else {
      }
    })  
    
    output$pathway_gsea_treeplot <- renderPlot ({
      message("Pathway: rendering ORA Treeplot")
      if (!is.null(pathway_gsea_treeplot())){
        pathway_gsea_treeplot()}
      else {
      }
    })
    
    output$pathway_gsea_cnetplot <- renderPlot ({
      message("Pathway: rendering GSEA Gene-Concept Network")
      if (!is.null(pathway_gsea_cnetplot())){
        pathway_gsea_cnetplot()}
      else {
      }
    })  
    
    output$pathway_gsea_enrichplot <- renderPlot ({
      message("Pathway: rendering GSEA Enrichplot")
      if (!is.null(pathway_gsea_enrichplot()) & pathway_gsea_enrichment()$parameters$analysis=="GSEA"){
        pathway_gsea_enrichplot()}
      else {
      }
  })
    
    pathway_combined_results <- reactive({
      ora_results <- pathway_ora_enrichment()
      gsea_results <- pathway_gsea_enrichment()
      
      common_cols <- c("ID", "Description", "pvalue", "p.adjust")
      
      if (!is.null(ora_results) && !is.null(gsea_results)) {
        
        ora_data <- ora_results$enrichment[, common_cols, drop = FALSE]
        gsea_data <- gsea_results$enrichment[, common_cols, drop = FALSE]
        
        ora_data$analysis_type <- "ORA"
        gsea_data$analysis_type <- "GSEA"
        
        combined_data <- rbind(ora_data, gsea_data)
        
        if (ora_results$parameters$DB == "KEGG") {
          choices <- setNames(combined_data$ID, 
                              paste0(combined_data$Description, " (", combined_data$analysis_type, ")"))
        } else { 
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
        ora_data <- ora_results$enrichment[, common_cols, drop = FALSE]
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
        gsea_data <- gsea_results$enrichment[, common_cols, drop = FALSE]
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
    
    output$pathway_image <- renderImage({
      combined_results <- pathway_combined_results()
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
              message("Pathway: displaying", image_path)
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
        ggsave("./out/empty.png", plot = empty_pathway_plot)
        return(list(src = "./out/empty.png", contentType = 'image/png', width="100%", alt="Please select at least one enriched pathway"))
      }
      
      return(list(src = "", alt = "No pathway image available"))
    }, deleteFile = FALSE)
    
  # -----------------------------------------
  # Save and export functions
  # -----------------------------------------

    output$export_options <- renderUI({
      if ((!is.null(preprocessed_data()) && !is.na(preprocessed_data()[1,1])) 
          || (!is.null(go_ora_results()$result)) 
          || (!is.null(go_gsea_results()))
          || (!is.null(pathway_ora_enrichment())) 
          || (!is.null(pathway_gsea_enrichment()))) {
        
        choices <- character(0)
        if (!is.null(preprocessed_data()) && !is.na(preprocessed_data()[1,1])) {
          choices <- "Whole Data Inspection (volcano_plot + Table)"
        }
        if (!is.null(go_ora_results()$result))    choices <- c(choices, "GO ORA")
        if (!is.null(go_gsea_results()))          choices <- c(choices, "GO GSEA")
        if (!is.null(pathway_ora_enrichment()))   choices <- c(choices, "Pathway ORA")
        if (!is.null(pathway_gsea_enrichment()))  choices <- c(choices, "Pathway GSEA")
      
        ui_elems <- list()
        
        if (length(choices) == 1 && choices == "Whole Data Inspection (volcano_plot + Table)") {
          ui_elems <- tagList(
            ui_elems,
            HTML("<em>No enrichment was carried out. 
              The current export will only contain Whole Data Inspection.</em><br/>")
          )
        } else {
          ui_elems <- tagList(
            ui_elems,
            checkboxGroupInput("export_choices", "Analysis to export:", choices = choices)
          )
        }

        ui_elems <- tagList(
          ui_elems,
          downloadButton("exportReport", "Download report")
        )
        ui_elems
        
      } else {
        HTML("<em>Only already-launched analyses can be exported.</em><br/>")
      }
    })

    output$exportReport <- downloadHandler(
      filename = function() {
        paste0("HEATraN_results_", Sys.Date(), ".html")
      },
      content = function(file) {
        withProgress(message = 'Generating report...', value = 0.2, {
          temp_report <- file.path(tempdir(), "template.Rmd")
          file.copy("www/template.Rmd", temp_report, overwrite = TRUE)
          
          rmarkdown::render(
            input = temp_report,
            output_file = file,
            params = list(
              # General parameters
              nGenes = 20,
              organismInfo = input$species,
              
              info_select = tryCatch({
                if (!is.null(preprocessed_data()) && !is.na(preprocessed_data()[1,1])) {
                  if (is.null(input$plot_brush)) {
                    paste("Current selection - p-value: [0 ; ", input$pval,"] log2(FoldChange): [-", input$Log2FC, " ; ", input$Log2FC,"]", sep="")
                  } else {
                    paste("Current selection - p-value: [", format(exp(-input$plot_brush$ymax), scientific=T, digits=3), " ; " , format(exp(-input$plot_brush$ymin), scientific=T, digits=3),"] log2(FoldChange): [", round(input$plot_brush$xmin, 3), " ; ", round(input$plot_brush$xmax, 3),"]", sep="")
                  }
                } else NULL
              }, error = function(e) NULL),
              
              volcano_plot = tryCatch({
                if (!is.null(preprocessed_data()) && !is.na(preprocessed_data()[1,1])) {
                  p <- plot()  
                  p  
                } else NULL
              }, error = function(e) NULL),
              
              table = tryCatch({
                if (!is.null(preprocessed_data()) && !is.na(preprocessed_data()[1,1])) {
                  data <- processed_data()
                  if (!is.null(data) && nrow(data) > 0) {
                    data[data$selected == TRUE, -c("selected", "minuslog10")]
                  } else NULL
                } else NULL
              }, error = function(e) NULL),
              
              # GO ORA outputs - Capture des fonctions réactives
              go_ora_barplot = tryCatch({
                if (!is.null(go_ora_results()) && !is.null(go_ora_results()$result)) {
                  barplot(go_ora_results()$result)
                } else NULL
              }, error = function(e) NULL),
              
              go_ora_dotplot = tryCatch({
                if (!is.null(go_ora_results()) && !is.null(go_ora_results()$result)) {
                  dotplot(go_ora_results()$result)
                } else NULL
              }, error = function(e) NULL),
              
              go_ora_emapplot = tryCatch({
                if (!is.null(go_ora_results()) && !is.null(go_ora_results()$result)) {
                  go_enrich <- pairwise_termsim(go_ora_results()$result)
                  tryCatch({emapplot(pairwise_termsim(dataR), layout.params = list(layout = "kk"), showCategory = 15)}, error = function(e) { not_enough_data_plot })
                } else NULL
              }, error = function(e) NULL),
              
              go_ora_netplot = tryCatch({
                if (!is.null(go_ora_results()) && !is.null(go_ora_results()$result)) {
                  go_enrich <- pairwise_termsim(go_ora_results()$result)
                  tryCatch({emapplot(pairwise_termsim(dataR), layout.params = list(layout = "kk"), showCategory = 15)}, error = function(e) { not_enough_data_plot })
                } else NULL
              }, error = function(e) NULL),
              
              go_ora_table = tryCatch({
                if (!is.null(go_ora_results()) && !is.null(go_ora_results()$result)) {
                  as.data.frame(go_ora_results()$result)
                } else NULL
              }, error = function(e) NULL),
              
              go_ora_upsetplot = tryCatch({
                if (!is.null(go_ora_results()) && !is.null(go_ora_results()$result)) {
                  upsetplot(go_ora_results()$result)
                } else NULL
              }, error = function(e) NULL),
              
              go_ora_goplot = tryCatch({
                if (!is.null(sig_go()) && !is.null(sig_go_result())) {
                  img_file <- tempfile(fileext = ".png")
                  png(filename = img_file, width = 800, height = 600)
                  showSigOfNodes(sig_go(), score(sig_go_result()), firstSigNodes = 3)
                  dev.off()
                  list(src = img_file, contentType = "image/png", width = "100%")
                } else {
                  NULL
                }
              }, error = function(e) {
                NULL
              }),
              
              go_gsea_upsetplot = tryCatch({
                if (!is.null(go_gsea_results())) {
                  upsetplot(go_gsea_results())
                } else NULL
              }, error = function(e) NULL),
              
              go_gsea_netplot = tryCatch({
                if (!is.null(go_gsea_results())) {
                  go_enrich <- pairwise_termsim(go_gsea_results())
                  cnetplot(go_enrich, layout.params = list(layout = "kk"), showCategory = 15)
                } else NULL
              }, error = function(e) NULL),
              
              go_gsea_emapPlot = tryCatch({
                if (!is.null(go_gsea_results())) {
                  go_enrich <- pairwise_termsim(go_gsea_results())
                  emapplot(go_enrich, layout.params = list(layout = "kk"), showCategory = 15)
                } else NULL
              }, error = function(e) NULL),
              
              go_gsea_dotplot = tryCatch({
                if (!is.null(go_gsea_results())) {
                  dotplot(go_gsea_results())
                } else NULL
              }, error = function(e) NULL),
              
              go_gsea_ridgeplot = tryCatch({
                if (!is.null(go_gsea_results())) {
                  ridgeplot(go_gsea_results(), showCategory = 13)
                } else NULL
              }, error = function(e) NULL),
              
              go_gsea_table = tryCatch({
                if (!is.null(go_gsea_results())) {
                  as.data.frame(go_gsea_results())
                } else NULL
              }, error = function(e) NULL),
              
              # Pathway ORA outputs
              pathway_ora_table = tryCatch({
                if (!is.null(pathway_ora_enrichment())) {
                  as.data.frame(pathway_ora_enrichment()$enrichment)
                } else NULL
              }, error = function(e) NULL),
              
              pathway_ora_upsetplot = tryCatch({
                if (!is.null(pathway_ora_enrichment())) {
                  p <- pathway_ora_upsetplot()
                  p
                } else NULL
              }, error = function(e) NULL),
              
              pathway_ora_emapplot = tryCatch({
                if (!is.null(pathway_ora_emapplot())) {
                  p <- pathway_ora_emapplot()
                  p
                } else NULL
              }, error = function(e) NULL),
              
              pathway_ora_dotplot = tryCatch({
                if (!is.null(pathway_ora_enrichment())) {
                  p <- pathway_ora_dotplot()
                  p
                } else NULL
              }, error = function(e) NULL),
              
              pathway_ora_treeplot = tryCatch({
                if (!is.null(pathway_ora_enrichment())) {
                  p <- pathway_ora_treeplot()
                  p
                } else NULL
              }, error = function(e) NULL),
              
              pathway_ora_cnetplot = tryCatch({
                if (!is.null(pathway_ora_enrichment())) {
                  p <- pathway_ora_cnetplot()
                  p
                } else NULL
              }, error = function(e) NULL),
              
              # Pathway GSEA outputs
              pathway_gsea_table = tryCatch({
                if (!is.null(pathway_gsea_enrichment())) {
                  as.data.frame(pathway_gsea_enrichment()$enrichment)
                } else NULL
              }, error = function(e) NULL),
              
              pathway_gsea_upsetplot = tryCatch({
                if (!is.null(pathway_gsea_enrichment())) {
                  p <- pathway_gsea_upsetplot()
                  p
                } else NULL
              }, error = function(e) NULL),
              
              pathway_gsea_dotplot = tryCatch({
                if (!is.null(pathway_gsea_enrichment())) {
                  p <- pathway_gsea_dotplot()
                  p
                } else NULL
              }, error = function(e) NULL),
              
              pathway_gsea_treeplot = tryCatch({
                if (!is.null(pathway_gsea_enrichment())) {
                  p <- pathway_gsea_treeplot()
                  p
                } else NULL
              }, error = function(e) NULL),
              
              pathway_gsea_cnetplot = tryCatch({
                if (!is.null(pathway_gsea_enrichment())) {
                  p <- pathway_gsea_cnetplot()
                  p
                } else NULL
              }, error = function(e) NULL),
              
              pathway_gsea_emapplot = tryCatch({
                if (!is.null(pathway_gsea_emapplot())) {
                  p <- pathway_gsea_emapplot()
                  p
                } else NULL
              }, error = function(e) NULL),
              
              go_ora_results = go_ora_results(),
              go_gsea_results = go_gsea_results(),
              pathway_gsea_enrichment = pathway_gsea_enrichment(),
              pathway_ora_enrichment = pathway_ora_enrichment(),
              preprocessed_data = preprocessed_data(),
              processed_data = processed_data(),
              export_choices     = input$export_choices,
              include_gsea_plots   = input$include_gsea_plots
            ),
            knit_root_dir = normalizePath(getwd()),
            envir = new.env(parent = globalenv())
          )
        })
      })
  
    session$onSessionEnded(function() {
      stopApp()
    })
    
}

