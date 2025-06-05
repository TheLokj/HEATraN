# Developed by DAHER Rayan, NAIT EL DJOUDI Lamia & LESAGE Louison (@thelokj).
# rayan.daher@univ-rouen.fr
# lamia.nait-el-djoudi@univ-rouen.fr
# louison.lesage@univ-rouen.fr
# Students at Rouen Normandy University
# Master of Bioinformatics, class M2.2 BIMS 2026 
# Last updated : 05/06/2025
# HEATraN version 1.0.0

# Theme definition
HEATraN_theme <- create_theme(
  adminlte_color(
    red = "#7e3535"
  )
)

fea_page <- readChar("./www/doc.html", file.info("./www/doc.html")$size)
stat_page <- readChar("./www/stat.html", file.info("./www/stat.html")$size)

spinner <- function(plot, type = 6, color = "#ef940b") {
  shinycssloaders::withSpinner(plot, type = type, color = color)
}

# -----------------------------------------
# Dashboard creation
# -----------------------------------------
ui <- dashboardPage(skin="red", header = dashboardHeader(title= HTML("<b style='font-size:26px; color:#ebb233; font-weight:900'>HEATraN</b>")), 
              
              # -----------------------------------------
              # Sidebar elements
              # -----------------------------------------
              sidebar = dashboardSidebar(
                
                sidebarMenu(menuItem("Home", tabName = "HOME", icon = icon("home"))),
                
                # Reference organism selection
                selectInput("species", "Select organism name:",
                            c("Homo sapiens",
                              "Mus musculus",
                              "Rattus norvegicus",
                              "Drosophila melanogaster",
                              "Escherichia coli (strain K12)",
                              "Saccharomyces cerevisiae",
                              "Danio rerio",
                              "Caenorhabditis elegans",
                              "Arabidopsis thaliana",
                              "Gallus gallus",
                              "Sus scrofa",
                              "Canis lupus familiaris",
                              "Bos taurus",
                              "Xenopus laevis")),
                
                # File importation: csv, tsv, xls and xlsx
                fileInput("table_input", 
                          "Import associated data:",
                          accept = c(
                            "text/csv",
                            "text/tsv",
                            "application/vnd.ms-excel",
                            "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                            "text/comma-separated-values,text/plain",
                            "text/tab-separated-values",
                            ".csv",
                            ".tsv")),
                
                # Reference organism selection
                
                # Menu
                sidebarMenu(
                  menuItem("Documentation", icon = icon('book', lib='glyphicon'), tabName = "DOC"),
                  menuItem("Whole Data Inspection", icon = icon('eye-open', lib='glyphicon'), tabName = "WDI"),
                  menuItem("GO Term Enrichment", icon = icon('text-background', lib='glyphicon'), tabName = "GO_analysis"),
                  menuItem("Pathways Enrichment", icon = icon('transfer', lib='glyphicon'),  tabName = "PATH"),
                  menuItem("Export results", icon = icon('share', lib='glyphicon'),  tabName = "EXPORT")
                  )
                ),
              
              # -----------------------------------------
              # Body tabs 
              # -----------------------------------------
              body = dashboardBody( 
                
                # Load style and theme
                use_theme(HEATraN_theme),
                includeCSS("./www/style.css"),
                
                tabItems(
                  tabItem(tabName = "HOME",
                          HTML('<div style="padding: 20px; max-width: 100%; box-sizing: border-box;">
                  <b style="font-size: 40px; color:#7e3535; font-weight:900; display: block; text-align: center;">HEATraN</b>
                  <br/>
                  <p style="font-size: 16px; text-align: justify;">
                      <b>HEATraN</b> (litteraly <i><b>H</b>yper-<b>E</b>xpression <b>A</b>nalysis <b>T</b>ool <b>ra</b>mpantly developed in <b>N</b>ormandy</i>) is a bioinformatics analysis tool dedicated to transcriptomic analysis. It was developed as part of a student project in the Bioinformatics Master of Rouen Normandy University under the supervision of Dauchel Hélène, Education Manager.
                      You can find its last version on its <a style="font-weight: bold;" href="https://github.com/TheLokj/HEATraN">GitHub</a>.<br/> 
                  </p>
                  <br/>
                  <div style="text-align: center;">
                       <img src="logo.png" style="max-width: 35%; height: auto;" alt="HEATraN logo">
                  </div>
                  <br/>
                  <div style="text-align: center; font-size: 16px;">
                      <em>HEATraN by DAHER Rayan, NAIT EL DJOUDI Lamia, LESAGE Louison. University of Rouen Normandy.</em>
                  </div>
                  <br/>
                  <div style="text-align:right; position: fixed; bottom: 20px; right: 20px; font-size: 16px; background-color: white; font-size: 12px; padding: 10px; border-radius: 5px; box-shadow: 0 0 10px rgba(0,0,0,0.1);">
                      <em>1.0.0</em>
                  </div>
              </div>'),
                          
                  ),
                  tabItem(tabName = "DOC",
                          tabsetPanel(id="doctab", 
                                      tabPanel(title="Functional Enrichment Analysis", HTML(fea_page)),
                                      tabPanel(title="Statistics", HTML(stat_page),
                                      fluidRow(
                                        box(title="Expert statistician's parameters", width=12,
                                          HTML(""),
                                          selectInput("ajust_method", label="Select p-value ajust method:",
                                                      choices=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"),
                                                      selected = read.ini("./conf.ini")$STAT$adjust_method),
                                          HTML("<br/>ORA enrichment is constrained both by the adjusted p-value threshold (using the method above) and by a q-value threshold.<br/>"),
                                          numericInput("q_val", label="Select a q-value threshold", min=0, max=1, value = read.ini("./conf.ini")$STAT$q_val, step=0.01))
                                        )))
                  ),
                  
                  # Whole Data Inspection tab
                  tabItem(tabName = "WDI",
                          h2("Whole Data Inspection"),
                          
                          fluidRow (
                            box(
                              title = HTML("<b>Volcano plot</b>"),
                              id = "volcano_plot", height = "500px",
                              plotOutput("volcano_plot", 
                                         height="425px",
                                         click = "plot_click",
                                         dblclick = "plot_dblclick",
                                         hover = "plot_hover",
                                         brush = brushOpts(id = "plot_brush", delay = 3000, delayType = "debounce", fill="#7e3535", stroke="#7e3535")
                              )
                            ),
                            box(
                              # Use ShinyJS to allow deactivation of download button
                              shinyjs::useShinyjs(),
                              title = HTML("<b>Selection options</b>"), side = "right",
                              id = "tabset2", height = "500px",
                              htmlOutput("info_select"),
                              sliderInput("pval", "Maximum ajusted p-value", min = 0, max = 1, value = 0.05, round=F, step=0.001),
                              sliderInput("Log2FC", "Minimum absolute log2(FoldChange)", min = 0, max = 1, value = 0),
                              actionButton("SelectAll", "Select all data", icon=icon('search', lib='glyphicon')),
                              actionButton("ResetButton", "Reset selection", icon=icon('refresh', lib='glyphicon')),
                              actionButton("ZoomButton", "Zoom in the selection", icon=icon('zoom-in', lib='glyphicon')),
                              HTML("<br/><br/><b>Save options</b><br/><br/>"),
                              downloadButton("Download", "Download plot", icon=icon('download-alt', lib='glyphicon')),
                              downloadButton("download_table", "Export table", icon=icon('share', lib='glyphicon'))
                            )),
                          fluidRow (
                            box(
                              title = HTML("<b>Selected points</b>"),
                              id = "table", 
                              width = 12,
                              # Charge the table returned by the server.R
                              DT::DTOutput('table')
                            ))
                  ),
                  
                 
                  tabItem(
                    tabName = "GO_analysis", h2("GO term enrichment"),
                    tabsetPanel(
                      id = "go_tabs",
                      #width = 12,
                      tabPanel("Analysis", 
                        fluidRow(
                           box(
                             title = HTML("<i>Analysis parameters</i>"),
                             id = "gotabset1",
                             width = 12,
                             status = NULL,
                             solidHeader = FALSE,
                             collapsible = FALSE,
                             
                             # Entrée pour la p-value et le q-value
                             numericInput("go_pval", "Select an adjusted p-value cutoff", 
                                         min = 0, max = 1, value = 0.05, step=0.001),
                             
                             # Sélection de la méthode d'analyse
                             checkboxGroupInput("go_analysisMethodChoice", "Analysis method",
                                                choices = c("Over Representation Analysis (ORA)", "Gene Set Enrichment Analysis (GSEA)"),
                                                selected = NULL,
                                                inline = TRUE),
                             
                             # Sélection des gènes d'intérêt pour ORA
                             conditionalPanel(
                               condition = "input.go_analysisMethodChoice.includes('Over Representation Analysis (ORA)')",
                               sliderInput("go_ora_fc", label="ORA - Minimum Absolute Fold-Change (×) to be considered over- or under-enriched", value=1, step=0.01, min=1, max=4),
                               checkboxGroupInput("go_ora_choice", "ORA - Interest",
                                                choices = c("Under expressed DEG", "Over expressed DEG"),
                                                selected = NULL,
                                                inline = TRUE),
                             ),
                             
                             # Sélection de l'ontologie GO
                             selectInput("inputGO", "Select a GO annotation",
                                         choices = list("Biological process" = "BP",
                                                        "Molecular function" = "MF",
                                                        "Cellular component" = "CC"),
                                         selected = "BP"),
                             
                             # Niveau de précision des termes GO
                             numericInput("goLevel", "GO level of precision | 1 (general) to 7 (ultra-specific)", 
                                          min = 1, max = 7, value = 1),
                             
                             # Bouton pour démarrer l'analyse
                             actionButton("go_analysisButton", "Start analysis", icon = icon('text-background', lib = 'glyphicon')
                          )
                        )
                      )),
                    
                  tabPanel("Results (ORA)", 
                      fluidRow(style = "height: 600px;",
                           box(title = HTML("Bar Plot"),
                               id = "box_go_ora_barplot", width = 6, height = "575px",
                               spinner(plotOutput("go_ora_barplot", height = "575px"))
                               ),
                           box(
                           title = HTML("Dot Plot"),
                           id = "box_go_ora_dotplot", width = 6, height = "575px",
                           spinner(plotOutput("go_ora_dotplot", height = "575px"))
                           )),
                        fluidRow(style = "height: 600px;",
                          box(
                            title = HTML("Gene–Concept Network (microscopic view)"),
                            id = "box_go_ora_netplot", width = 6, height = "575px",
                            spinner(plotOutput("go_ora_netplot", height = "575px"))
                          ),
                          box(
                            title = HTML("Enrichment Map (macroscopic overview)"),
                            id = "box_go_ora_emapplot", width = 6, height = "575px",
                            spinner(plotOutput("go_ora_emapplot", height = "575px"))
                          )),
                      fluidRow(style = "height: 600px;",
                        box(title = HTML("Upset Plot"),
                            id = "box_go_ora_upsetplot", width = 12, height = "575px",
                            spinner(plotOutput("go_ora_upsetplot", height = "575px"))
                        )),
                        fluidRow(
                          box(
                            title = HTML("Results Table"),
                            id = "box_go_ora_table", width = 12,
                            spinner(DT::dataTableOutput("go_ora_table"))
                          )
                        ),
                      fluidRow(
                        box(
                          title = HTML("Hierarchical network of enriched GO terms"),
                          id = "box_go_ora_goplot", width = 12,
                          spinner(plotOutput("go_ora_goplot",height = "1200px"))
                        )
                      )
                                ),
                                
                         tabPanel("Results (GSEA)",
                                  fluidRow(style = "height: 600px;",
                                    box(
                                      title = HTML("Ridge Plot"),
                                      id = "box_go_gsea_ridgeplot", width = 6, height = "575px",
                                      spinner(plotOutput("go_gsea_ridgeplot", height = "575px"))
                                    ),
                                    box(
                                      title = HTML("Dot Plot"),
                                      id = "box_go_gsea_dotplot", width = 6, height = "575px",
                                      spinner(plotOutput("go_gsea_dotplot", height = "575px"))
                                    )
                                  ),
                                  fluidRow(style = "height: 600px;",
                                    box(
                                      title = HTML("Gene–Concept Network (microscopic view)"),
                                      id = "box_go_gsea_netplot", width = 6, height = "575px",
                                      spinner(plotOutput("go_gsea_netplot", height = "575px"))
                                    ),
                                    box(
                                      title = HTML("Enrichment Map (macroscopic overview)"),
                                      id = "box_go_gsea_emapPlot", width = 6, height = "575px",
                                      spinner(plotOutput("go_gsea_emapPlot", height = "575px"))
                                    )),
                                  fluidRow(style = "height: 600px;",
                                    box(title = HTML("Upset Plot"),
                                       id = "box_go_gsea_upsetplot", width = 12, height = "575px",
                                       spinner(plotOutput("go_gsea_upsetplot", height = "575px"))
                                  )),
                                  fluidRow(
                                    box(
                                      title = HTML("GSEA Enrichment Plot"),
                                      id = "box_go_gsea_enrichplot", width = 12,
                                      selectInput("go_selected_gsea", "Select GO terms:", 
                                                  choices = c("None"), 
                                                  selected = "None", 
                                                  multiple = TRUE, 
                                                  width = "100%"),
                                      spinner(plotOutput("go_gsea_enrichplot", height = "425px"))
                                    )),
                                  fluidRow(
                                    box(
                                      title = HTML("Results Table"),
                                      id = "box_go_gsea_table", width = 12,
                                      spinner(DT::dataTableOutput("go_gsea_table"))
                                    )
                                  )
                        ),
                    )),
                  
                  
                  
                  # Pathways Enrichment tab
                  tabItem(tabName = "PATH", h2("Pathway enrichment"),
                          tabsetPanel(
                                     id = "pathway_tabs",
                                     tabPanel("Analysis",
                                     fluidRow (
                                       box(
                                         title = HTML("<i>Analysis parameters</i>"),
                                         id = "GO_analysis", width = 12,
                                         numericInput("pvalPathway", "Select an adjusted p-value cutoff", min = 0, max = 1, value = 0.05, step=0.001),
                                         checkboxGroupInput("analysisMethodChoice", "Analysis method", choices = list("Over Representation Analysis (ORA)"="ORA", "Gene Set Enrichment Analysis (GSEA)"="GSEA"), selected = NULL,
                                                      inline = TRUE, width = NULL),
                                         conditionalPanel(
                                           condition = "input.analysisMethodChoice.includes('ORA')",
                                           sliderInput("pathway_ora_fc", label="ORA - Minimum Absolute Fold-Change (×) to be considered over- or under-enriched", value=1, step=0.01, min=1, max=4),
                                          checkboxGroupInput("pathway_ora_choice", "ORA - Interest", choices = list("Under expressed DEG"="down", "Over expressed DEG"="up"), selected = NULL,
                                                            inline = TRUE, width = NULL)),
                                          
                                         radioButtons("dbPathwaychoice", "Database choice", choices = c("Reactome", "KEGG"), selected = NULL,
                                                      inline = TRUE, width = NULL),
                                         actionButton("pathway_analysis_button", "Start analysis", icon=icon('transfer', lib='glyphicon'))
                                       ),
                                     )),
                                     tabPanel("Results (ORA)",
                                              fluidRow (style = "height: 600px;",
                                                box(title = "Tree plot", width = 7, height="575px",
                                                      spinner(plotOutput("pathway_ora_treeplot", height="575px"))),
                                                box(title = "Dot plot", width = 5, height="575px",
                                                    spinner(plotOutput("pathway_ora_dotplot", height="575px")))
                                                ),
                                                fluidRow(style = "height: 600px;",
                                                         box(title = "Gene–Concept Network (microscopic view)", width = 6, height="575px",
                                                             spinner(plotOutput("pathway_ora_cnetplot", height="575px"))),
                                                          box(title = "Enrichment Map (macroscopic overview)", width = 6, height="575px",
                                                              spinner(plotOutput("pathway_ora_emapplot", height="575px")))),
                                                fluidRow(style = "height: 600px;",
                                                          box(title = "Upset Plot", width = 12, height="575px",
                                                              spinner(plotOutput("pathway_ora_upsetplot",  height="575px")))),
                                                fluidRow (
                                                  box(title = "Table", width = 12,
                                                      spinner(DT::dataTableOutput("pathway_ora_table"))))
                                     ),
                                     tabPanel("Results (GSEA)",
                                                fluidRow (style = "height: 600px;",
                                                  box(title = "Tree plot", width = 7, height="575px", 
                                                      spinner(plotOutput("pathway_gsea_treeplot", height="575px"))),
                                                  box(title = "Dot plot", width = 5, height="575px",
                                                      spinner(plotOutput("pathway_gsea_dotplot",  height="575px")))
                                                  ),
                                                  fluidRow(style = "height: 600px;",
                                                           box(title = "Gene–Concept Network (microscopic view)", width = 6, height="575px",
                                                               spinner(plotOutput("pathway_gsea_cnetplot",  height="575px"))),
                                                           box(title = "Enrichment Map (macroscopic overview)", width = 6, height="575px",
                                                               spinner(plotOutput("pathway_gsea_emapplot", height="575px")))),
                                                fluidRow(style = "height: 600px;",
                                                       box(title = "Upset Plot", width = 12, height="575px",
                                                           spinner(plotOutput("pathway_gsea_upsetplot",  height="575px")))),
                                                  fluidRow (
                                                    box(title = "EnrichPlot", width = 12, height="850px",
                                                        selectInput("pathway_selected_gsea", "Select one or more pathways:", choices = c("None"), selected = "None", multiple=TRUE),
                                                        spinner(plotOutput("pathway_gsea_enrichplot", height = "425px")))),
                                                  fluidRow (
                                                    box(title = "Table", width = 12,
                                                        spinner(DT::dataTableOutput("pathway_gsea_table"))),
                                                  )
                                     ),
                                     tabPanel("View Pathway",
                                              fluidRow(
                                                box(title = "Pathway visualisation", width = 12,
                                                column(12,
                                                       selectInput("pathway", "Select a pathway:", choices = NULL))),
                                                       spinner(imageOutput("pathway_image", height = "600px")
                                                )
                                              )
                                     ))),
                
                  tabItem(tabName = "EXPORT",
                          fluidRow(
                            box(id = "export",  title="Export results", width = 12,
                                uiOutput("export_options"),      
                                conditionalPanel(
                                  condition = "input.export_choices && ((Array.isArray(input.export_choices) && (input.export_choices.includes('GO GSEA') || input.export_choices.includes('Pathway GSEA'))) || input.export_choices == 'GO GSEA' || input.export_choices == 'Pathway GSEA')",
                                  checkboxInput("include_gsea_plots", "Include all GSEA enrichplots", value = TRUE)
                                )
                                
                            )
                          )
                  ))
              ),
              
              #Webpage title
              title="HEATraN"        
)
