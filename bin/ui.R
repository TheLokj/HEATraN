# Developed by DAHER Rayan, NAIT EL DJOUDI Lamia & LESAGE Louison (@thelokj).
# rayan.daher@univ-rouen.fr
# lamia.nait-el-djoudi@univ-rouen.fr
# louison.lesage@univ-rouen.fr
# Students at Rouen Normandy University
# Master of Bioinformatics, class M2.2 BIMS 2026 
# Last updated : 13/05/2025
# HEATraN version 0.3.0

# Theme definition
HEATraN_theme <- create_theme(
  adminlte_color(
    red = "#7e3535"
  )
)

fea_page <- readChar("../www/doc.html", file.info("../www/doc.html")$size)
stat_page <- readChar("../www/stat.html", file.info("../www/stat.html")$size)

# -----------------------------------------
# Dashboard creation
# -----------------------------------------
dashboardPage(skin="red", header <- dashboardHeader(title= HTML("<b style='font-size:26px; color:#ebb233; font-weight:900'>HEATraN</b>")), 
              
              # -----------------------------------------
              # Sidebar elements
              # -----------------------------------------
              sidebar <- dashboardSidebar(
                
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
                fileInput("tableInput", 
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
                  menuItem("GO Term Enrichment", icon = icon('text-background', lib='glyphicon'), tabName = "GO",
                    menuSubItem("Analysis", tabName = "GO_analysis", icon = icon('flash', lib='glyphicon')),
                    menuSubItem("Results", tabName = "GO_results", icon = icon('search', lib='glyphicon'))),
                  menuItem("Pathways Enrichment", icon = icon('transfer', lib='glyphicon'),  tabName = "PATH", 
                           menuSubItem("Analysis", tabName = "PATH_analysis", icon = icon('flash', lib='glyphicon')),
                           menuSubItem("Results", tabName = "PATH_results", icon = icon('search', lib='glyphicon'))),
                  menuItem("Export results", icon = icon('share', lib='glyphicon'),  tabName = "EXPORT")
                  )
                ),
              
              # -----------------------------------------
              # Body tabs 
              # -----------------------------------------
              body <- dashboardBody( 
                
                # Load style and theme
                use_theme(HEATraN_theme),
                includeCSS("../www/style.css"),
                
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
                      <em>0.2.0-a.5</em>
                  </div>
              </div>'),
                          
                  ),
                  tabItem(tabName = "DOC",
                          tabsetPanel(id="doctab", 
                                      tabPanel(title="Functional Enrichment Analysis", HTML(fea_page)),
                                      tabPanel(title="Statistics", HTML(stat_page)))
                  ),
                  
                  # Whole Data Inspection tab
                  tabItem(tabName = "WDI",
                          h2("Whole Data Inspection"),
                          
                          fluidRow (
                            box(
                              title = HTML("<b>Volcano plot</b>"),
                              id = "volcanoplot", height = "500px",
                              plotOutput("volcanoPlot", 
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
                              htmlOutput("InfoSelect"),
                              sliderInput("pval", "Maximum ajusted p-value", min = 0, max = 1, value = 0.05, round=F),
                              sliderInput("Log2FC", "Minimum absolute log2(FoldChange)", min = 0, max = 1, value = 0),
                              actionButton("SelectAll", "Select all data", icon=icon('search', lib='glyphicon')),
                              actionButton("ResetButton", "Reset selection", icon=icon('refresh', lib='glyphicon')),
                              actionButton("ZoomButton", "Zoom in the selection", icon=icon('zoom-in', lib='glyphicon')),
                              HTML("<br/><br/><b>Save options</b><br/><br/>"),
                              downloadButton("Download", "Download plot", icon=icon('download-alt', lib='glyphicon')),
                              downloadButton("DownloadTable", "Export table", icon=icon('share', lib='glyphicon'))
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
                    tabName = "GO_analysis",
                    h2("GO Term Enrichment"),
                    
                    # TabsetPanel pour séparer les paramètres et les résultats
                    tabsetPanel(
                      
                      # Onglet 1: Paramètres de l'analyse GO
                      tabPanel("Analysis parameters", tabName = "GO_analysis", 
                               fluidRow(
                                 box(
                                   title = HTML("<i>Analysis parameters</i>"),
                                   id = "gotabset1",
                                   width = 12,
                                   status = NULL,
                                   solidHeader = FALSE,
                                   collapsible = FALSE,
                                   
                                   # Entrée pour la p-value et le q-value
                                   sliderInput("go_pval", "Select an adjusted p-value cutoff", 
                                               min = 0, max = 1, value = 0.05, round = FALSE),
                                   numericInput("qvalueCutoff", "q-value cutoff", 
                                                min = 0, max = 1, value = 0.2, step = 0.01),
                                   
                                   # Sélection de la méthode d'analyse
                                   checkboxGroupInput("go_analysisMethodChoice", "Analysis method",
                                                      choices = c("Over Representation Analysis (ORA)", "Gene Set Enrichment Analysis (GSEA)"),
                                                      selected = NULL,
                                                      inline = TRUE),
                                   
                                   # Sélection des gènes d'intérêt pour ORA
                                   checkboxGroupInput("go_oraChoice", "Interest for ORA method",
                                                      choices = c("Under expressed DEG", "Over expressed DEG"),
                                                      selected = NULL,
                                                      inline = TRUE),
                                   
                                   # Sélection de l'ontologie GO
                                   selectInput("inputGO", "Select a GO annotation",
                                               choices = list("Biological process" = "BP",
                                                              "Molecular function" = "MF",
                                                              "Cellular component" = "CC"),
                                               selected = "BP"),
                                   
                                   # Niveau de précision des termes GO
                                   numericInput("goLevel", "GO level of precision (from 1 to 7)", 
                                                min = 1, max = 7, value = 1),
                                   
                                   # Bouton pour démarrer l'analyse
                                   actionButton("go_analysisButton", "Start analysis", icon = icon('text-background', lib = 'glyphicon'))
                                 )
                               )
                      ),
                      
                      tabPanel(
                        "Results",
                        fluidRow(
                          conditionalPanel(
                            condition = "input.go_analysisMethodChoice.includes('Over Representation Analysis (ORA)')",
                            tabBox(
                              title = "ORA results",
                              width = 12,
                              id = "ora_results_tabs",
                              tabPanel("Up-regulated",
                                       tabsetPanel(
                                         tabPanel("Barplot", plotOutput("goBarplotUp")),
                                         tabPanel("Dotplot", plotOutput("goDotplotUp")),
                                         tabPanel("Network", plotOutput("goNetplotUp")),
                                         tabPanel("Table", DT::dataTableOutput("goTableUp"))
                                       )),
                              tabPanel("Down-regulated",
                                       tabsetPanel(
                                         tabPanel("Barplot", plotOutput("goBarplotDown")),
                                         tabPanel("Dotplot", plotOutput("goDotplotDown")),
                                         tabPanel("Network", plotOutput("goNetplotDown")),
                                         tabPanel("Table", DT::dataTableOutput("goTableDown"))
                                       )),
                              tabPanel("Both-regulated",
                                       tabsetPanel(
                                         tabPanel("Barplot", plotOutput("goBarplotBoth")),
                                         tabPanel("Dotplot", plotOutput("goDotplotBoth")),
                                         tabPanel("Network", plotOutput("goNetplotBoth")),
                                         tabPanel("Table", DT::dataTableOutput("goTableBoth"))
                                       ))
                            )
                          ),
                          conditionalPanel(
                            condition = "input.go_analysisMethodChoice.includes('Gene Set Enrichment Analysis (GSEA)')",
                            tabBox(
                              title = "GSEA results",
                              width = 12,
                              id = "gsea_results_tab",
                              tabPanel("Dotplot", plotOutput("goGseaDotplot")),
                              tabPanel("Enrichment Plot", plotOutput("goGseaEnrichmentPlot")),
                              tabPanel("Ridgeplot", plotOutput("goGseaRidgeplot")),
                              tabPanel("Table", DT::dataTableOutput("goGseaTable"))
                            )
                          )
                        )
                      ),
                      
                    )
                  ),
                  
                  # Pathways Enrichment tab
                  tabItem(tabName = "PATH_analysis", h2("Pathway enrichment"),
                                     fluidRow (
                                       box(
                                         title = HTML("<i>Analysis parameters</i>"),
                                         id = "GO_analysis", width = 12,
                                         sliderInput("pvalPathway", "Select an adjusted p-value cutoff", min = 0, max = 1, value = 0.05, round=F),
                                         radioButtons("analysisMethodChoice", "Analysis method", choices = list("Over Representation Analysis (ORA)"="ORA", "Gene Set Enrichment Analysis (GSEA)"="GSEA"), selected = NULL,
                                                      inline = TRUE, width = NULL),
                                         conditionalPanel(
                                           condition = "input.analysisMethodChoice == 'ORA'",
                                          checkboxGroupInput("oraChoice", "Interest for ORA method", choices = list("Under expressed DEG"="down", "Over expressed DEG"="up"), selected = NULL,
                                                            inline = TRUE, width = NULL)),
                                         radioButtons("dbPathwaychoice", "Database choice", choices = c("Reactome", "KEGG"), selected = NULL,
                                                      inline = TRUE, width = NULL),
                                         actionButton("analysisPathwayButton", "Start analysis", icon=icon('transfer', lib='glyphicon'))
                                       ),
                                     )),
                  tabItem(tabName = "PATH_results", 
                          tabsetPanel(id="pathwaytab1", tabPanel(
                          title="Global results",
                          fluidRow (
                            box(
                            title = HTML("<b>Analysis summary</b>"),
                            id = "analysisDescBox", width = 12,
                            htmlOutput("analysisDesc"))
                          ),fluidRow (
                          box(
                            id = "pathwaytablebox", width = 12,
                            sliderInput("qval", "Select an adjusted p-value to reduce paths", min = 0, max = 1, value = 0.05, round=F),
                            DT::DTOutput('pathwaytable')))),
                          tabPanel(title="Results visualisation",
                                   uiOutput("conditional_gsea_row"),
                                   fluidRow (
                                     box(
                                       height="850px",
                                       title = HTML("<b>Dot plot</b>"),
                                       id = "pathwayplot1",
                                       shinycssloaders::withSpinner(plotOutput("pathwayplotout", 
                                                  height="800px",
                                                  click = "plot_click",
                                                  dblclick = "plot_dblclick",
                                                  hover = "plot_hover",
                                                  brush = brushOpts(id = "plot_brush", delay = 3000, delayType = "debounce", fill="#7e3535", stroke="#7e3535")), type = 6, color = "#ef940b")
                                     ),
                                     box(
                                       title = HTML("<b>Tree plot</b>"),
                                       id = "pathwayplot2",
                                       shinycssloaders::withSpinner(plotOutput("pathwayplotout2", 
                                                  height="800px",
                                                  click = "plot_click",
                                                  dblclick = "plot_dblclick",
                                                  hover = "plot_hover",
                                                  brush = brushOpts(id = "plot_brush", delay = 3000, delayType = "debounce", fill="#7e3535", stroke="#7e3535")), type = 6, color = "#ef940b")
                                     ),
                                   ),
                                   fluidRow (
                                     box(
                                       title = HTML("<b>Gene-Concept Network</b>"),
                                       id = "pathwayplot3", width = 12,
                                       shinycssloaders::withSpinner(plotOutput("pathwayplotout3", 
                                                  height="850px",
                                                  click = "plot_click",
                                                  dblclick = "plot_dblclick",
                                                  hover = "plot_hover",
                                                  brush = brushOpts(id = "plot_brush", delay = 3000, delayType = "debounce", fill="#7e3535", stroke="#7e3535")), type = 6, color = "#ef940b")
                                     ),
                                   )
                                     ),
                            tabPanel(title="Pathway exploration",
                                     box(id = "pathwayselection", height = "100px", width = 12,
                                     selectInput("pathway", "Select a pathway:", choices="None")),
                                     box(width = 12, uiOutput("pathway"), height = "800px")),
                                  )
                                ),
                
                  tabItem(tabName = "EXPORT",
                                fluidRow(
                                  box(id="export", width = 12,
                                    uiOutput("exportOptions"),
                                    downloadButton("exportReport",    
                                                   "Download report")
                                      )
                                  )
                           ))
              ),
              
              #Webpage title
              title="HEATraN"        
)
