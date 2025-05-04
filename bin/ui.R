# Developed by LESAGE Louison (@thelokj).
# louison.lesage@univ-rouen.fr
# Student at Rouen Normandy University
# University project 2024-2025
# Last updated : 18/11/2024
# HEATraN version 0.2.0-a.5

# Theme definition
HEATraN_theme <- create_theme(
  adminlte_color(
    red = "#7e3535"
  )
)

# -----------------------------------------
# Dashboard creation
# -----------------------------------------
dashboardPage(skin="red", header <- dashboardHeader(title= HTML("<b style='font-size:26px; color:#ebb233; font-weight:900'>HEATraN</b>")), 
              
              # -----------------------------------------
              # Sidebar elements
              # -----------------------------------------
              sidebar <- dashboardSidebar(
                
                sidebarMenu(
                  menuItem("Home", tabName = "HOME", icon = icon("home"))),
                
                # File importation: csv, tsv, xls and xlsx
                fileInput("tableInput", 
                          "Import data",
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
                # Ajouter des entrées pour l'ontologie et l'organisme
                selectInput("ontology", "Ontology:", choices = c("BP", "MF", "CC")),
                selectInput("organism", "Organism:", choices = c("org.Mm.eg.db", "org.Hs.eg.db", "org.Sc.sgd.db")),
                
                # Menu
                sidebarMenu(
                  menuItem("Whole Data Inspection", icon = icon('eye-open', lib='glyphicon'), tabName = "WDI"),
                  menuItem("GO Term Enrichment", icon = icon('text-background', lib='glyphicon'), tabName = "GO"),
                  menuItem("Pathways Enrichment", icon = icon('transfer', lib='glyphicon'),  tabName = "PATH"),
                  menuItem("About", icon = icon('info-sign', lib='glyphicon'), tabName = "ABOUT")
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
                          h2("Project home"),
                          p("Click on the Whole Data Inspection menu to access to the second model!"),
                  ),
                  
                  # Whole Data Inspection tab
                  tabItem(tabName = "WDI",
                          h2("Whole Data Inspection"),
                          
                          fluidRow (
                            box(
                              title = HTML("<b>Volcano plot</b>"),
                              id = "volcanoplot", 
                              height = "500px",
                              width = 12,  # Add width argument
                              status = "primary",  # Add status argument
                              solidHeader = TRUE,  # Add solidHeader argument
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
                              sliderInput("pval", "Maximum p-value", min = 0, max = 1, value = 0.05, round=F),
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
                  
                  # GO Term Enrichment tab, WIP
                  tabItem(tabName = "GO",
                          h2("GO Term Enrichment"),
                          fluidRow(
                            box(
                              title = "GO Analysis Parameters",
                              width = 4,
                              selectInput("goOntology", "Ontology Type:",
                                          choices = c("Biological Process" = "BP",
                                                      "Molecular Function" = "MF",
                                                      "Cellular Component" = "CC")),
                              numericInput("pvalueCutoff", "P-value Cutoff:", 0.05, min = 0, max = 1, step = 0.01),
                              numericInput("qvalueCutoff", "Q-value Cutoff:", 0.10, min = 0, max = 1, step = 0.01),
                              actionButton("runGO", "Run GO Analysis", icon = icon("play"))
                            ),
                            box(
                              title = "Visualization Options",
                              width = 8,
                              tabsetPanel(
                                tabPanel("Word Cloud (Up-regulated)",
                                         sliderInput("maxWordsUp", "Maximum Words:", min = 5, max = 50, value = 25),
                                         plotOutput("wordcloudPlotUp", height = "400px")
                                ),
                                tabPanel("Word Cloud (Down-regulated)",
                                         sliderInput("maxWordsDown", "Maximum Words:", min = 5, max = 50, value = 25),
                                         plotOutput("wordcloudPlotDown", height = "400px")
                                ),
                                tabPanel("Bar Plot (Up-regulated)",
                                         sliderInput("topCategoriesUp", "Top Categories:", min = 5, max = 30, value = 10),
                                         plotOutput("goBarplotUp", height = "400px")
                                ),
                                tabPanel("Bar Plot (Down-regulated)",
                                         sliderInput("topCategoriesDown", "Top Categories:", min = 5, max = 30, value = 10),
                                         plotOutput("goBarplotDown", height = "400px")
                                ),
                                tabPanel("Dot Plot (Up-regulated)",
                                         plotOutput("goDotplotUp", height = "400px")
                                ),
                                tabPanel("Dot Plot (Down-regulated)",
                                         plotOutput("goDotplotDown", height = "400px")
                                ),
                                tabPanel("Network Plot (Up-regulated)",
                                         plotOutput("goNetplotUp", height = "400px")
                                ),
                                tabPanel("Network Plot (Down-regulated)",
                                         plotOutput("goNetplotDown", height = "400px")
                                )
                              )
                            )
                          ),
                          fluidRow(
                            box(
                              title = "GO Enrichment Results (Up-regulated)",
                              width = 12,
                              DT::DTOutput("goTableUp")
                            ),
                            box(
                              title = "GO Enrichment Results (Down-regulated)",
                              width = 12,
                              DT::DTOutput("goTableDown")
                            )
                          )
                  ),
                  
                  # Pathways Enrichment tab, WIP
                  tabItem(tabName = "PATH",
                          h2("Pathways Enrichment")
                  ),
                  
                  # Information tab
                  tabItem(tabName = "ABOUT",
                          h2("About"),
                          HTML("<p><b>HEATraN</b> (litteraly <i><b>H</b>yper-<b>E</b>xpression <b>A</b>nalysis <b>T</b>ool <b>ra</b>mpantly developed in <b>N</b>ormandy</i>) is a bioinformatics analysis tool dedicated to transcriptomic analysis. It was developed as part of a student project in the Bioinformatics Master of Rouen Normandy University.</p>
                               <br/><img src='logo.png' class='center' width='512' alt='HEATraN logo'>
                               <br/><p>You can find its last version on its <a style='font-weight: bold;', href='https://github.com/TheLokj/HEATraN'>GitHub</a>.<br/></p>
                               <i style='text-align:right'> Current version : 0.2.0-a.5.</i>")
                  )
                )),
              
              #Webpage title
              title="HEATraN"        
)
