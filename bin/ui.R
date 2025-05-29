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
                  
                  
                  tabItem(
                    tabName = "GO",
                    h2("GO Term Enrichment"),
                    
                    # TabsetPanel pour séparer les paramètres et les résultats
                    tabsetPanel(
                      
                      # Onglet 1: Paramètres de l'analyse GO
                      tabPanel(
                        "Paramètres de l'analyse GO",
                        fluidRow(
                          box(
                            title = HTML("<i>Analysis parameters</i>"),
                            id = "boxPathwayparamers",
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
                            
                            # Niveau de pour le graph GO
                            selectInput("groupGO", "Gene group", choices = c("up", "down", "both")),
                                       numericInput("nbTopGO", "Nombre de noeuds à afficher", value = 10, min = 1),
                                     
                            # Bouton pour démarrer l'analyse
                            actionButton("go_analysisButton", "Start analysis", icon = icon('text-background', lib = 'glyphicon'))
                          )
                        )
                      ),
                      
                      # Onglet 2: Résultats Analyse GO - séparer Up-regulated, Down-regulated et Both
                      tabPanel(
                        "Résultats Analyse GO",
                        fluidRow(
                          tabBox(
                            title = "Résultats Analyse GO",
                            width = 12,
                            id = "go_results_tabs",
                            
                            # Onglet Up-regulated
                            tabPanel("Up-regulated",
                                     fluidRow(
                                       column(
                                         width = 6,
                                         h4("Up-regulated ORA results"),
                                         tabsetPanel(
                                           tabPanel("Barplot", plotOutput("goBarplotUp")),
                                           tabPanel("Dotplot", plotOutput("goDotplotUp")),
                                           tabPanel("Network", plotOutput("goNetplotUp")),
                                           tabPanel("Cneplot", plotOutput("cneUp")),
                                          
                                           tabPanel("Graphe", plotOutput("plotTopGOUp", height = "600px", width = "100%")),
                                           
                                           tabPanel("UpsetPlot", plotOutput("upsetplotUp")),
                                           
                                           tabPanel("Table", DT::dataTableOutput("goTableUp"))
                                         )
                                       ),
                                       column(
                                         width = 6,
                                         h4("Up-regulated GSEA results"),
                                         plotOutput("gseaPlotUp"),  # À implémenter dans server.R
                                         DT::dataTableOutput("gseaTableUp")
                                       )
                                     )
                            ),
                            
                            # Onglet Down-regulated
                            tabPanel("Down-regulated",
                                     fluidRow(
                                       column(
                                         width = 6,
                                         h4("Down-regulated ORA results"),
                                         tabsetPanel(
                                           tabPanel("Barplot", plotOutput("goBarplotDown")),
                                           tabPanel("Dotplot", plotOutput("goDotplotDown")),
                                           tabPanel("Network", plotOutput("goNetplotDown")),
                                           tabPanel("UpsetPlot", plotOutput("upsetplotDown")),
                                           tabPanel("Cneplot", plotOutput("cneDown")),
                                           tabPanel("Graphe", plotOutput("plotTopGODown", height = "600px", width = "100%")),
                                           
                                           tabPanel("Table", DT::dataTableOutput("goTableDown"))
                                         )
                                       ),
                                       column(
                                         width = 6,
                                         h4("Down-regulated GSEA results"),
                                         plotOutput("gseaPlotDown"),  # À implémenter dans server.R
                                         DT::dataTableOutput("gseaTableDown")
                                       )
                                     )
                            ),
                            
                            # Onglet Both-regulated
                            tabPanel("Both-regulated",
                                     fluidRow(
                                       column(
                                         width = 6,
                                         h4("Both-regulated ORA results"),
                                         tabsetPanel(
                                           tabPanel("Barplot", plotOutput("goBarplotBoth")),
                                           tabPanel("Dotplot", plotOutput("goDotplotBoth")),
                                           tabPanel("Network", plotOutput("goNetplotBoth")),
                                           tabPanel("Cneplot", plotOutput("cneBoth")),
                                          
                                           tabPanel("UpsetPlot", plotOutput("upsetplotBoth")),
                                       
                                           tabPanel("Graphe", plotOutput("plotTopGOBoth", height = "600px", width = "100%")),
                                           
                                           tabPanel("Table", DT::dataTableOutput("goTableBoth"))
                                         )
                                       ),
                                       column(
                                         width = 6,
                                         h4("Both-regulated GSEA results"),
                                         plotOutput("gseaPlotBoth"),  # À implémenter dans server.R
                                         DT::dataTableOutput("gseaTableBoth")
                                       )
                                     )
                            )
                          )
                        )
                      )
                    )
                  )
                  
                  ,
                  
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