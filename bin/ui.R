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
                  menuItem("GO Term Enrichment", icon = icon('text-background', lib='glyphicon'), tabName = "GO_analysis"),
                  menuItem("Pathways Enrichment", icon = icon('transfer', lib='glyphicon'),  tabName = "PATH"),
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
                             conditionalPanel(
                               condition = "input.go_analysisMethodChoice.includes('Over Representation Analysis (ORA)')",
                             checkboxGroupInput("go_oraChoice", "Interest for ORA method",
                                                choices = c("Under expressed DEG", "Over expressed DEG"),
                                                selected = NULL,
                                                inline = TRUE)),
                             
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
                             actionButton("go_analysisButton", "Start analysis", icon = icon('text-background', lib = 'glyphicon')
                          )
                        )
                      )),
                    
                  tabPanel("Results (ORA)", 
                      fluidRow(
                           box(title = HTML("Bar Plot"),
                               id = "goOraBarplot", width = 6,
                               plotOutput("goBarplot", height = "425px")
                               ),
                               box(
                               title = HTML("Dot Plot"),
                               id = "goOraDotplot", width = 6,
                               plotOutput("goDotplot", height = "425px")
                               )),
                      fluidRow(
                        box(title = HTML("Upset Plot"),
                            id = "goOraUpsetplot", width = 6,
                            plotOutput("goUpsetplot", height = "425px")
                        ),
                        box(
                          title = HTML("Gene-Concept Network"),
                          id = "goOraNetplot", width = 6,
                          plotOutput("goNetplot", height = "425px")
                        )),
                      
                        fluidRow(
                          box(
                            title = HTML("Results Table"),
                            id = "goOraTable", width = 12,
                            DT::dataTableOutput("goTable")
                          )
                        ),
                      fluidRow(
                        box(
                          title = HTML("<h4><b>Hierarchical network of enriched GO terms</b></h4>"),
                          id = "goOrasigOfnodesplot", width = 12,
                          plotOutput("gosigOfnodesplot",height = "1200px")
                        )
                      )
                                ),
                                
                                # Onglet GSEA
                         tabPanel("Results (GSEA)",
                                  # Graphiques GSEA dans des box séparées
                                  fluidRow(
                                    box(
                                      title = HTML("Ridge Plot"),
                                      id = "goGseaRidgeplotBox", width = 6,
                                      plotOutput("goGseaRidgeplot", height = "425px")
                                    ),
                                    box(
                                      title = HTML("Dot Plot"),
                                      id = "goGseaDotplotBox", width = 6,
                                      plotOutput("goGseaDotplot", height = "425px")
                                    )
                                  ),
                                  fluidRow(
                                    box(title = HTML("Upset Plot"),
                                        id = "goGseaUpsetplot", width = 6,
                                        plotOutput("goGseaUpsetplot", height = "425px")
                                    ),
                                    box(
                                      title = HTML("Gene-Concept Network"),
                                      id = "goGseaNetplot", width = 6,
                                      plotOutput("goGseaNetplot", height = "425px")
                                    )),
                                  # Box conditionnelle pour les graphiques GSEA détaillés
                                  fluidRow(
                                    box(
                                      title = HTML("GSEA Enrichment Plot"),
                                      id = "goGseaEnrichmentBox", width = 12,
                                      selectInput("goGSEA", "Select GO terms:", 
                                                  choices = c("None"), 
                                                  selected = "None", 
                                                  multiple = TRUE, 
                                                  width = "100%"),
                                      plotOutput("goGseaPlot2", 
                                                 height = "425px",
                                                 click = "plot_click",
                                                 dblclick = "plot_dblclick",
                                                 hover = "plot_hover",
                                                 brush = brushOpts(id = "plot_brush", 
                                                                   delay = 3000, 
                                                                   delayType = "debounce", 
                                                                   fill = "#7e3535", 
                                                                   stroke = "#7e3535"))
                                    )),
                                  fluidRow(
                                    box(
                                      title = HTML("Results Table"),
                                      id = "goGseaTable", width = 12,
                                      DT::dataTableOutput("goGseaTable")
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
                                         sliderInput("pvalPathway", "Select an adjusted p-value cutoff", min = 0, max = 1, value = 0.05, round=F),
                                         checkboxGroupInput("analysisMethodChoice", "Analysis method", choices = list("Over Representation Analysis (ORA)"="ORA", "Gene Set Enrichment Analysis (GSEA)"="GSEA"), selected = NULL,
                                                      inline = TRUE, width = NULL),
                                         conditionalPanel(
                                           condition = "input.analysisMethodChoice.includes('ORA')",
                                          checkboxGroupInput("oraChoice", "Interest for ORA method", choices = list("Under expressed DEG"="down", "Over expressed DEG"="up"), selected = NULL,
                                                            inline = TRUE, width = NULL)),
                                         radioButtons("dbPathwaychoice", "Database choice", choices = c("Reactome", "KEGG"), selected = NULL,
                                                      inline = TRUE, width = NULL),
                                         actionButton("analysisPathwayButton", "Start analysis", icon=icon('transfer', lib='glyphicon'))
                                       ),
                                     )),
                                     tabPanel("Results (ORA)",
                                              # Plots spécifiques à ORA
                                              fluidRow (
                                                box(title = "Tree plot", width = 6, height="850px",
                                                    plotOutput("pathway_ora_treeplot")),
                                                box(title = "Dot plot", width = 6, height="850px",
                                                    plotOutput("pathway_ora_dotplot")),
                                                fluidRow (height="875px",
                                                          box(title = "Upset Plot", width = 6, height="850px",
                                                              plotOutput("pathway_ora_upsetplot")),
                                                          box(title = "Gene-Concept Network", width = 6, height="850px",
                                                              plotOutput("pathway_ora_cnetplot")),
                                                ), 
                                                fluidRow (
                                                  box(title = "Table", width = 12, height="850px",
                                                      DT::dataTableOutput("pathway_ora_table")))
                                              )
                                     ),
                                     tabPanel("Results (GSEA)",
                                              # Plots spécifiques à GSEA
                                                fluidRow (
                                                  box(title = "Tree plot", width = 6, height="850px",
                                                      plotOutput("pathway_gsea_treeplot")),
                                                  box(title = "Dot plot", width = 6, height="850px",
                                                    plotOutput("pathway_gsea_dotplot")),
                                                  fluidRow (height="875px",
                                                    box(title = "Upset Plot", width = 6, height="850px",
                                                        plotOutput("pathway_gsea_upsetplot")),
                                                    box(title = "Gene-Concept Network", width = 6, height="850px",
                                                        plotOutput("pathway_gsea_cnetplot")),
                                                  ), 
                                                  fluidRow (
                                                    box(title = "EnrichPlot", width = 12, height="850px",
                                                        selectInput("pathwayGSEA", "Select one or more pathways:", choices = c("None"), selected = "None", multiple=TRUE),
                                                        plotOutput("pathway_gsea_enrichplot", height = "425px")),
                                                    box(title = "Table", width = 12, height="850px",
                                                        DT::dataTableOutput("pathway_gsea_table")),
                                                  )
                                                )
                                     ),
                                     tabPanel("View Pathway",
                                              fluidRow(
                                                box(width = 12,
                                                column(12,
                                                       selectInput("pathway", "Select pathway:", choices = NULL),
                                                       imageOutput("pathwayImage", height = "600px")
                                                ))
                                              )
                                     ))),
                
                  tabItem(tabName = "EXPORT",
                                fluidRow(
                                  box(id="export", width = 12,
                                    uiOutput("exportOptions"),
                                    downloadButton("exportReport", "Download report")
                                      )
                                  )
                           ))
              ),
              
              #Webpage title
              title="HEATraN"        
)
