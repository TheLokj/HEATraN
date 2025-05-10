# Developed by DAHER Rayan, NAIT EL DJOUDI Lamia & LESAGE Louison (@thelokj).
# rayan.daher@univ-rouen.fr
# lamia.nait-el-djoudi@univ-rouen.fr
# louison.lesage@univ-rouen.fr
# Students at Rouen Normandy University
# Master of Bioinformatics, class M2.2 BIMS 2026 
# Last updated : 08/05/2025
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
                
                # Menu
                sidebarMenu(
                  menuItem("Documentation", icon = icon('book', lib='glyphicon'), tabName = "DOC"),
                  menuItem("Whole Data Inspection", icon = icon('eye-open', lib='glyphicon'), tabName = "WDI"),
                  menuItem("GO Term Enrichment", icon = icon('text-background', lib='glyphicon'), tabName = "GO"),
                  menuItem("Pathways Enrichment", icon = icon('transfer', lib='glyphicon'),  tabName = "PATH", 
                           menuSubItem("Analysis", tabName = "PATH_analysis", icon = icon('flash', lib='glyphicon')),
                           menuSubItem("Results", tabName = "PATH_results", icon = icon('search', lib='glyphicon')))
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
                          HTML("<b style='font-size:40px; color:#7e3535; font-weight:900'>HEATraN</b>
                                <br/><br/><p><b>HEATraN</b> (litteraly <i><b>H</b>yper-<b>E</b>xpression <b>A</b>nalysis <b>T</b>ool <b>ra</b>mpantly developed in <b>N</b>ormandy</i>) is a bioinformatics analysis tool dedicated to transcriptomic analysis. It was developed as part of a student project in the Bioinformatics Master of Rouen Normandy University under the supervision of Dauchel Hélène, Education Manager.</p>
                               <br/><br/>
                               <br/><img src='logo.png' class='center' width='712' alt='HEATraN logo'>
                               <br/><div style='text-align:center;'><em>HEATraN by DAHER Rayan, NAIT EL DJOUDI Lamia, LESAGE Louison. University of Rouen Normandy. </em></div>
                               <br/><br/><br/><br/><br/><br/><br/><div style='text-align:right;'><em> You can find its last version on its <a style='font-weight: bold;', href='https://github.com/TheLokj/HEATraN'>GitHub</a>. Current version : 0.2.0-a.5. </em></div>
                               "),
                          
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
                              ),
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
                  
                  # GO Term Enrichment tab, WIP
                  tabItem(tabName = "GO",
                          h2("GO Term Enrichment"),
                          fluidRow (
                            box(
                              title = HTML("<i>Analysis parameters</i>"),
                              id = "boxPathwayparamers", width = 12,
                              sliderInput("pval", "Select an adjusted p-value cutoff", min = 0, max = 1, value = 0.05, round=F),
                              checkboxGroupInput("analysisMethodChoiceGO", "Analysis method", choices = c("Over Representation Analysis (ORA)", "Gene Set Enrichment Analysis (GSEA)"), selected = NULL,
                                                 inline = TRUE, width = NULL, choiceNames = NULL, choiceValues = NULL),
                              checkboxGroupInput("oraChoiceGO", "Interest for ORA method", choices = c("Under expressed DEG", "Over expressed DEG"), selected = NULL,
                                                 inline = TRUE, width = NULL, choiceNames = NULL, choiceValues = NULL),
                              selectInput("inputGO","Select a GO annotation",c("Biological process", "Molecular function", "Cellular component")),
                              numericInput("goLevel", "GO level of precision (from 1 to 7)", min = 1, max = 7, value = 1),
                              actionButton("analysisGOButton", "Start analysis", icon=icon('text-background', lib='glyphicon')),
                            ),
                          ),
                  ),
                  
                  # Pathways Enrichment tab, WIP
                  tabItem(tabName = "PATH_analysis", h2("Pathway enrichment"),
                                     fluidRow (
                                       box(
                                         title = HTML("<i>Analysis parameters</i>"),
                                         id = "boxPathwayparamers", width = 12,
                                         sliderInput("pvalPathway", "Select an adjusted p-value cutoff", min = 0, max = 1, value = 0.05, round=F),
                                         radioButtons("analysisMethodChoice", "Analysis method", choices = list("Over Representation Analysis (ORA)"="ORA", "Gene Set Enrichment Analysis (GSEA)"="GSEA"), selected = NULL,
                                                      inline = TRUE, width = NULL),
                                         checkboxGroupInput("oraChoice", "Interest for ORA method", choices = list("Under expressed DEG"="down", "Over expressed DEG"="up"), selected = NULL,
                                                            inline = TRUE, width = NULL),
                                         radioButtons("dbPathwaychoice", "Database choice", choices = c("Reactome", "KEGG"), selected = NULL,
                                                      inline = TRUE, width = NULL),
                                         actionButton("analysisPathwayButton", "Start analysis", icon=icon('transfer', lib='glyphicon')),
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
                                       title = HTML("<b>Dotplot</b>"),
                                       id = "pathwayplot1",
                                       plotOutput("pathwayplotout", 
                                                  height="800px",
                                                  click = "plot_click",
                                                  dblclick = "plot_dblclick",
                                                  hover = "plot_hover",
                                                  brush = brushOpts(id = "plot_brush", delay = 3000, delayType = "debounce", fill="#7e3535", stroke="#7e3535"))
                                     ),
                                     box(
                                       title = HTML("<b>TreePlot</b>"),
                                       id = "pathwayplot2",
                                       plotOutput("pathwayplotout2", 
                                                  height="800px",
                                                  click = "plot_click",
                                                  dblclick = "plot_dblclick",
                                                  hover = "plot_hover",
                                                  brush = brushOpts(id = "plot_brush", delay = 3000, delayType = "debounce", fill="#7e3535", stroke="#7e3535"))
                                     ),
                                   ),
                                   fluidRow (
                                     box(
                                       title = HTML("<b>Cnetplot</b>"),
                                       id = "pathwayplot3", width = 12,
                                       plotOutput("pathwayplotout3", 
                                                  height="850px",
                                                  click = "plot_click",
                                                  dblclick = "plot_dblclick",
                                                  hover = "plot_hover",
                                                  brush = brushOpts(id = "plot_brush", delay = 3000, delayType = "debounce", fill="#7e3535", stroke="#7e3535"))
                                     ),
                                   )
                                     ),
                            tabPanel(title="Pathway exploration",
                                     box(id = "pathwayselection", height = "100px", width = 12,
                                     selectInput("pathway", "Select a pathway:", choices="None")),
                                     box(width = 12, uiOutput("pathway"), height = "800px")),
                                  )
                                )
                
                )),
              
              #Webpage title
              title="HEATraN"        
)
