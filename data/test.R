#install.packages("wordcloud")
library(clusterProfiler)
library(wordcloud)


library(clusterProfiler)
library(patchwork)
library(igraph)
library(ggraph)

gc()
organism = "org.Mm.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
#reading in input from deseq2
setwd("/home/naitelam/Alternance_1/git_tuto/HEATraN/data")
df = read.csv("/home/naitelam/exemple.csv", header=TRUE,sep = ";")
head(df)
original_gene_list <- df$log2FC

names(original_gene_list) <- df$ID

gene_list<-na.omit(original_gene_list)
#sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
#Exctract significant results (padj < 0.05)
sig_genes_df = subset(df, padj < 0.05)
#From significant results, we want to filter on log2fold change
genes <- sig_genes_df$log2FC
#Name the vector
names(genes) <- sig_genes_df$ID
#omit NA values
genes <- na.omit(genes)
head(genes)
# 
# filter on min log2fold change (log2FoldChange > 2)
# 
# genes <- names(genes)[abs(genes) > 2]
# length(genes)

#Ontology Options: [“BP”, “MF”, “CC”]

organism <- "org.Mm.eg.db"
genes <- names(genes)
universe <- names(gene_list)
go_enrich <- enrichGO(
  gene = genes,
  universe = universe,
  OrgDb = organism,
  keyType = "ENSEMBL",
  readable = TRUE,
  ont = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.10
)


go_enrich <- enrichGO(gene = genes,
                      universe = names(gene_list),
                      OrgDb = organism,
                      keyType = 'ENSEMBL',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10)

head(gene_list)
View(as.data.frame(go_enrich))
#BiocManager::install("enrichplot")
library(enrichplot)
upsetplot(go_enrich)

wcdf<-read.table(text=go_enrich$GeneRatio, sep = "/")[1]
wcdf$term<-go_enrich[,2]
wordcloud(words = wcdf$term, freq = wcdf$V1, scale=(c(4, .1)), colors=brewer.pal(8, "Dark2"), max.words = 25)
wcdf$term <- substr(wcdf$term, 1, 25)  # Truncate terms to 25 characters

wordcloud(
  words = wcdf$term,
  freq = wcdf$V1,
  scale = c(3, 0.5),
  colors = brewer.pal(8, "Dark2"),
  max.words = 25
)
wordcloud(
  words = wcdf$term,
  freq = wcdf$V1,
  scale = c(4, 0.1),
  colors = brewer.pal(8, "Dark2"),
  max.words = 15
)
png("wordcloud.png", width = 1200, height = 1200)
wordcloud(
  words = wcdf$term,
  freq = wcdf$V1,
  scale = c(4, 0.1),
  colors = brewer.pal(8, "Dark2"),
  max.words = 25
)
dev.off()

#Barplot

barplot(go_enrich,
        drop = TRUE,
        showCategory = 10,
        title = "GO Biological Pathways",
        font.size = 8)

dotplot(go_enrich)

#Encrichment map:

emapplot(go_enrich)

go_enrich <- pairwise_termsim(go_enrich)
emapplot(go_enrich, layout.params = list(layout = "kk"), showCategory = 15)

#Enriched GO induced graph:

      
      up_genes <- names(genes)[genes > 0]
      down_genes <- names(genes)[genes < 0]
      df <- as.data.frame(go_enrich)
      head(df,n=1)
      treemap::treemap(df, index = "Description", vSize = "Count", vColor = "p.adjust", 
                       type = "value", palette = "RdYlBu")
    # Correction :
    df <- as.data.frame(goResults()$up)
    if (nrow(df) == 0) {
      plot.new(); title("Aucun terme GO significatif")
      return()
    }
    treemap::treemap(df, index = "Description", vSize = "Count", vColor = "p.adjust", 
                     type = "value", palette = "RdYlBu")
    df_test <- df[1:50, ]  # 10 premières lignes
    df_test$Count <- as.numeric(as.character(df_test$Count))
    df_test$p.adjust <- as.numeric(as.character(df_test$p.adjust))
    
    treemap::treemap(df_test,
                     index = "Description",
                     vSize = "Count",
                     vColor = "p.adjust",
                     type = "value",
                     palette = "RdYlBu")
    library(dplyr)
    
    df_clean <- df %>%
      group_by(Description) %>%
      summarise(
        Count = sum(Count, na.rm = TRUE),
        p.adjust = min(p.adjust, na.rm = TRUE)  # ou mean(p.adjust)
      ) %>%
      ungroup()
    dt=df_clean[1:8,]
      treemap::treemap(dt,
                       index = "Description",
                       vSize = "Count",
                       vColor = "p.adjust",
                       type = "value",
                       palette = "RdYlBu")
    
      dt <- df_clean[1:8,]
      
      # Abbreviation des descriptions longues
      dt$Description <- ifelse(nchar(dt$Description) > 10,
                               paste0(substr(dt$Description, 1, 5), "..."),
                               dt$Description)
      
      # Treemap
      treemap::treemap(dt,
                       index = "Description",
                       vSize = "Count",
                       vColor = "p.adjust",
                       type = "value",
                       palette = "RdYlBu")
      
A=simplify(go_enrich, cutoff = 0.1, by = "p.adjust", select_fun = min)
A 
  
# Cnetplot (réseau gènes-termes)
cnetplot(go_enrich_simple, showCategory = 10, foldChange = gene_list)

# Treeplot
treeplot(go_enrich)



  #########################"
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
  
  
  
  
  
  ##########################################################################################"
  
  BiocManager::install("GOSemSim", force = T)
  BiocManager::install("topGO")
  library(clusterProfiler)
  library(GOSemSim)
  library(topGO)
  # Calculer les niveaux hiérarchiques des termes GO
  
  go_df$GO_level <- godata(organism, ont="BP")@level[go_df$ID]
  
  # Supprimer les NAs et garder GO_level ≤ 4
  go_df_filtered <- go_df[!is.na(go_df$GO_level) & go_df$GO_level <= 4, ]
  
  # Créer un nouvel enrichResult filtré
  go_enrich_filtered <- go_enrich
  go_enrich_filtered@result <- go_df_filtered
 
  plotGOgraph(go_enrich,firstSigNodes=2)
  
  plotGOgraph(go_enrich,firstSigNodes=5)
  showSigOfNodes(go_enrich,  firstSigNodes = 5, useInfo = 'all')
  # printGraph(GOdata, resultFisher, firstSigNodes = 5, fn.prefix = "sampleFile", useInfo = "all", pdfSW = TRUE)
  # ## End(Not run)
  library(topGO)
  GOdata <- new("topGOdata",
                ontology = "BP",
                allGenes = gene_list,
                geneSelectionFun = function(x)(x < 0.05),
                annot = annFUN.org,
                mapping = "org.Mm.eg.db",
                ID = "ensembl")
  
  
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 5)
  
  showSigOfNodes(GOdata, firstSigNodes = 5)
  
  
   printGraph(GOdata)
   # 
   # plotEnrich(go_enrich, plot_type = "gomap", wrap_length = 25,
   #            up_color = '#a32a31',down_color = '#3665a6')
   ##########################""""
 
  
   