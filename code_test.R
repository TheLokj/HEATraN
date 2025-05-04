#install.packages("wordcloud")
library(clusterProfiler)
library(wordcloud)
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


'''

emapplot(go_enrich)
# Erreur dans has_pairsim(x) :
#   Term similarity matrix not available. Please use pairwise_termsim function to deal with the results of enrichment analysis.
# 
    go_enrich <- pairwise_termsim(go_enrich)
    emapplot(go_enrich)
    emapplot(go_enrich, layout = "kk", showCategory = 15)

#Enriched GO induced graph:

goplot(go_enrich, showCategory = 10)

#Category Netplot
categorySize can be either 'pvalue' or 'geneNum'

cnetplot(go_enrich, categorySize="pvalue", foldChange=gene_list)
cnetplot(
go_enrich,
categorySize = "pvalue",
color.params = list(foldChange = gene_list),
layout = "kk"
)

#KEGG Pathway Enrichment

