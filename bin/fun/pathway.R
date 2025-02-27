# Developed by LESAGE Louison (@thelokj).
# louison.lesage@univ-rouen.fr
# Student at Rouen Normandy University
# University project 2024-2025
# Last updated : 18/11/2024
# HEATraN version 0.2.0-a.5

library(clusterProfiler)
library(ReactomePA)
library(pathview)
library(DOSE)
library(enrichplot)
library(dplyr)
setwd("~/Bureau/Fac/Rshiny/HEATraN")

# Global variables definition
GeneID = read.table("~/Bureau/Fac/Rshiny/HEATraN/data/exemple.csv", sep=";", h=T)$GeneID
Log2FC = read.table("~/Bureau/Fac/Rshiny/HEATraN/data/exemple.csv", sep=";", h=T)$Log2FC
organism = "Mus musculus"
data = as.data.frame(cbind(GeneID, Log2FC))
data$Log2FC = as.numeric(data$Log2FC)

# Organism dataframe
orgs = data.frame(organism="Homo sapiens", db="org.Hs.eg.db", commonName="human", TLname="hsa")
orgs = rbind(orgs, data.frame(organism="Mus musculus", db="org.Mm.eg.db", commonName="mouse", TLname="mmu"))

# Functions
rankFC = function(data, colId){
  rankedVec = data$Log2FC
  names(rankedVec) = data[[colId]]
  rankedVec = rankedVec[!is.na(names(rankedVec))]
  rankedVec = tapply(rankedVec, names(rankedVec), mean)
  rankedVec = sort(rankedVec, decreasing=T)
  namesrV = names(rankedVec)
  rankedVec = as.numeric(rankedVec)
  names(rankedVec) = namesrV
  return(rankedVec)
}

getKEGGpathway = function(geneList, pathwayID, organism, local=T){
  if (local){
    wd = getwd()
    setwd(paste(wd, "/out", sep=""))
    pathview(gene.data  = geneList, pathway.id = pathwayID, species = organism, kegg.dir=".")
    setwd(wd)
    } else {
    browseKEGG(geneList, pathwayID)
  } 
}
    
getReactomePathway = function(rankedLog2FC, pathwayDesc, organism){
  print(pathwayDesc)
  print(rankedLog2FC)
  for (i in 1:length(pathwayDesc)){
    name = paste(gsub("[\\\\ ]", "_", pathwayDesc[i]),".png",sep="")
    if (! name %in% list.files("./out/")){
      try(ggplot2::ggsave(paste(getwd(), "/out/", name, sep=""), viewPathway(pathwayDesc[i], 
                                                                          organism="mouse",
                                                                          readable = TRUE, 
                                                                          foldChange = rankedLog2FC[names(table(names(rankedLog2FC))[table(names(rankedLog2FC))==1])])))}
    else {
      message(paste(name, "already saved."))
    }}
  
}

preprocessPathway = function(data, organism, DB){
  if (DB == "KEGG"){
    keggOrganism = as.character(search_kegg_organism(organism, by='scientific_name')['kegg_code'])
    entrezIds = bitr(data$GeneID, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = orgs[orgs$organism==organism, "db"], drop=T)
    keggIds = bitr_kegg(entrezIds$ENTREZID, fromType="ncbi-geneid", toType="kegg", organism=keggOrganism, drop=T)
    colnames(entrezIds) = c("GeneID", "ENTREZID")
    colnames(keggIds) = c("ENTREZID", "KeggIds")
    getids = merge(entrezIds, keggIds)
    data = merge(data, getids)
    return(list(organism=keggOrganism, genes=data, ranked=rankFC(data, "KeggIds")))
  } else if (DB == "Reactome"){
    entrezIds = bitr(data$GeneID, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = orgs[orgs$organism==organism, "db"], drop=T)
    colnames(entrezIds) = c("GeneID", "ENTREZID")
    data = merge(data, entrezIds)
    return(list(organism=organism, genes=data, ranked=rankFC(data, "ENTREZID")))
  }
}

pathway = function(data, organism, DB, analysis="GSEA", pAdjustMethod="BH", threshold=1, oraInterest = c("up", "down"), pval= 0.05){
  parameters = list(organism=organism, DB=DB, analysis=analysis, pAdjustMethod=pAdjustMethod, threshold=threshold, oraInterest=oraInterest, pval=pval)
  if (analysis == "ORA"){
    message("Doing an over representation analysis on differentially expressed genes regarding their pathway.")
    #--- Preprocess ---#
    if ("up" %in% oraInterest & "down" %in% oraInterest) {
      message(paste("Analyzing DEG with a Log2FC <", threshold, " and >", threshold, ".", sep=""))
      data = data[abs(data$Log2FC)>threshold,]}
    else if (oraInterest == "down"){
      message(paste("Analyzing DEG with a Log2FC <", threshold, ".", sep=""))
      data = data[data$Log2FC<threshold,]}
    else if (oraInterest == "up"){
      message(paste("Analyzing DEG with a Log2FC >", threshold, ".", sep=""))
      data = data[data$Log2FC>threshold,]}
    data = preprocessPathway(data, organism, DB)
    genesList = names(data$ranked)[abs(data$ranked)]
    #--- Analysis ---#
    if (DB == "KEGG"){
      enrichment = enrichKEGG(genesList, organism=orgs[orgs$organism==organism, "TLname"], universe=data$entrezIds, pvalueCutoff = pval)
      getKEGGpathway(enrichment@gene, enrichment$ID, organism=orgs[orgs$organism==organism, "TLname"], local=T)
    }
    if (DB == "Reactome"){
      enrichment = enrichPathway(genesList, organism=orgs[orgs$organism==organism, "commonName"], universe=data$entrezIds, pvalueCutoff = pval)
      print(enrichment)
      getReactomePathway(data$ranked, enrichment$Description, organism=organism)
       }
    return(list(enrichment=enrichment, processedData = data, parameters=parameters))
    }
  if (analysis == "GSEA"){
    message("Doing Gene Set Enrichment Analysis on differentially expressed genes regarding their pathway.")
    #--- Preprocess ---#
    data = preprocessPathway(data, organism, DB)
    #--- Analysis ---#
    if (DB == "KEGG"){
      enrichment = gseKEGG(data$ranked, organism=orgs[orgs$organism==organism, "TLname"], 
                                               pvalueCutoff = pval,
                                               minGSSize    = 10, #Previously 120
                                               verbose      = FALSE)
      getKEGGpathway(enrichment$core_enrichment, enrichment$ID, organism=orgs[orgs$organism==organism, "TLname"], local=T)}
    if (DB == "Reactome"){
      enrichment = gsePathway(data$ranked, 
                               organism=orgs[orgs$organism==organism, "commonName"],
                               pvalueCutoff = pval,
                               minGSSize    = 10, #Previously 120
                               verbose      = FALSE)
    getReactomePathway(data$ranked, enrichment$Description, organism=organism)}
    return(list(enrichment=enrichment, processedData = data, parameters=parameters))
    }
}


#test = pathway(data, organism, DB="Reactome")
#
#ggplot(test$enrichment, aes(x = reorder(Description, setSize), y = setSize, fill=p.adjust)) +
# geom_bar(stat = "identity") +
#  coord_flip() +  # Inverser les axes pour une meilleure lisibilité
    # labs(title = "Résultats de l'enrichissement GSEA",
#      x = "Voies",
       #      y = "Set size") +
#  theme_minimal() +
# scale_fill_gradient2(low = "royalblue", mid="royalblue", high = "red3")




#cnetplot(
#  x = test$enrichment, category = "Description",  gene = "core_enrichment", pvalue = "p.adjust",  foldchange = test$processedData$ranked,  # Si vous avez des fold changes, vous pouvez les inclure ici
# pointSize = 3,
# labelSize = 3)


#geneSetsSignificant = test$enrichment@geneSets[names(test$enrichment@geneSets)%in%test$enrichment$ID]

#geneSetsSignificant = lapply(geneSetsSignificant, as.vector)

#gene_sets <- data.frame(geneSetsSignificant) %>%
# select(geneID, Description) %>%
# table() > 0


#print(gsea$plot$plot)
##print(gsea$plot$dotplot)
#print(gsea$enrichment)
#print(as.data.frame(gsea$enrichment))


#Plots = function(id){
  #--- Plotting ---#
  # print(barplot(enrichmentORA, showCategory=20))
  # print(dotplot(enrichmentORA, showCategory=30) + ggtitle("Dotplot for ORA"))
  #--- Plotting ---#
  # gseaplot = gseaplot(enrichmentGSEA, 1)
  # gseabarplot = barplot(enrichmentGSEA, showCategory=20)
  # gseadotplot = dotplot(enrichmentGSEA, showCategory=30) + ggtitle("dotplot for GSEA")
  
#}

