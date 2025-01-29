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

# Global variables definition
pAdjustMethod = "BH"
organism = "Mus musculus"
pathwayMethod = "KEGG"
geneID = read.table("data/exemple.csv", sep=";", h=T)$ID
log2FC = read.table("data/exemple.csv", sep=";", h=T)$log2FC
data = as.data.frame(cbind(geneID, log2FC))
data$log2FC = as.numeric(data$log2FC)

# Organism dataframe
orgs = data.frame(organism="Homo sapiens", db="org.Hs.eg.db", commonName="human", TLname="hsa")
orgs = rbind(orgs, data.frame(organism="Mus musculus", db="org.Mm.eg.db", commonName="mouse", TLname="mmu"))

# Functions
rankFC = function(data, colId){
  rankedVec = data$log2FC
  names(rankedVec) = data[,colId]
  rankedVec = sort(rankedVec, decreasing=T)
  rankedVec = rankedVec[!is.na(names(rankedVec))]
  return(rankedVec)
}

getKEGGpathway = function(geneList, pathwayID, organism, local=T){
  if (local){
    wd = getwd()
    setwd(paste(wd, "/out", sep=""))
    pathview(gene.data  = geneList, pathway.id = pathwayID, species = organism, kegg.dir=".")
    setwd(wd)
    } else {
    browseKEGG(gsea$enrichment, pathwayID)
  } 
}
    
preprocessPathway = function(data, organism, DB){
  if (DB == "KEGG"){
    keggOrganism = as.character(search_kegg_organism(organism, by='scientific_name')['kegg_code'])
    entrezIds = bitr(data$geneID, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = orgs[orgs$organism==organism, "db"], drop=F)
    keggIds = bitr_kegg(entrezIds$ENTREZID, fromType="ncbi-geneid", toType="kegg", organism=keggOrganism, drop=F)
    colnames(entrezIds) = c("geneID", "ENTREZID")
    colnames(keggIds) = c("ENTREZID", "KeggIds")
    getids = merge(entrezIds, keggIds)
    data = merge(data, getids)
    return(list(organism=keggOrganism, genes=data, ranked=rankFC(data, "KeggIds")))
  } else if (DB == "Reactome"){
    entrezIds = bitr(data$geneID, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = orgs[orgs$organism==organism, "db"], drop=F)
    colnames(entrezIds) = c("geneID", "ENTREZID")
    data = merge(data, entrezIds)
    print(data)
    return(list(organism=organism, genes=data, ranked=rankFC(data, "ENTREZID")))
  }
}

pathway = function(data, organism, DB, analysis=c("ORA", "GSEA"), pAdjustMethod="BH", threshold=1, oraInterest = "both"){
  if ("ORA" %in% analysis){
    message("Doing an over representation analysis on differentially expressed genes regarding their pathway.")
    #--- Preprocess ---#
    if (oraInterest == "up") {
      message(paste("Analyzing DEG with a log2FC >", threshold, ".", sep=""))
      data = data[data$log2FC>threshold,]}
    else if (oraInterest =="down"){
      message(paste("Analyzing DEG with a log2FC <", threshold, ".", sep=""))
      data = data[data$log2FC<threshold,]
    } else {
      message(paste("Analyzing DEG with a log2FC <", threshold, " and >", threshold, ".", sep=""))
      data = data[abs(data$log2FC)>threshold,]}
    data = preprocessPathway(data, organism, DB)
    genesList = names(data$ranked)[abs(data$ranked)]
    #--- Analysis ---#
    if (DB == "KEGG"){
      enrichmentORA = enrichKEGG(genesList, organism=orgs[orgs$organism==organism, "commonName"], universe=data, pvalueCutoff = 0.05)
    }
    if (DB == "Reactome"){
      enrichmentORA = enrichPathway(genesList, organism=orgs[orgs$organism==organism, "commonName"], universe=data, pvalueCutoff = 1)
    }
    #--- Plotting ---#
    print(barplot(enrichmentORA, showCategory=20))
    print(dotplot(enrichmentORA, showCategory=30) + ggtitle("Dotplot for ORA"))
    return(enrichmentORA)
    }
  if ("GSEA" %in% analysis){
    message("Doing Gene Set Enrichment Analysis on differentially expressed genes regarding their pathway.")
    #--- Preprocess ---#
    data = preprocessPathway(data, organism, DB)
    #--- Analysis ---#
    if (DB == "KEGG"){
      enrichmentGSEA = gseKEGG(data$ranked, organism=orgs[orgs$organism==organism, "TLname"], 
                                               pvalueCutoff = 0.05,
                                               minGSSize    = 10, #Previously 120
                                               verbose      = FALSE)}
    
    if (DB == "Reactome"){
      enrichmentGSEA = gsePathway(data$ranked, 
                               organism=orgs[orgs$organism==organism, "commonName"],
                               pvalueCutoff = 1,
                               minGSSize    = 10, #Previously 120
                               verbose      = FALSE)}
    #--- Plotting ---#
    gseaplot = gseaplot(enrichmentGSEA, 1)
    #gseabarplot = barplot(enrichmentGSEA, showCategory=20)
    gseadotplot = dotplot(enrichmentGSEA, showCategory=30) + ggtitle("dotplot for GSEA")
    return(list(enrichment=enrichmentGSEA, plot=list(plot=gseaplot, dotplot=gseadotplot)))
    }
}



gsea = pathway(data, organism, DB="Reactome", analysis="GSEA")

print(gsea$plot$plot)
print(gsea$plot$dotplot)
print(gsea$enrichment)
print(as.data.frame(gsea$enrichment))

#getKEGGpathway(gsea$enrichment$core_enrichment, gsea$enrichment$ID, organism=orgs[orgs$organism==organism, "TLname"], local=T)

ora = pathway(data, organism, DB="Reactome", analysis="ORA")

convertID = function(id){
  
}
