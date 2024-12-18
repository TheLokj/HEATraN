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

orgs = data.frame(organism="Homo sapiens", db="org.Hs.eg.db", commonName="human", TLname="hsa")
orgs = rbind(orgs, data.frame(organism="Mus musculus", db="org.Mm.eg.db", commonName="mouse", TLname="mmu"))

preprocessKEGG = function(data, organism, pathwayMethod){
  if (pathwayMethod == "KEGG"){
    keggOrganism = as.character(search_kegg_organism(organism, by='scientific_name')['kegg_code'])
    entrezIds = bitr(data$geneID, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = orgs[orgs$organism==organism, "db"], drop=F)
    keggIds = bitr_kegg(entrezIds$ENTREZID, fromType="ncbi-geneid", toType="kegg", organism=keggOrganism, drop=F)
    colnames(entrezIds) = c("geneID", "ENTREZID")
    colnames(keggIds) = c("ENTREZID", "KeggIds")
    getids = merge(entrezIds, keggIds)
    data = merge(data, getids)
    return(list(organism=keggOrganism, genes=data))
  }
}

pathway = function(data, organism, DB, analysis=c("ORA", "GSEA"), pAdjustMethod="BH"){
  if (DB=="KEGG"){
    data = preprocessKEGG(data, organism, pathwayMethod)
    if ("ORA" %in% analysis){
      enrichmentORA = enrichKEGG(data$genes$KeggIds, organism=orgs[orgs$organism==organism, "commonName"], pvalueCutoff = 0.05)
      print(barplot(enrichmentORA, showCategory=20))
      print(dotplot(enrichmentORA, showCategory=30) + ggtitle("dotplot for ORA"))
      }
    if ("GSEA" %in% analysis){
      rankedVec = data$genes$log2FC
      names(rankedVec) = data$genes$KeggIds
      rankedVec = sort(rankedVec, decreasing=T)
      rankedVec= rankedVec[-is.na(names(rankedVec))]
      enrichmentGSEA = gseKEGG(rankedVec, organism=orgs[orgs$organism==organism, "TLname"], 
                               pvalueCutoff = 0.05,
                               minGSSize    = 120,
                               verbose      = FALSE)
      return(enrichmentGSEA)
    }
    if (DB=="Reactome"){
    }
  }
}

gsea = pathway(data, organism, DB="KEGG", analysis="GSEA")

gsea
barplot(gsea, showCategory=20)
dotplot(gsea, showCategory=30) + ggtitle("dotplot for GSEA")
