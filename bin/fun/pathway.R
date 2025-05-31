# Developed by DAHER Rayan, NAIT EL DJOUDI Lamia & LESAGE Louison (@thelokj).
# rayan.daher@univ-rouen.fr
# lamia.nait-el-djoudi@univ-rouen.fr
# louison.lesage@univ-rouen.fr
# Students at Rouen Normandy University
# Master of Bioinformatics, class M2.2 BIMS 2026 
# Last updated : 08/05/2025
# HEATraN version 0.3.0
options(clusterProfiler.download.method = "wget")

library(clusterProfiler)
library(ReactomePA)
library(pathview)
library(DOSE)
library(enrichplot)
library(dplyr)

# Global variables definition
#GeneID = read.table("~/Bureau/Fac/Rshiny/HEATraN/data/exemple.csv", sep=";", h=T)$GeneID
#Log2FC = read.table("~/Bureau/Fac/Rshiny/HEATraN/data/exemple.csv", sep=";", h=T)$Log2FC
#organism = "Mus musculus"
#data = as.data.frame(cbind(GeneID, Log2FC))
#data$Log2FC = as.numeric(data$Log2FC)

# Organism dataframe
orgs = data.frame(organism="Homo sapiens", db="org.Hs.eg.db", commonName="human", TLname="hsa")
orgs = rbind(orgs, data.frame(organism="Mus musculus", db="org.Mm.eg.db", commonName="mouse", TLname="mmu"))
orgs = rbind(orgs, data.frame(organism="Arabidopsis thaliana", db="org.At.tair.db", commonName="arabidopsis", TLname="ath"))
orgs = rbind(orgs, data.frame(organism="Bos taurus", db="org.Bt.eg.db", commonName="bovine", TLname="bta"))
orgs = rbind(orgs, data.frame(organism="Canis lupus familiaris", db="org.Cf.eg.db", commonName="canine", TLname="cfa"))
orgs = rbind(orgs, data.frame(organism="Gallus gallus", db="org.Gg.eg.db", commonName="chicken", TLname="gga"))
orgs = rbind(orgs, data.frame(organism="Escherichia coli (strain K12)", db="org.EcK12.eg.db", commonName="ecolik12", TLname="eco"))
orgs = rbind(orgs, data.frame(organism="Drosophila melanogaster", db="org.Dm.eg.db", commonName="fly", TLname="dme"))
orgs = rbind(orgs, data.frame(organism="Sus scrofa", db="org.Ss.eg.db", commonName="pig", TLname="ssc"))
orgs = rbind(orgs, data.frame(organism="Rattus norvegicus", db="org.Rn.eg.db", commonName="rat", TLname="rno"))
orgs = rbind(orgs, data.frame(organism="Caenorhabditis elegans", db="org.Ce.eg.db", commonName="celegans", TLname="cel"))
orgs = rbind(orgs, data.frame(organism="Xenopus laevis", db="org.Xl.eg.db", commonName="xenopus", TLname="xla"))
orgs = rbind(orgs, data.frame(organism="Saccharomyces cerevisiae", db="org.Sc.sgd.db", commonName="yeast", TLname="sce"))
orgs = rbind(orgs, data.frame(organism="Danio rerio", db="org.Dr.eg.db", commonName="zebrafish", TLname="dre"))

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
    
    # Générer l'image avec pathview
    tryCatch({
      pathview(gene.data = geneList, pathway.id = pathwayID, species = organism, kegg.dir=".")
      
      base_filename <- paste(pathwayID, ".pathview", sep="")
      png_file <- paste0(base_filename, ".png")
      multi_png_file <- paste0(base_filename, ".multi.png")
      
      generated_files <- c()
      if (file.exists(png_file)) {
        generated_files <- c(generated_files, file.path(wd, "out", png_file))
      } else if (file.exists(multi_png_file)) {
        generated_files <- c(generated_files, file.path(wd, "out", multi_png_file))
      }
      
      setwd(wd)
      return(generated_files)
      
    }, error = function(e) {
      setwd(wd)
      message(paste("Erreur lors de la génération du pathway:", e$message))
      return(NULL)
    })
    
  } else {
    browseKEGG(geneList, pathwayID)
    return(NULL)  # Pas de fichier généré en mode navigateur
  }
}

    
getReactomePathway <- function(rankedLog2FC, pathwayIDs, pathwayDesc, organism) {
  if (!dir.exists("./out/")) {
    dir.create("./out/", recursive = TRUE)
  }
  
  pathway_images <- list()
  
  for (i in seq_along(pathwayIDs)) {
    pathway_id   <- pathwayIDs[i]
    pathway_name <- pathwayDesc[i]
    
    safe_name <- gsub("[^a-zA-Z0-9]", "_", pathway_name)
    file_name <- paste0("./out/", safe_name, ".png")
    
    message(sprintf("Processing pathway: %s (%s)", pathway_name, pathway_id))
    
    # Si le fichier existe déjà, on le réutilise et on passe au suivant
    if (file.exists(file_name)) {
      message(sprintf("File already exists, skipping download: %s", file_name))
      pathway_images[[i]] <- file_name
      next
    }
    
    # Tentative de téléchargement via l'API Reactome
    tryCatch({
      api_url <- paste0(
        "https://reactome.org/ContentService/exporter/diagram/",
        pathway_id, ".png?quality=7"
      )
      
      download.file(api_url, file_name, mode = "wb", quiet = TRUE)
      message(sprintf("Saved pathway image to: %s", file_name))
      pathway_images[[i]] <- file_name
      
    }, error = function(e) {
      message(sprintf(
        "Failed to download from API, trying embedded method for %s", 
        pathway_name
      ))
      
      tryCatch({
        p <- viewPathway(
          pathway_id,
          organism = organism,
          readable = TRUE,
          foldChange = rankedLog2FC
        )
        
        ggplot2::ggsave(file_name, p, width = 10, height = 8)
        message(sprintf("Saved pathway image using viewPathway to: %s", file_name))
        pathway_images[[i]] <- file_name
        
      }, error = function(e2) {
        message(sprintf(
          "Failed to generate pathway image for %s: %s",
          pathway_name, e2$message
        ))
        pathway_images[[i]] <- NULL
      })
    })
  }
  
  return(pathway_images)
}


preprocessPathway = function(data, organism, DB){
  entrezIds = bitr(data$GeneID, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = orgs[orgs$organism==organism, "db"], drop=T)
  if (DB == "KEGG"){
    keggOrganism = as.character(search_kegg_organism(organism, by='scientific_name')['kegg_code'])
    keggIds = bitr_kegg(entrezIds$ENTREZID, fromType="ncbi-geneid", toType="kegg", organism=keggOrganism, drop=T)
    colnames(entrezIds) = c("GeneID", "ENTREZID")
    colnames(keggIds) = c("ENTREZID", "KeggIds")
    getids = merge(entrezIds, keggIds)
    data = merge(data, getids)
    return(list(organism=keggOrganism, genes=data, ranked=rankFC(data, "KeggIds")))
  } else if (DB == "Reactome"){
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
    data <- tryCatch({preprocessPathway(data, organism, DB)}, error = function(e) {return(NULL)})
    if (is.null(data)){
      shinyalert("Database issue", text="Unrecognized ids, please check the selected species or convert your ids using an online converter.", type = "error")
      return()
    }
    genesList = names(data$ranked)[abs(data$ranked)]
    #--- Analysis ---#
    if (DB == "KEGG"){
      enrichment = enrichKEGG(genesList, organism=orgs[orgs$organism==organism, "TLname"], universe=data$entrezIds, pvalueCutoff = pval)
      pathway_images <- NULL
      if(nrow(enrichment@result) > 0) {
        pathway_images <- mapply(function(pathway_id) {
          getKEGGpathway(enrichment@result[enrichment@result$ID == pathway_id, "geneID"], 
                         pathway_id, 
                         organism=orgs[orgs$organism==organism, "TLname"], 
                         local=T)
        }, enrichment@result$ID, SIMPLIFY=FALSE)
      }
    }
    if (DB == "Reactome"){
      enrichment = enrichPathway(genesList, organism=orgs[orgs$organism==organism, "commonName"], universe=data$entrezIds, pvalueCutoff = pval)
      pathway_images <- NULL
        if(nrow(enrichment@result) > 0) {
          pathway_images <- getReactomePathway(data$ranked,
                                               enrichment$ID,
                                               enrichment$Description,
                                               organism=orgs[orgs$organism==organism, "commonName"])
        }
    }
    return(list(enrichment=enrichment, processedData = data, parameters=parameters, pathway_images=pathway_images))
  }
  if (analysis == "GSEA"){
    message("Doing Gene Set Enrichment Analysis on differentially expressed genes regarding their pathway.")
    #--- Preprocess ---#
    data <- tryCatch({preprocessPathway(data, organism, DB)}, error = function(e) {return(NULL)})
    if (is.null(data)){
      shinyalert("Database issue", text="Unrecognized ids, please check the selected species or convert your ids using an online converter.", type = "error")
      return()
    }
    #--- Analysis ---#
    if (DB == "KEGG"){
      enrichment = gseKEGG(data$ranked, organism=orgs[orgs$organism==organism, "TLname"], 
                                               pvalueCutoff = pval,
                                               minGSSize    = 10, #Previously 120
                                               verbose      = FALSE)
      pathway_images <- NULL
      if(nrow(enrichment@result) > 0) {
        pathway_images <- mapply(function(pathway_id) {
          getKEGGpathway(enrichment@result[enrichment@result$ID == pathway_id, "core_enrichment"], 
                         pathway_id, 
                         organism=orgs[orgs$organism==organism, "TLname"], 
                         local=T)
        }, enrichment@result$ID, SIMPLIFY=FALSE)
      }
    }
    if (DB == "Reactome"){
      enrichment = gsePathway(data$ranked, 
                               organism=orgs[orgs$organism==organism, "commonName"],
                               pvalueCutoff = pval,
                               minGSSize    = 10, #Previously 120
                               verbose      = FALSE)
      pathway_images <- NULL
        if(nrow(enrichment@result) > 0) {
          pathway_images <- getReactomePathway(data$ranked,
                                               enrichment$ID,
                                               enrichment$Description,
                                               organism=orgs[orgs$organism==organism, "commonName"])
        }
    }
    
    return(list(enrichment=enrichment, processedData = data, parameters=parameters, pathway_images=pathway_images))
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
#labelSize = 3)


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

dir_size <- function(path, recursive = TRUE) {
  stopifnot(is.character(path))
  files <- list.files(path, full.names = T, recursive = recursive)
  vect_size <- sapply(files, function(x) file.size(x))
  size_files <- sum(vect_size)
  size_files
}
