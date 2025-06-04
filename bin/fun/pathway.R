# Developed by DAHER Rayan, NAIT EL DJOUDI Lamia & LESAGE Louison (@thelokj).
# rayan.daher@univ-rouen.fr
# lamia.nait-el-djoudi@univ-rouen.fr
# louison.lesage@univ-rouen.fr
# Students at Rouen Normandy University
# Master of Bioinformatics, class M2.2 BIMS 2026 
# Last updated : 04/06/2025
# HEATraN version 1.0.0
options(clusterProfiler.download.method = "wget")

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
rankFC = function(data, col_id){
  ranked_vec = data$Log2FC
  names(ranked_vec) = data[[col_id]]
  ranked_vec = ranked_vec[!is.na(names(ranked_vec))]
  ranked_vec = tapply(ranked_vec, names(ranked_vec), mean)
  ranked_vec = sort(ranked_vec, decreasing=T)
  namesrV = names(ranked_vec)
  ranked_vec = as.numeric(ranked_vec)
  names(ranked_vec) = namesrV
  return(ranked_vec)
}

get_kegg_pathway <- function(gene_list, pathway_id, organism, local = TRUE) {
  out_dir <- file.path(getwd(), "out")
  if (local && !dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  base_name     <- paste0(pathway_id, ".pathview")
  png_file      <- file.path(out_dir, paste0(base_name, ".png"))
  multi_png_file<- file.path(out_dir, paste0(base_name, ".multi.png"))
  
  if (local) {
    existing <- c(png_file, multi_png_file)[file.exists(c(png_file, multi_png_file))]
    if (length(existing) > 0) {
      message("File already exists, skipping download:", paste(existing))
      return(existing)
    }
    
    wd <- getwd()
    setwd(out_dir)
    on.exit(setwd(wd), add = TRUE)
    
    tryCatch({
      pathview(
        gene.data   = gene_list,
        pathway.id  = pathway_id,
        species     = organism,
        kegg.dir    = ".",
        gene.idtype = "KEGG"
      )
      generated <- c()
      if (file.exists(basename(png_file))) {
        generated <- c(generated, png_file)
      } 
      if (file.exists(basename(multi_png_file))) {
        generated <- c(generated, multi_png_file)
      }
      if (length(generated) == 0) {
        warning("Something went wrong with generation")
        return(NULL)
      }
      return(generated)
      
    }, error = function(e) {
      message("Error: ", e$message)
      return(NULL)
    })
    
  } else {
    browseKEGG(gene_list, pathway_id)
    return(NULL)
  }
}

get_reactome_pathway <- function(ranked_log2fc, pathway_ids, pathway_descs, organism) {
  if (!dir.exists("./out/")) {
    dir.create("./out/", recursive = TRUE)
  }
  
  pathway_images <- list()
  
  for (i in seq_along(pathway_ids)) {
    pathway_id   <- pathway_ids[i]
    pathway_name <- pathway_descs[i]
    
    safe_name <- gsub("[^a-zA-Z0-9]", "_", pathway_name)
    file_name <- paste0("./out/", safe_name, ".png")
    
    message(sprintf("Processing pathway: %s (%s)", pathway_name, pathway_id))
    
    if (file.exists(file_name)) {
      message(sprintf("File already exists, skipping download: %s", file_name))
      pathway_images[[i]] <- file_name
      next
    }
    
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
          foldChange = ranked_log2fc
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
  # Convertir TOUS les gènes originaux pour l'univers
  entrezIds = bitr(data$GeneID, fromType = "ENSEMBL", toType = "ENTREZID", 
                   OrgDb = orgs[orgs$organism==organism, "db"], drop=T)
  
  if (DB == "KEGG"){
    keggOrganism = as.character(search_kegg_organism(organism, by='scientific_name')['kegg_code'])
    keggIds = bitr_kegg(entrezIds$ENTREZID, fromType="ncbi-geneid", toType="kegg", 
                        organism=keggOrganism, drop=T)
    
    colnames(entrezIds) = c("GeneID", "ENTREZID")
    colnames(keggIds) = c("ENTREZID", "KeggIds")
    
    # Univers = TOUS les KEGG IDs convertibles depuis les données originales
    universe_kegg = merge(entrezIds, keggIds)
    
    # Données pour l'analyse = merge avec les données originales
    getids = merge(entrezIds, keggIds)
    processed_data = merge(data, getids)
    
    return(list(
      organism = keggOrganism, 
      genes = processed_data, 
      ranked = rankFC(processed_data, "KeggIds"),
      universe = universe_kegg$KeggIds[!is.na(universe_kegg$KeggIds)]  
    ))
    
  } else if (DB == "Reactome"){
    colnames(entrezIds) = c("GeneID", "ENTREZID")
    
    universe_entrez = entrezIds
    
    processed_data = merge(data, entrezIds)
    
    return(list(
      organism = organism, 
      genes = processed_data, 
      ranked = rankFC(processed_data, "ENTREZID"),
      universe = universe_entrez$ENTREZID[!is.na(universe_entrez$ENTREZID)]  # Vecteur d'IDs
    ))
  }
}


pathway = function(data, organism, DB, analysis="GSEA", pAdjustMethod="BH", threshold=0, ora_interest = c("up", "down"), pval= 0.05){
  parameters = list(organism=organism, DB=DB, analysis=analysis, pAdjustMethod=pAdjustMethod, threshold=threshold, ora_interest=ora_interest, pval=pval)
  
  if (analysis == "ORA"){
    message("Doing an over representation analysis on differentially expressed genes regarding their pathway.")
    
    #--- Preprocess ---#
    background_data <- tryCatch({preprocessPathway(data, organism, DB)}, error = function(e) {return(NULL)})
    if (is.null(background_data)){
      shinyalert("Database issue", text="Unrecognized ids, please check the selected species or convert your ids using an online converter.", type = "error")
      return()
    }
    
    sig_genes_df <- as.data.frame(na.omit(background_data$genes[background_data$genes$padj < pval,]))
    #--- Filtration ---#
    if ("up" %in% ora_interest & "down" %in% ora_interest) {
      message(paste("Analyzing DEG with a Log2FC < -", log2(threshold), " and > ", log2(threshold), ".", sep=""))
      filtered_data = sig_genes_df[abs(sig_genes_df$Log2FC) > log2(threshold),]
    } else if (length(ora_interest) == 1 && ora_interest == "down"){
      message(paste("Analyzing DEG with a Log2FC < -", log2(threshold), ".", sep=""))
      filtered_data = sig_genes_df[sig_genes_df$Log2FC < -log2(threshold),]
    } else if (length(ora_interest) == 1 && ora_interest == "up"){
      message(paste("Analyzing DEG with a Log2FC >", log2(threshold), ".", sep=""))
      filtered_data = sig_genes_df[sig_genes_df$Log2FC > log2(threshold),]
    } else {
      return()
    }
    genes_list = if(DB == "KEGG") filtered_data$KeggIds else filtered_data$ENTREZID

    ratio = (length(genes_list) / length(background_data$universe))
    message(paste("ORA: selected genes represent ", round(ratio*100,2), "% of Universe", sep=""))
    if (ratio > 0.1) {
      shinyalert("Statistical issue", text=paste("The genes selected for the ORA analysis represent more than 10% of the initial universe (", round(ratio*100, 2), "%), which can greatly distort and negatively impact the statistical analysis. Please increase the minimum fold change threshold or decrease the maximum p-value threshold.", sep=""), type = "error")
      return()
    }
    #--- Analysis ---#
    if (DB == "KEGG"){
      enrichment = enrichKEGG(genes_list, 
                              organism=orgs[orgs$organism==organism, "TLname"], 
                              universe=background_data$universe, 
                              pvalueCutoff = pval, 
                              pAdjustMethod = pAdjustMethod, 
                              qvalueCutoff = as.numeric(read.ini("./conf.ini")$STAT$q_val))
      setProgress(message="Downloading pathway visualisations related to ORA enrichment...", value = 0.66)
      pathway_images <- NULL
      if ((!is.null(enrichment)) && (!is.null(enrichment@result)) && (nrow(enrichment@result) > 0)) {
        message(paste(nrow(enrichment), "enriched pathways"))
        pathway_images <- mapply(function(pathway_id) {
          get_kegg_pathway(enrichment@result[enrichment$ID == pathway_id, "geneID"], 
                         pathway_id, 
                         organism=orgs[orgs$organism==organism, "TLname"], 
                         local=T)
        }, enrichment$ID, SIMPLIFY=FALSE)
      }
    } else if (DB == "Reactome"){
      enrichment = enrichPathway(genes_list, 
                                 organism=orgs[orgs$organism==organism, "commonName"], 
                                 universe=background_data$universe, 
                                 pvalueCutoff = pval, 
                                 pAdjustMethod = pAdjustMethod,
                                 qvalueCutoff = as.numeric(read.ini("./conf.ini")$STAT$q_val))
      setProgress(message="Downloading pathway visualisations related to ORA enrichment...")
      pathway_images <- NULL
      if ((!is.null(enrichment)) && (!is.null(enrichment@result)) && (nrow(enrichment@result) > 0)) {
        message(paste(nrow(enrichment), "enriched pathways"))
        pathway_images <- get_reactome_pathway(background_data$ranked,
                                             enrichment$ID,
                                             enrichment$Description,
                                             organism=orgs[orgs$organism==organism, "commonName"])
      }
    }
    return(list(enrichment=enrichment, processed_data = background_data, parameters=parameters, pathway_images=pathway_images))}
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
                                               verbose      = FALSE, 
                           pAdjustMethod = pAdjustMethod)
      setProgress(message="Downloading pathway visualisations related to GSEA enrichment...", value = 0.66)
      pathway_images <- NULL
      if ((!is.null(enrichment)) && (!is.null(enrichment@result)) && (nrow(enrichment@result) > 0)) {
        message(paste(nrow(enrichment), "enriched pathways"))
        pathway_images <- mapply(function(pathway_id) {
          get_kegg_pathway(enrichment@result[enrichment$ID == pathway_id, "core_enrichment"], 
                         pathway_id, 
                         organism=orgs[orgs$organism==organism, "TLname"], 
                         local=T)
        }, enrichment$ID, SIMPLIFY=FALSE)
      }
    } else if (DB == "Reactome"){
      enrichment = gsePathway(data$ranked, 
                               organism=orgs[orgs$organism==organism, "commonName"],
                               pvalueCutoff = pval,
                               verbose      = FALSE, 
                              pAdjustMethod = pAdjustMethod)
      setProgress(message="Downloading pathway visualisations related to GSEA enrichment...")
      pathway_images <- NULL
      if ((!is.null(enrichment)) && (!is.null(enrichment@result)) && (nrow(enrichment@result) > 0)) {
        message(paste(nrow(enrichment), "enriched pathways"))  
        pathway_images <- get_reactome_pathway(data$ranked,
                                               enrichment$ID,
                                               enrichment$Description,
                                               organism=orgs[orgs$organism==organism, "commonName"])
        }
    }
    return(list(enrichment=enrichment, processed_data = data, parameters=parameters, pathway_images=pathway_images))
    }
}

dir_size <- function(path, recursive = TRUE) {
  stopifnot(is.character(path))
  files <- list.files(path, full.names = T, recursive = recursive)
  vect_size <- sapply(files, function(x) file.size(x))
  size_files <- sum(vect_size)
  size_files
}
