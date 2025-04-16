#===========================#
# Pathway Analysis Functions
#===========================#
# External packages: fgsea, DOSE, ReactomePA, AnnotationHub (for GO terms), biomaRt (for ID mapping)
# External packages: clusterProfiler

# Scope: Perform pathway analyses (Over-representation, GSEA, etc.) on results from RWR_pipeline() or AMEND()

# Inputs: AMENDresult object, RWR_pipeline result (vector or list or matrix),
#   Pathway database (e.g., GO, Reactome, KEGG), molecular IDs, Enrichment method (GSEA or ORA), data for GSEA (seed value, original data value, Z-score)

# Other input arguments: which iteration of AMEND to analyze (best, iteration k, or all)

# If AMENDresult, coerce method to ORA
#   If ORA, apply ORA to one or several subnetworks, depending on 'iter'
# If RWRresult, coerce method to GSEA
#   For GSEA, apply GSEA to one or several numeric vectors
# If matrix, coerce method to GSEA
#   Apply GSEA to each column of matrix independently

# supported Databases: GO, (M)KEGG, Reactome, WikiPathways, DO(Lite), NCG, DisGeNET, snpDisGeNET

# supported organisms for KEGG: http://www.genome.jp/kegg/catalog/org_list.html
# supported organisms for HMNM package (also for ReactomePA package): "human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly"

## keyType parameter...
# supported identifiers for WikiPathways: Entrez IDs
# supported identifiers for Reactome: Entrez IDs
# supported identifiers for DO(Lite), NCG, DisGeNET, snpDisGeNET: Entrez IDs
# supported ids for (M)KEGG: "kegg", "ncbi-geneid", "ncbi-proteinid", "uniprot"
# common supported keytypes for GO: ENTREZID, ENSEMBL, ENSEMBLPROT, ENSEMBLTRANS, REFSEQ, SYMBOL, UNIGENE, UNIPROT... and many others
#   See ?AnnotationDbi::ACCNUM


# *** IDEA *** (1/24/25): Be able to jointly analyze node sets of different molecules (i.e., different keyTypes). Currently unable to.
# This is a critical need, since there are databases like KEGG & Reactome that have annotations containing genes/proteins and metabolites.
# To do this, I cannot rely on clusterProfiler. Perhaps borrow from multiGSEA.
# Resources: 
#   multiGSEA package: https://bioconductor.org/packages/release/bioc/html/multiGSEA.html
#   metaboliteIDmapping package: https://bioconductor.org/packages/release/data/annotation/html/metaboliteIDmapping.html

# Goal: To provide pathway analysis capabilites that are as general as possible-in terms of annotation databases & terms, ID mapping, and enrichment method- but at the same time takes advantage of already-existing packages.
# Potentially Useful Packages: clusterProfiler, multiGSEA, metaboliteIDmapping, fgsea
# Maybe move away from clusterProfiler code and rely instead on multiGSEA, metaboliteIDmapping, graphite, and fgsea....

# IDEA (2/7/25): Allow users to specify the node attribute of the network to use as node IDs, since the function as-is might not be able to perform the mapping the user needs. 
# So the user can perform the necessary mapping prior to pathway_analysis() and then specify the node attribute which stores those node IDs to be used
# Currently the function only accesses 'original_name' vertex attribute.

#=== Preliminary Steps ===#
# 1. Define the acceptable ontologies (e.g., KEGG, GO, Reactome)
# 2. Be able to access these ontologies
# 3. Determine which organism are valid for each ontology
# 4. Define the acceptable organisms (based on these ontologies)
# 5. Determine the keytypes present in each ontology
# 6. Categorize these keytypes into keyclasses
# 7. Choose mapping databases (e.g., AnnotationDbi, metaboliteIDmapping)
# 8. From the chosen mapping databases, determine acceptable keytypes
# 9. Categorize these keytypes into keyclasses

#=== Steps for Pathway Analysis ===#
## [User]
# 1. Define layer sets to analyze (layer set = set of layer(s) to jointly analyze)
# 2. Define the organism and ontology for each layer set
# 3. Define keytypes for each individual layer in each layer set
## [Code]
# For each layer set...
# 1. Access the specified ontology
# 2. Determine keytypes and keyclasses contained in ontology terms
# 3. For each layer in the set... 
#   3.a. Determine mapping database to use based on the layer keytype (and ontology keytype?)
#   3.b. Determine the ontology keytype to map layer keytype to based on the keyclasses of layer and ontology keytypes
#   3.c. Map from layer keytype to appropriate ontology keytype
#     NB: For ORA, this mapping must be done for the full network and the module.
# 4. For each term in ontology, concatenate IDs from different keyclasses into a single vector
# 5. Concatenate IDs from different layers within a layer set together into a single vector
# 6. For each layer set, perform pathway analysis (ORA or GSEA) using fgsea package



#=========================================================#
#=========================================================#

#' @title Internal function for pathway analysis
#'
#' @description na
#'
#' @details na
#' 
#' @returns na
#' 
#' @noRd
#'
pathway_analysis_internal <- function(x, db = c("GO", "KEGG", "MKEGG", "Reactome", "WikiPathways", "DO", "DOLite", "NCG", "DisGeNET", "snpDisGeNET"), 
                                      organism = c("human", "mouse", "rat", "fly", "zebrafish", "yeast", "celegans"), ont = c("BP", "MF", "CC", "ALL"), 
                                      keyType, iter = NULL, pvalueCutoff = 0.05, qvalueCutoff = 0.2, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, layers = NULL, ...) {
  if(0) { # TEST
    # devtools::load_all()
    # library(clusterProfiler)
    uniq_names <- unique(c(simple_network[,1], simple_network[,2]))
    p_values <- runif(length(uniq_names))
    names(p_values) <- uniq_names
    x <- AMEND(network_layers = simple_network,
               n = 50,
               data = p_values,
               FUN = "p_value",
               directed = FALSE,
               normalize = "degree",
               in_parallel = TRUE,
               verbose = TRUE)
    db <- "GO"
    OrgDb <- org.Hs.eg.db
    ont <- "BP"
    keyType <- "ENSEMBLPROT"
    iter <- NULL
    pvalueCutoff = 0.05
    qvalueCutoff = 0.2
    pAdjustMethod = "BH"
    minGSSize = 10
    maxGSSize = 500
    layers = NULL
  }
  
  #=== Script Settings ===#
  KEGG_KEYTYPES <- c("kegg", "ncbi-geneid", "ncbi-proteinid", "uniprot")
  #=======================#
  
  if(!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("The clusterProfiler package is required for this function. Please install it using BiocManager::install('clusterProfiler').")
  }
  
  
  db <- match.arg(db)
  ont <- match.arg(ont)
  organism <- match.arg(organism)
  
  
  if(db == "GO") {
    OrgDb <- load_annotation_db(organism)
  } else if(db %in% c("KEGG", "MKEGG")) {
    if(!keyType %in% KEGG_KEYTYPES) stop("For db=KEGG or db=MKEGG, keyType must be one of: ", paste(KEGG_KEYTYPES, collapse = ", "))
    organism <- kegg_organism_mapper(organism)
    if(db == "KEGG") {
      if(method == "ORA") {
        enrich_FUN <- get("enrichKEGG", asNamespace("clusterProfiler"))
      } else enrich_FUN <- get("gseKEGG", asNamespace("clusterProfiler"))
    } else {
      if(method == "ORA") {
        enrich_FUN <- get("enrichMKEGG", asNamespace("clusterProfiler"))
      } else enrich_FUN <- get("gseMKEGG", asNamespace("clusterProfiler"))
    }
  } else if(db == "WikiPathways") {
    organism <- wikipath_organism_mapper(organism)
  }
  
  
  if(inherits(x, "AMENDresult")) {
    method <- "ORA"
    
    if("aggregated_network" %in% names(x)) {
      net_name <- "aggregated_network"
      layer_name <- "agg_layer"
    } else {
      net_name <- "network"
      layer_name <- "layer"
    }
    if(is.character(layers)) {
      universe <- unique(V(x[[net_name]])$original_name[igraph::vertex_attr(x[[net_name]], layer_name) %in% layers])
    } else universe <- unique(V(x[[net_name]])$original_name)
    
    
    if(is.null(iter)) {
      
      if(is.character(layers)) {
        nodeList <- list(unique(V(x$module)$original_name[igraph::vertex_attr(x$module, layer_name) %in% layers]))
      } else nodeList <- list(unique(V(x$module)$original_name))
      
    } else if(is.numeric(iter)) {
      if(!iter %in% seq_along(x$subnetworks)) stop("'iter' must be within the iteration range of 'x' when it is a AMENDresult object.")
      
      g_tmp <- igraph::induced_subgraph(x[[net_name]], which(V(x[[net_name]])$name %in% x$subnetworks[[iter]]))
      
      if(is.character(layers)) {
        nodeList <- list(unique(V(g_tmp)$original_name[igraph::vertex_attr(g_tmp, layer_name) %in% layers]))
      } else nodeList <- list(unique(V(g_tmp)$original_name))
      
    } else if(iter == "all") {
      nodeList <- vector("list", length(x$subnetworks))
      for(i in seq_along(nodeList)) {
        g_tmp <- igraph::induced_subgraph(x[[net_name]], which(V(x[[net_name]])$name %in% x$subnetworks[[i]]))
        if(is.character(layers)) {
          nodeList[[i]] <- unique(V(g_tmp)$original_name[igraph::vertex_attr(g_tmp, layer_name) %in% layers])
        } else nodeList[[i]] <- unique(V(g_tmp)$original_name)
      }
      
    } else stop("'iter' must be numeric, character string 'all', or NULL.")
  } else if(inherits(x, "RWRresult")) {
    method <- "GSEA"
    
    if(is.numeric(x$scores)) {
      nodeList <- list(x$scores)
    } else if(is.list(x$scores)) {
      nodeList <- x$scores
    }
    
    nodeList <- lapply(nodeList, function(x) {
      l <- extract_string(names(x), "\\|", 2)
      names(x) <- extract_string(names(x), "\\|", 1)
      if(is.character(layers)) {
        x[l %in% layers]
      } else x
    })
    
  } else if(is.matrix(x) && is.numeric(x)) {
    method <- "GSEA"
    nodeList <- vector("list", ncol(x))
    for(i in seq_along(nodeList)) {
      nodeList[[i]] <- x[,i]
    }
    nodeList <- lapply(nodeList, function(x) {
      l <- extract_string(names(x), "\\|", 2)
      names(x) <- extract_string(names(x), "\\|", 1)
      if(is.character(layers)) {
        x[l %in% layers]
      } else x
    })
  } else stop("Unrecognized input for 'x'. Must be an 'AMENDresult' object, a 'RWRresult' object, or a numeric matrix.")
  
  
  #========================#
  # Convert IDs in nodeList 
  #========================#
  # WikiPathways, Reactome, DO(Lite), NCG, DisGeNET, snpDisGeNET: Entrez IDs only!
  # (M)KEGG: "kegg", "ncbi-geneid" (i.e., Entrez ID), "ncbi-proteinid", "uniprot"
  # GO: "ENTREZID", "ENSEMBL", "ENSEMBLPROT", "ENSEMBLTRANS", "REFSEQ", "SYMBOL", "UNIGENE", "UNIPROT",...
  # 
  # NB: keyType tells you the type of IDs being used in nodeList.
  # 
  # 1) Determine if the corresponding values of db and keyType are compatible.
  # 2) If not, attempt to map IDs in nodeList from key given by keyType to appropriate key types.
  
  # 1.
  if(!db %in% c("GO", "KEGG", "MKEGG")) {
    if(keyType != "ENTREZID") {
      bitr_fun <- getFromNamespace("bitr", "clusterProfiler")
      # 2.
      OrgDb <- load_annotation_db(organism)
      # nodeList
      for(i in seq_along(nodeList)) {
        map2entrez <- bitr_fun(geneID = nodeList[[i]], fromType = keyType, toType = "ENTREZID", OrgDb = OrgDb, drop = TRUE)
        if(nrow(map2entrez) == 0) stop("No IDs mapped from '", keyType, "' to 'ENTREZID'.")
        map_ids <- match(nodeList[[i]], map2entrez[, keyType]) # Take the first match 
        map_ids <- map_ids[!is.na(map_ids)]
        nodeList[[i]] <- unique(map2entrez[map_ids, "ENTREZID"])
      }
      # universe
      map2entrez <- bitr_fun(geneID = universe, fromType = keyType, toType = "ENTREZID", OrgDb = OrgDb, drop = TRUE)
      if(nrow(map2entrez) == 0) stop("No IDs mapped from '", keyType, "' to 'ENTREZID'.")
      map_ids <- match(universe, map2entrez[, keyType]) # Take the first match 
      map_ids <- map_ids[!is.na(map_ids)]
      universe <- unique(map2entrez[map_ids, "ENTREZID"])
    }
  }
  
  #=================#
  # Pathway Analysis
  #=================#
  if(method == "ORA") {
    # geneList, pvalueCutoff, pAdjustMethod, qvalueCutoff, minGSSize, maxGSSize
    # GO: ont, keyType
    # KEGG: organism, keyType
    # WikiPathways: organism
    # Reactome: organism
    # DOSE: ontology
    enrich_res <- vector("list", length(nodeList))
    for(i in seq_along(enrich_res)) {
      if(db == "GO") {
        enrich_FUN <- getFromNamespace("enrichGO", "clusterProfiler")
        enrich_res[[i]] <- enrich_FUN(gene = nodeList[[i]],
                                                     OrgDb = OrgDb,
                                                     keyType = keyType,
                                                     ont = ont,
                                                     pvalueCutoff = pvalueCutoff,
                                                     pAdjustMethod = pAdjustMethod,
                                                     universe = universe,
                                                     qvalueCutoff = qvalueCutoff,
                                                     minGSSize = minGSSize,
                                                     maxGSSize = maxGSSize,
                                                     ...)
      } else if(db %in% c("KEGG", "MKEGG")) {
        enrich_res[[i]] <- enrich_FUN(gene = nodeList[[i]],
                                      organism = organism,
                                      keyType = keyType,
                                      pvalueCutoff = pvalueCutoff,
                                      pAdjustMethod = pAdjustMethod,
                                      universe = universe,
                                      qvalueCutoff = qvalueCutoff,
                                      minGSSize = minGSSize,
                                      maxGSSize = maxGSSize,
                                      ...)
      } else if(db == "Reactome") {
        enrich_FUN <- getFromNamespace("enrichPathway", "ReactomePA")
        enrich_res[[i]] <- enrich_FUN(gene = nodeList[[i]],
                                                     organism = organism,
                                                     pvalueCutoff = pvalueCutoff,
                                                     pAdjustMethod = pAdjustMethod,
                                                     qvalueCutoff = qvalueCutoff,
                                                     universe = universe,
                                                     minGSSize = minGSSize,
                                                     maxGSSize = maxGSSize,
                                                     ...)
      }else if(db == "WikiPathways") {
        enrich_FUN <- getFromNamespace("enrichWP", "clusterProfiler")
        enrich_res[[i]] <- enrich_FUN(gene = nodeList[[i]],
                                                     organism = organism,
                                                     pvalueCutoff = pvalueCutoff,
                                                     pAdjustMethod = pAdjustMethod,
                                                     universe = universe,
                                                     minGSSize = minGSSize,
                                                     maxGSSize = maxGSSize,
                                                     qvalueCutoff = qvalueCutoff,
                                                     ...)
      }else if(db %in% c("DO", "DOLite", "NCG", "DisGeNET", "snpDisGeNET")) {
        enrich_FUN <- getFromNamespace("enrichDisease", "DOSE")
        enrich_res[[i]] <- enrich_FUN(gene = nodeList[[i]],
                                               pvalueCutoff = pvalueCutoff,
                                               pAdjustMethod = pAdjustMethod,
                                               universe = universe,
                                               minGSSize = minGSSize,
                                               maxGSSize = maxGSSize,
                                               qvalueCutoff = qvalueCutoff,
                                               ontology = db,
                                               ...)
      }
    }
    
    
  } else if(method == "GSEA") {
    enrich_res <- vector("list", length(nodeList))
    for(i in seq_along(enrich_res)) {
      if(db == "GO") {
        enrich_FUN <- getFromNamespace("gseGO", "clusterProfiler")
        enrich_res[[i]] <- enrich_FUN(geneList = nodeList[[i]],
                                                  ont = ont,
                                                  OrgDb = OrgDb,
                                                  keyType = keyType,
                                                  minGSSize = minGSSize,
                                                  maxGSSize = maxGSSize,
                                                  pvalueCutoff = pvalueCutoff,
                                                  pAdjustMethod = pAdjustMethod,
                                                  ...)
      } else if(db %in% c("KEGG", "MKEGG")) {
        enrich_res[[i]] <- enrich_FUN(geneList = nodeList[[i]],
                                      organism = organism,
                                      keyType = keyType,
                                      minGSSize = minGSSize,
                                      maxGSSize = maxGSSize,
                                      pvalueCutoff = pvalueCutoff,
                                      pAdjustMethod = pAdjustMethod,
                                      ...)
      } else if(db == "Reactome") {
        enrich_FUN <- getFromNamespace("gsePathway", "ReactomePA")
        enrich_res[[i]] <- enrich_FUN(geneList = nodeList[[i]],
                                                  organism = organism,
                                                  minGSSize = minGSSize,
                                                  maxGSSize = maxGSSize,
                                                  pvalueCutoff = pvalueCutoff,
                                                  pAdjustMethod = pAdjustMethod,
                                                  ...)
      }else if(db == "WikiPathways") {
        enrich_FUN <- getFromNamespace("gseWP", "clusterProfiler")
        enrich_res[[i]] <- enrich_FUN(geneList = nodeList[[i]], 
                                                  organism = organism, 
                                                  minGSSize = minGSSize,
                                                  maxGSSize = maxGSSize,
                                                  pvalueCutoff = pvalueCutoff,
                                                  pAdjustMethod = pAdjustMethod,
                                                  ...)
      }else if(db %in% c("DO", "DOLite", "NCG", "DisGeNET", "snpDisGeNET")) {
        enrich_FUN <- getFromNamespace("gseDisease", "DOSE")
        enrich_res[[i]] <- enrich_FUN(geneList = nodeList[[i]],
                                            minGSSize = minGSSize,
                                            maxGSSize = maxGSSize,
                                            pvalueCutoff = pvalueCutoff,
                                            pAdjustMethod = pAdjustMethod,
                                            ontology = db,
                                            ...)
      }
    }
  }
  
  return(enrich_res)
}


load_annotation_db <- function(organism) {
  # Create the Bioconductor package name based on organism
  organism_to_pkg <- list("human" = "org.Hs.eg.db",
                         "mouse" = "org.Mm.eg.db",
                         "zebrafish" = "org.Dr.eg.db",
                         "rat" = "org.Rn.eg.db",
                         "celegans" = "org.Ce.eg.db",
                         "yeast" = "org.Sc.sgd.db",
                         "fly" = "org.Dm.eg.db")
  
  # Check if the organism is supported
  if(!organism %in% names(organism_to_pkg)) {
    stop("Annotation database for this organism is not available in this package.")
  }
  
  # Get the corresponding package name
  pkg <- organism_to_pkg[[organism]]
  
  # Dynamically load the package
  if(!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("Please install", pkg, "from Bioconductor."))
  }
  
  # Return the database object (e.g., org.Hs.eg.db)
  return(get(pkg, envir = asNamespace(pkg)))
}


kegg_organism_mapper <- function(organism) {
  organism_to_kegg <- list("human" = "hsa",
                           "mouse" = "mmu",
                           "zebrafish" = "dre",
                           "rat" = "rno",
                           "celegans" = "cel",
                           "yeast" = "sce",
                           "fly" = "dme")
  
  # Check if the organism is supported
  if(!organism %in% names(organism_to_kegg)) {
    stop("Pathway analysis for this organism is not available in this package.")
  }
  
  # Get the corresponding package name
  kegg_name <- organism_to_kegg[[organism]]
  
  return(kegg_name)
}


wikipath_organism_mapper <- function(organism) {
  organism_to_wiki <- list("human" = "Homo sapiens",
                           "mouse" = "Mus musculus",
                           "zebrafish" = "Danio rerio",
                           "rat" = "Rattus norvegicus",
                           "celegans" = "Caenorhabditis elegans",
                           "yeast" = "Saccharomyces cerevisiae",
                           "fly" = "Drosophila melanogaster")
  
  # Check if the organism is supported
  if(!organism %in% names(organism_to_wiki)) {
    stop("Pathway analysis for this organism is not available in this package.")
  }
  
  # Get the corresponding package name
  wiki_name <- organism_to_wiki[[organism]]
  
  return(wiki_name)
}


#' @title Perform Pathway Analysis
#'
#' @description
#' Perform pathway analysis on results from `AMEND()`/`module_identification()` (which will use Over-representation Analysis) and `RWR_pipeline()` (which will use Gene Set Enrichment Analysis).
#'
#' @details
#' Pathway analysis can be performed on specific layers or sets of layers. Additionally, pathway analysis can be performed for sets of parameter settings for the same layer or layer set. This is determinged by the _db_, _organism_, _ont_, _keyType_, and _groups_ arguments. For the supplied db/organism/ont/keyType parameters, it's assumed that if one is a single unnamed character string, it applies to all elements of the vectors. Pathway analysis will then be performed for all parameter sets associated with each set of layers.
#' 
#' The given key type must be compatible with the given database. The GO database can accept a large number of key types, depending on the specific organism being analyzed. See `?AnnotationDbi::ACCNUM` for a complete list. KEGG and MKEGG can accept "kegg", "ncbi-geneid", "ncbi-proteinid", or "uniprot". For the remaining databases, the IDs must be Entrez IDs ("ENTREZID"). If the given database is neither GO, KEGG, or MKEGG, and the given key type is not "ENTREZID", an effort will be made to map the IDs from the given key type to Entrez IDs using `bitr()` from _clusterProfiler_, in which case the given key type must be among those from `?AnnotationDbi::ACCNUM`.   
#' 
#' This function relies heavily on the _clusterProfiler_ package.
#' 
#' @param x a 'AMENDresult' or 'RWRresult' object.
#' @param db named character vector (names = layers), or unnamed character string from the following: "GO", "KEGG", "MKEGG", "Reactome", "WikiPathways", "DO", "DOLite", "NCG", "DisGeNET", "snpDisGeNET"
#' @param organism named character vector (names = layers), or unnamed character string from the following: "human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly". 
#' @param ont named character vector (names = layers), or unnamed character string from the following GO ontologies: "BP" (default), "MF", "CC", or "ALL" for all three.
#' @param keyType named character vector (names = layers), or unnamed character string. The type of identifier for nodes in the given layers. See Details.
#' @param groups A list of character vectors of network layers to be analyzed jointly, or NULL (default), in which case all layers are analyzed separately. 
#' @param iter iteration number of the AMEND algorithm corresponding to the subnetwork to analyze. NULL (default) analyzes the highest scoring subnetwork. Set to "all" to analyze the subnetwork from each iteration of AMEND. 
#' @param pvalueCutoff p-value cutoff
#' @param qvalueCutoff q-value cutoff
#' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param minGSSize minimal size of each pathway
#' @param maxGSSize maximal size of each pathway
#' @param `...` other parameters passed to either `enrich_` (for ORA) or `gse_` (for GSEA) functions in the _clusterProfiler_, _ReactomePA_, or _DOSE_ packages.
#' 
#' @returns A list (for layer sets associated with an organism & keyType) of lists (for db/ont parameter sets) of either 1) lists of enrichment results (for numeric _iter_ argument), or 2) enrichment results.
#' 
# @export
#' @noRd
#'
pathway_analysis <- function(x, db, organism, ont = "BP", keyType = "ENTREZID", groups = NULL, iter = NULL, pvalueCutoff = 0.05, qvalueCutoff = 0.2, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, ...) {
  # db / organism / ont / keyType: named character vector (names = layers), or unnamed character string
  # If some are scalar and some vector, assume the scalar applies to all elements of the vectors
  # For each layer present, gather all relevant arguments in a call to expand.grid().
  # Make calls to pathway_analysis_internal() for each set of arguments for each layer.
  # Each call to pathway_analysis_internal() returns a list.
  # The result of this function will be a list (for different layers / sets of layers, each associated with an organism & keyType) of lists (for different db/ont arguments) of either 1) lists (for different iterations) of enrichment results or 2) of enrichment result 
  
  if(0) { # TEST
    # devtools::load_all()
    # library(clusterProfiler)
    uniq_names <- unique(c(simple_network[,1], simple_network[,2]))
    p_values <- runif(length(uniq_names))
    names(p_values) <- uniq_names
    x = AMEND(network_layers = simple_network,
               n = 50,
               data = p_values,
               FUN = "p_value",
               directed = FALSE,
               normalize = "degree",
               in_parallel = TRUE,
               verbose = TRUE)
    
    # Ex 1
    db = c(prot = "GO", protdat = "KEGG")
    organism = "human"
    ont = c(protdat = "BP", phosphdat = "MF")
    keyType = "ENSEMBLPROT"
    groups = list("prot", "meta")
    net_hier = create_network_hierarchy(multilayer_hierarchy)
    # Ex 2
    db = c(prot = "GO", metadat = "KEGG")
    organism = "human"
    ont = c(prot = "BP")
    keyType = c(prot = "ENSEMBLPROT", metadat = "SYMBOL")
    groups = list("prot", "meta")
    uniq_layers = c("prot", "metadat", "protdat2", "metadat2")
    # Ex 3
    db = c("GO", "Reactome")
    organism = "human"
    ont = c("BP", "MF")
    # keyType = c("SYMBOL", "ENSEMBLPROT")
    keyType = "ENSEMBLPROT"
    groups = NULL
    # Ex 4
    db = "KEGG"
    organism = "human"
    ont = "BP"
    keyType = "ENSEMBLPROT"
    groups = NULL
    
    iter = NULL
    pvalueCutoff = 0.05
    qvalueCutoff = 0.2
    pAdjustMethod = "BH"
    minGSSize = 10
    maxGSSize = 500
  }
  
  
  if(inherits(x, "AMENDresult") || inherits(x, "RWRresult")) {
    if("aggregated_network" %in% names(x)) {
      uniq_layers <- unique(V(x$aggregated_network)$agg_layer)
    } else {
      uniq_layers <- unique(V(x$network)$layer)
    }
  } else stop("Unrecognized input for 'x'. Must be an 'AMENDresult' object or a 'RWRresult' object.")
  
  
  #=== db / organism / ont / keyType processing ===#
  # Ensure that each layer present has exactly 1 of organism and keyType.
  fan <- c("db", "organism", "ont", "keyType")
  FUN_args <- list(db, organism, ont, keyType); names(FUN_args) <- fan
  FUN_args <- lapply(names(FUN_args), function(x) {
    if(is.null(names(FUN_args[[x]]))) {
      uniq_x <- unique(FUN_args[[x]])
      FUN_args[[x]] <- rep(uniq_x, each = length(uniq_layers))
      names(FUN_args[[x]]) <- rep(uniq_layers, times = length(uniq_x))
    } else {
      FUN_args[[x]] <- FUN_args[[x]][names(FUN_args[[x]]) %in% uniq_layers]
    }
    if(x %in% c("organism", "keyType")) {
      tmp_nm <- unique(names(FUN_args[[x]]))
      new_arg <- character(length(tmp_nm)); names(new_arg) <- tmp_nm
      for(i in seq_along(new_arg)) {
        new_arg[i] <- FUN_args[[x]][tmp_nm[i]]
      }
      FUN_args[[x]] <- new_arg
    }
    FUN_args[[x]]
  })
  names(FUN_args) <- fan
  # FUN_args <- lapply(list(db = db, organism = organism, ont = ont, keyType = keyType), function(x) {
  #   if(is.null(names(x))) {
  #     uniq_x <- unique(x)
  #     x <- rep(uniq_x, each = length(uniq_layers))
  #     names(x) <- rep(uniq_layers, times = length(uniq_x))
  #   } else {
  #     x <- x[names(x) %in% uniq_layers]
  #   }
  #   x
  # })
  
  uniq_cats <- unique(unlist(lapply(FUN_args, names)))
  
  arg_sets <- vector("list", length(uniq_cats)); names(arg_sets) <- uniq_cats
  for(i in seq_along(arg_sets)) {
    tmp_FUN_args <- lapply(names(FUN_args), function(x) {
      FUN_args[[x]][names(FUN_args[[x]]) %in% names(arg_sets)[i]]
    })
    names(tmp_FUN_args) <- names(FUN_args)
    
    if(length(unlist(tmp_FUN_args)) < 4) next
    # Only need a special case for when db includes GO and at least one other, and ont has more than one.
    tmp_df <- expand.grid(db = tmp_FUN_args[["db"]],
                          organism = tmp_FUN_args[["organism"]],
                          keyType = tmp_FUN_args[["keyType"]],
                          ont = tmp_FUN_args[["ont"]], stringsAsFactors = FALSE)
    if("GO" %in% tmp_df$db && length(setdiff(tmp_df$db, "GO")) > 0 && length(unique(tmp_df$ont)) > 1) {
      tmp_df1 <- expand.grid(db = "GO",
                             organism = tmp_FUN_args[["organism"]],
                             keyType = tmp_FUN_args[["keyType"]],
                             ont = tmp_FUN_args[["ont"]], stringsAsFactors = FALSE)
      tmp_df2 <- expand.grid(db = tmp_FUN_args[["db"]][tmp_FUN_args[["db"]] != "GO"],
                             organism = tmp_FUN_args[["organism"]],
                             keyType = tmp_FUN_args[["keyType"]],
                             ont = "BP", stringsAsFactors = FALSE)
      tmp_df <- rbind(tmp_df1, tmp_df2)
    }
    arg_sets[[i]] <- tmp_df
  }
  arg_sets <- arg_sets[!sapply(arg_sets, is.null)]
  if(length(arg_sets) == 0) stop("Incomplete information for Pathway Analysis. Please fully specify db, organism, keyType, and ont (for db='GO').")
  
  #=== group processing ===#
  if(is.null(groups)) {
    layer_sets <- names(arg_sets)
  } else if(is.list(groups)) {
    layer_sets <- sapply(groups, function(x) {
      x <- x[x %in% names(arg_sets)]
      ifelse(length(x) == 0, NA, paste(x, collapse = "|"))
    })
    layer_sets <- layer_sets[!is.na(layer_sets)]
    if(length(layer_sets) == 0) {
      layer_sets <- names(arg_sets)
    } else {
      for(i in seq_along(layer_sets)) {
        ids <- which(names(arg_sets) %in% extract_string(layer_sets[i], "\\|", 0))
        arg_sets[[length(arg_sets) + 1]] <- unique(do.call(rbind, arg_sets[ids]))
        names(arg_sets)[length(arg_sets)] <- layer_sets[i]
        arg_sets <- arg_sets[-ids]
        if(length(unique(arg_sets[[length(arg_sets)]]$organism)) > 1 || length(unique(arg_sets[[length(arg_sets)]]$keyType)) > 1) {
          stop("Attmepting to jointly analyze layers corresponding to different organisms or key types.")
        }
      }
    }
  } else stop("Unrecognized input for 'groups'. Must be a list or NULL.")
  
  
  enrich_res <- vector("list", length(arg_sets)); names(enrich_res) <- names(arg_sets)
  for(i in seq_along(enrich_res)) {
    enrich_res[[i]] <- vector("list", nrow(arg_sets[[i]]))
    for(j in seq_along(enrich_res[[i]])) {
      enrich_res[[i]][[j]] <- pathway_analysis_internal(x = x,
                                                        db = arg_sets[[i]][j,"db"],
                                                        organism = arg_sets[[i]][j,"organism"],
                                                        ont = arg_sets[[i]][j,"ont"],
                                                        keyType = arg_sets[[i]][j,"keyType"],
                                                        iter = iter,
                                                        pvalueCutoff = pvalueCutoff,
                                                        qvalueCutoff = qvalueCutoff,
                                                        pAdjustMethod = pAdjustMethod,
                                                        minGSSize = minGSSize,
                                                        maxGSSize = maxGSSize,
                                                        layers = extract_string(names(arg_sets)[i], "\\|", 0))
      attr(enrich_res[[i]][[j]], "settings") <- arg_sets[[i]][j,]
    }
  }
  
  return(enrich_res)
}

# TEST
if(0) {
  set.seed(593)
  # devtools::load_all()
  # library(clusterProfiler)
  uniq_names <- unique(c(simple_network[,1], simple_network[,2]))
  p_values <- runif(length(uniq_names))
  names(p_values) <- uniq_names
  mod <- AMEND(network_layers = simple_network,
               n = 50,
               data = p_values,
               FUN = "p_value",
               directed = FALSE,
               normalize = "degree",
               in_parallel = TRUE,
               verbose = TRUE)
  
  pa <- pathway_analysis(x = mod,
                         db = "GO",
                         organism = "human",
                         ont = "BP",
                         keyType = "ENSEMBLPROT",
                         iter = 5,
                         pvalueCutoff = 0.5,
                         qvalueCutoff = 0.9,
                         pAdjustMethod = "BH",
                         minGSSize = 10,
                         maxGSSize = 500,
                         layers = NULL)
  
  pa <- pathway_analysis(x = mod,
                         db = "Reactome",
                         organism = "human",
                         keyType = "ENSEMBLPROT",
                         iter = NULL,
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2,
                         pAdjustMethod = "BH",
                         minGSSize = 10,
                         maxGSSize = 500,
                         layers = NULL)
}









