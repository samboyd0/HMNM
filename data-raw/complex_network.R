#================#
#================#
# Complex Network
#================#
#================#

rm(list=ls())

# Attach required libraries
library(igraph)
library(dplyr)

#==========================#
# Network Utility Functions
#==========================#
# adapted from Python code provided at stringdb website
new_combined_score <- function(ppi, evidence_channels, prior = 0.041) {
  compute_prior_away <- function(score, prior = 0.041) {
    score <- ifelse(score < prior, prior, score)
    score_no_prior <- (score - prior) / (1 - prior)
    return(score_no_prior)
  }
  multiply_elements <- function(lst) {
    # Base case: if the list has only one element, return it
    if(length(lst) == 1) {
      return(1 - lst[[1]])
    }else {
      # Recursive case: multiply the first element with the result of the function on the rest of the list
      return((1 - lst[[1]]) * multiply_elements(lst[-1]))
    }
  }
  
  new_score <- vector('list', length(evidence_channels))
  names(new_score) <- evidence_channels
  for(i in seq_along(new_score)) {
    # compute prior away
    new_score[[i]] <- compute_prior_away(score = ppi[,names(new_score)[i]] / 1000, prior = prior)
  }
  
  # then, combine the direct and transferred scores for each category
  channels <- c('experiments', 'database', 'textmining', 'experimental', 'prediction')
  for(i in seq_along(channels)) {
    if(all(c(channels[i], paste0(channels[i], '_transferred')) %in% evidence_channels)) {
      new_score[[channels[i]]] <- 1 - (1 - new_score[[channels[i]]]) * (1 - new_score[[paste0(channels[i], '_transferred')]])
    }else if(all(c(paste0(channels[i], '_direct'), paste0(channels[i], '_transferred')) %in% evidence_channels)) {
      new_score[[paste0(channels[i], '_direct')]] <- 1 - (1 - new_score[[paste0(channels[i], '_direct')]]) * (1 - new_score[[paste0(channels[i], '_transferred')]])
    }
  }
  new_score <- new_score[!grepl('transferred', names(new_score))]
  
  ## next, do the 1 - multiplication:
  combined_score_one_minus <- multiply_elements(new_score)
  
  ## and lastly, do the 1 - conversion again, and put back the prior *exactly once*
  
  combined_score <- (1.0 - combined_score_one_minus)            ## 1- conversion
  combined_score <- combined_score * (1.0 - prior)              ## scale down
  combined_score <- combined_score + prior                      ## and add prior.
  
  ## round
  
  # combined_score <- round(combined_score * 1000, 0)
  combined_score <- round(combined_score, 3)
  
  return(combined_score)
}

# Remove duplicates and Aggregate edge scores
remove_duplicate_edges <- function(el) {
  tmp <- character(nrow(el))
  for(i in seq_along(tmp)) {
    tmp[i] <- paste(sort(c(el[i,1], el[i,2])), collapse = "|")
  }
  scores <- as.numeric(el[,3])
  dat <- stats::aggregate(scores, by = list(tmp), FUN = min)
  el <- matrix(c(extract_string(dat[,1], "\\|", 1), extract_string(dat[,1], "\\|", 2), dat[,2]), ncol = 3)
  return(el)
}

extract_string <- function(x, k, pos) {
  if(pos == 0) {
    res <- unlist(lapply(strsplit(x, k), function(y) y[1:length(y)]))
  } else {
    res <- unlist(lapply(strsplit(x, k), function(y) y[min(length(y), pos)]))
  }
  return(res)
} 

set.seed(20)

# Example 1 ----
# Including 'root' node to help with processing in the cases with one level
# May generalize better to multi-scale networks... this 'root' node could be the name(s) of leaf node for the scale above to which this hierarchy applies...
# This root node must be the only node with in-degree = 0... 
network_hierarchy <- matrix(c('root', 'blood',
                              'root', 'liver',
                              'liver', 'meta',
                              'liver', 'prot',
                              'meta', 'metadat',
                              'prot', 'protdat',
                              'prot', 'phosphdat',
                              'blood', 'meta2',
                              'blood', 'prot2',
                              'meta2', 'metadat2',
                              'prot2', 'protdat2'), ncol = 2, byrow = TRUE)
#=================#
# Practice Network 
#=================#
edge_threshold <- 0.7

# prior of 0.041 from STRING v12
prior <- 0.041

# Load in edgelists
ppi <- data.table::fread(file = '/Users/samboyd/Documents/HMNM/R package/Data/9606.protein.physical.links.full.v12.0.txt.gz', header = TRUE, data.table = FALSE)
ppi <- ppi %>%
  dplyr::mutate(combined_score = new_combined_score(., evidence_channels = c('experiments', 'experiments_transferred', 
                                                                             'database', 'database_transferred'), prior = prior)) %>%
  dplyr::filter(combined_score >= edge_threshold) %>%
  dplyr::mutate(protein1 = extract_string(protein1, "\\.", 2), # remove taxonomy identifier
                protein2 = extract_string(protein2, "\\.", 2)) %>%
  dplyr::select(protein1, protein2, combined_score) %>%
  dplyr::distinct() %>%
  as.matrix() %>%
  remove_duplicate_edges()

mmi <- data.table::fread(file = '/Users/samboyd/Documents/HMNM/R package/Data/chemical_chemical.links.detailed.v5.0.tsv.gz', header = TRUE, data.table = FALSE)
mmi <- mmi %>%
  dplyr::mutate(combined_score = new_combined_score(., evidence_channels = c('similarity', 'experimental', 'database'), prior = prior)) %>%
  dplyr::filter(combined_score >= edge_threshold) %>%
  dplyr::mutate(chemical1 = chemical1,
                chemical2 = chemical2) %>% # remove taxonomy identifier
  dplyr::select(chemical1, chemical2, combined_score) %>%
  dplyr::distinct() %>%
  as.matrix() %>%
  remove_duplicate_edges()

bp <- data.table::fread(file = '/Users/samboyd/Documents/HMNM/R package/Data/9606.protein_chemical.links.transfer.v5.0.tsv.gz', header = TRUE, data.table = FALSE)
bp <- bp %>%
  dplyr::mutate(combined_score = new_combined_score(., evidence_channels = c('experimental_direct', 'experimental_transferred', 
                                                                             'prediction_direct', 'prediction_transferred', 
                                                                             'database_direct', 'database_transferred'), prior = prior)) %>%
  dplyr::filter(combined_score >= edge_threshold) %>%
  dplyr::mutate(chemical = chemical,
                protein = extract_string(protein, "\\.", 2)) %>% # remove taxonomy identifier
  dplyr::select(chemical, protein, combined_score) %>%
  dplyr::distinct() %>%
  as.matrix() %>%
  remove_duplicate_edges()

# Create igraph objects
g.mmi <- igraph::graph_from_edgelist(el = mmi[,1:2], directed = FALSE)
E(g.mmi)$weight <- as.numeric(mmi[,3])

g.ppi <- igraph::graph_from_edgelist(el = ppi[,1:2], directed = FALSE)
E(g.ppi)$weight <- as.numeric(ppi[,3])

g.bp1 <- igraph::graph_from_edgelist(el = bp[,1:2], directed = FALSE)
E(g.bp1)$weight <- as.numeric(bp[,3])
V(g.bp1)$layer <- ifelse(substr(V(g.bp1)$name,1,3) == 'CID', 'meta', 'prot')

g.bp2 <- g.bp1
V(g.bp2)$layer <- ifelse(substr(V(g.bp2)$name,1,3) == 'CID', 'meta2', 'prot2')

# Create two separate networks for each molecular type
clust <- cluster_louvain(g.ppi, resolution = 0.2)
c.ids <- as.numeric(names(sort(table(clust$membership), decreasing=T)))
g.ppi1 <- induced_subgraph(g.ppi, which(clust$membership == c.ids[1]))
# g.ppi2 <- induced_subgraph(g.ppi, which(clust$membership == c.ids[2]))
clust <- cluster_louvain(g.ppi1, resolution = 0.2)
c.ids <- as.numeric(names(sort(table(clust$membership), decreasing=T)))
g.ppi2 <- induced_subgraph(g.ppi1, which(clust$membership == c.ids[1]))

clust <- cluster_louvain(g.mmi, resolution = 1.5)
c.ids <- as.numeric(names(sort(table(clust$membership), decreasing=T)))
g.mmi1 <- induced_subgraph(g.mmi, which(clust$membership == c.ids[1]))
# g.mmi2 <- induced_subgraph(g.mmi, which(clust$membership == c.ids[2]))
clust <- cluster_louvain(g.mmi1, resolution = 1)
c.ids <- as.numeric(names(sort(table(clust$membership), decreasing=T)))
g.mmi2 <- induced_subgraph(g.mmi1, which(clust$membership == c.ids[1]))

# Create a integrated graph of proteins and metabolites
ppi.tmp <- cbind(igraph::as_edgelist(g.ppi2), E(g.ppi2)$weight)
mmi.tmp <- cbind(igraph::as_edgelist(g.mmi2), E(g.mmi2)$weight)
ppi.names <- V(g.ppi2)$name
mmi.names <- V(g.mmi2)$name
bp.tmp <- bp[apply(bp[,1:2], 1, function(x) all(x %in% c(ppi.names, mmi.names))),]
el <- rbind(ppi.tmp, mmi.tmp, bp.tmp)
g.int <- igraph::graph_from_edgelist(el[,1:2], directed = FALSE)
E(g.int)$weight <- el[,3]
V(g.int)$layer <- ifelse(V(g.int)$name %in% ppi.names, 'protdat2', 'metadat2')

# All of these are equivalent
if(1) {
  multilayer_network <- list(metadat = g.mmi1, 
                         protdat = g.ppi1, 
                         phosphdat = g.ppi1, 
                         metadat2 = g.mmi2, 
                         protdat2 = g.ppi2)
  interlayer_links <- list("protdat|phosphdat" = "common",
                             "prot|meta" = g.bp1,
                             "prot2|meta2" = g.bp2,
                             "metadat|metadat2" = "common",
                             "protdat|protdat2" = "common")
} 
if(0) { 
  # Or... for an example of merged networks in the layer input list
  # There will be duplicate edges for metadat2 <--> protdat2 mapping. Duplicate edges will be removed during processing
  multilayer_network <- list(metadat = g.mmi1, 
                         protdat = g.ppi1, 
                         phosphdat = g.ppi1, 
                         "metadat2|protdat2" = g.int)
  interlayer_links <- list("protdat|phosphdat" = "common",
                             "prot|meta" = g.bp1,
                             "prot2|meta2" = g.bp2,
                             "metadat|metadat2" = "common",
                             "protdat|protdat2" = "common")
} 
if(0) {
  # Same as the above
  multilayer_network <- list(metadat = g.mmi1, 
                         protdat = g.ppi1, 
                         phosphdat = g.ppi1, 
                         "blood" = g.int)
  interlayer_links <- list("protdat|phosphdat" = "common",
                             "prot|meta" = g.bp1,
                             "prot2|meta2" = g.bp2,
                             "metadat|metadat2" = "common",
                             "protdat|protdat2" = "common")
  # Equivalent
  if(0){
    interlayer_links <- list("protdat|phosphdat" = "common",
                               "prot|meta" = g.bp1,
                               "blood" = g.bp2,
                               "metadat|metadat2" = "common",
                               "protdat|protdat2" = "common")
  }
}


#=== Other Function Arguments ===#
if(0){
  data_type <- c("normal", "exponential", "uniform")[1]
  if(data_type == "exponential") {
    data <- list(protdat = rexp(vcount(g.ppi1)),
                 phosphdat = rexp(vcount(g.ppi1)),
                 metadat = rexp(vcount(g.mmi1)),
                 protdat2 = rexp(vcount(g.ppi2)),
                 metadat2 = rexp(vcount(g.mmi2)))
    names(data$protdat) <- names(data$phosphdat) <- V(g.ppi1)$name
    names(data$protdat2) <- V(g.ppi2)$name
    names(data$metadat) <- V(g.mmi1)$name
    names(data$metadat2) <- V(g.mmi2)$name
    
    # exponential alternative for 'data'
    if(0) { 
      data <- list(prot = rexp(vcount(g.ppi1)),
                   metadat = rexp(vcount(g.mmi1)),
                   protdat2 = rexp(vcount(g.ppi2)),
                   metadat2 = rexp(vcount(g.mmi2)))
      names(data$prot) <- V(g.ppi1)$name
      names(data$protdat2) <- V(g.ppi2)$name
      names(data$metadat) <- V(g.mmi1)$name
      names(data$metadat2) <- V(g.mmi2)$name
    }
    
    FUN <- list(prot = function(x, d) d * sqrt(x),
                protdat2 = function(x, d) d * sqrt(x),
                metadat = NULL,
                metadat2 = NULL)
    FUN_params <- list(protdat = list(d = 2), 
                       phosphdat = list(d = 1),
                       protdat2 = list(d = 0.5))
    
  } else if(data_type == "normal") {
    data <- list(protdat = rnorm(vcount(g.ppi1)),
                phosphdat = rnorm(vcount(g.ppi1)),
                metadat = rnorm(vcount(g.mmi1)),
                protdat2 = rnorm(vcount(g.ppi2)),
                metadat2 = rnorm(vcount(g.mmi2)))
    names(data$protdat) <- names(data$phosphdat) <- V(g.ppi1)$name
    names(data$protdat2) <- V(g.ppi2)$name
    names(data$metadat) <- V(g.mmi1)$name
    names(data$metadat2) <- V(g.mmi2)$name
    
    FUN <- "shift_scale"
    FUN_params <- list(protdat = list(DOI = 1, w = 0.2),
                       phosphdat = list(DOI = -1, w = 0.5),
                       metadat = list(DOI = 1, w = 0.2),
                       protdat2 = list(DOI = 1, w = 0.4),
                       metadat2 = list(DOI = 1, w = 0.5))
    
    
  } else if(data_type == "uniform") {
    data <- list(protdat = runif(vcount(g.ppi1)),
                 phosphdat = runif(vcount(g.ppi1)),
                 metadat = runif(vcount(g.mmi1)),
                 protdat2 = runif(vcount(g.ppi2)),
                 metadat2 = runif(vcount(g.mmi2)))
    names(data$protdat) <- names(data$phosphdat) <- V(g.ppi1)$name
    names(data$protdat2) <- V(g.ppi2)$name
    names(data$metadat) <- V(g.mmi1)$name
    names(data$metadat2) <- V(g.mmi2)$name
    
    multilayer_data_runif <- data
    
    FUN <- "p_value"
  }
  
  #=== brw_attr ===#
  brw_attr <- list(protdat = rexp(vcount(g.ppi1)),
                   phosphdat = rexp(vcount(g.ppi1)),
                   metadat = rexp(vcount(g.mmi1)),
                   protdat2 = rexp(vcount(g.ppi2)),
                   metadat2 = rexp(vcount(g.mmi2)))
  names(brw_attr$protdat) <- names(brw_attr$phosphdat) <- V(g.ppi1)$name
  names(brw_attr$protdat2) <- V(g.ppi2)$name
  names(brw_attr$metadat) <- V(g.mmi1)$name
  names(brw_attr$metadat2) <- V(g.mmi2)$name
  
  directed <- FALSE
  lcc <- FALSE
  in_parallel <- FALSE
  n_cores <- NULL
  
  n <- 50
  normalize <- "degree"
  k <- 0.5
  
  # crosstalk_params <- c(protdat = 0.8, phosphdat = 0.4, meta = 0.2)
  crosstalk_params <- c(protdat = 0.2, phosphdat = 0.4, meta = 0.6)
  
  # degree_bias <- list(layers = "prot", gamma = 0.15)
  degree_bias <- list(layers = "meta", gamma = 0.15)
  
  aggregate_layers <- list(c("protdat", "phosphdat"), c("metadat", "metadat2"), method = "mean")
  # aggregate_layers <- NULL
  
  seed_weights <- list(c(protdat = 0.4, phosphdat = 0.6),
                       c(prot = 0.7, meta = 0.3),
                       c(prot2 = 0.7, meta2 = 0.3))
  
  verbose <- FALSE
  eta <- NULL
  
}


# Simulate p-values from uniform distribution
data <- list(protdat = runif(vcount(g.ppi1)),
             phosphdat = runif(vcount(g.ppi1)),
             metadat = runif(vcount(g.mmi1)),
             protdat2 = runif(vcount(g.ppi2)),
             metadat2 = runif(vcount(g.mmi2)))
names(data$protdat) <- names(data$phosphdat) <- V(g.ppi1)$name
names(data$protdat2) <- V(g.ppi2)$name
names(data$metadat) <- V(g.mmi1)$name
names(data$metadat2) <- V(g.mmi2)$name

multilayer_data_runif <- data

multilayer_hierarchy <- network_hierarchy

usethis::use_data(multilayer_data_runif, overwrite = TRUE)
usethis::use_data(multilayer_hierarchy, overwrite = TRUE)
usethis::use_data(multilayer_network, overwrite = TRUE)
usethis::use_data(interlayer_links, overwrite = TRUE)

