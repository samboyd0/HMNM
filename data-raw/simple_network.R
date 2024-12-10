#================#
#================#
# Simple Network
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

# Example 5 ----
# 1 level, 1 category (i.e., monoplex-homogeneous graph, using old terminology) 
# network_hierarchy <- matrix(c('root', 'prot'), ncol = 2, byrow = TRUE)
network_hierarchy <- NULL

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
  as.matrix() 

ppi <- ppi[sample(1:nrow(ppi), 0.4 * nrow(ppi)),]

simple_network <- ppi

# Other function arguments
if(0){
  uniq_names <- unique(c(ppi[,1], ppi[,2]))
  bipartite_networks <- NULL
  directed <- TRUE
  
  data <- runif(length(uniq_names))
  names(data) <- uniq_names
  
  brw_attr <- rexp(length(uniq_names))
  names(brw_attr) <- uniq_names
  
  FUN <- "p_value"
  FUN_params <- NULL
  
  n <- 50
  aggregate_layers <- NULL
  normalize <- "degree"
  k <- 0.5
  degree_bias <- NULL
  crosstalk_params <- NULL
  seed_weights <- NULL
  verbose <- FALSE
  eta <- NULL
  
  in_parallel <- TRUE
  n_cores <- NULL
}

usethis::use_data(simple_network, overwrite = TRUE)
