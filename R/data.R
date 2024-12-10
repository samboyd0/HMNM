#==================#
# Data Descriptions
#==================#


#' @title Multilayer Network
#'
#' @description
#' `multilayer_network` is a list of 5 igraph objects representing layers of a multilayer network. Protein-protein interactions (PPIs) were obtained from STRING (v12.0), and chemical-chemical & chemcial-protein interactions were obtained from STITCH (v5.0).
#' 
#' The 5 igraph objects were obtained by getting clusters (Louvain clustering method) for STRING and STITCH networks and taking the two largest clusters from each network. One cluster is included twice to represent another layer.
#'
#' @format list of 5 igraph objects with names 'metadat', 'protdat', 'phosphdat', 'metadat2' and 'protdat2'.
#'
#' @source \url{https://string-db.org/}, \url{http://stitch.embl.de/}
"multilayer_network"


#' @title Interlayer links for multilayer network
#'
#' @description
#' `interlayer_links` is a list of 3 character scalars of keyword "common" and 2 igraph objects containing interlayer links. Chemcial-protein interactions (i.e., interlayer links) were obtained from STITCH (v5.0).
#'
#' @format list with names 'protdat|phosphdat', 'prot|meta', 'prot2|meta2', 'metadat|metadat2' and 'protdat|protdat2'.
#'
#' @source \url{http://stitch.embl.de/}
"interlayer_links"


#' @title Hierarchy for multilayer network
#'
#' @description
#' `multilayer_hierarchy` is a 2-column matrix representing a directed edgelist for the hierarchy of the multilayer network given in `multilayer_network`.
#'
#' @format 2-column matrix of a directed edgelist
#'
"multilayer_hierarchy"


#' @title Node-wise data for multilayer network
#'
#' @description
#' `multilayer_data_runif` is a list containing vectors of data randomly generated from a Uniform(0,1) distribution representing p-values from some 'omic experiment. 
#'
#' @format list of numeric vectors with names 'metadat', 'protdat', 'phosphdat', 'metadat2' and 'protdat2'.
#'
"multilayer_data_runif"


#' @title Simple Network
#'
#' @description
#' `simple_network` is a 3-column matrix represnting an edgelist with an additional columm of edge weights. These edges are protein-protein interactions (PPIs) obtained from STRING (v12.0).
#'
#' @format 3-column matrix represnting an edgelist with an additional columm of edge weights.
#'
#' @source \url{https://string-db.org/}
"simple_network"

