#' @importFrom igraph vcount V V<- E E<- vertex_attr vertex_attr<-
#' @importFrom stats sd quantile


#' @title Identify modules from a network
#' 
#' @description 
#' Given a network (or several networks, to be connected) with node-wise data, `module_identification()` uses the __AMEND__ algorithm to identify a single module containing highly-connected nodes with large data values. `AMEND()` is an alias.
#'
#' @section Algorithm Overview: 
#' 
#' The AMEND algorithm iteratively performs Random Walk with Restart (RWR) to obtain node scores which are used to heuristically find a maximum-weight connected subgraph. This subgraph is input for the next iteration, and the process continues until the subgraph size is close to the user-specified size __n__ or there is no change in network between iterations. See the [manuscript](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10324253/) for more details.
#'
#' @section Network Hierarchy: 
#' 
#' __AMEND__ is powered by a generalized RWR algorithm that can analyze multilayer networks where nodes are categorized by overlapping factors organized hierarchically. This hierarchy is supplied by the user as a directed edeglist and governs the construction of the transition matrix for RWR, with information being shared more readily between layers that lie closer in the hierarchy. 
#' 
#' The hierarchy can have an arbitrary number of levels, where each level is associated with a factor (e.g., tissue, molecule, data type). Each leaf in the hierarchy (i.e., category with no child categories, bottom level) corresponds to a layer in the multilayer network, and the nodes in a layer are defined by the categories of its ancestors in the hierarchy. See [create_network_hierarchy()] for more details.
#' 
#' @section Network Inputs: 
#' 
#' Users can supply layers as separate network-like objects (igraph, adjacency matrix, or edgelist) in a named list, where names correspond to categories in hierarchy. Similar inputs can be given for bipartite mappings between layers, where list names denote the layer names to which the mapping will be applied, with multiple layer names separated by "|". Networks can be (un)weighted and/or (un)directed. A single, connected network will be constructed from the __network_layers__ and __bipartite_networks__ inputs (see [create_integrated_network()]).
#' 
#' @section Seed Value and Transition Matrix Calculation: 
#' 
#' Seed values for RWR are computed from __data__, __FUN__, and __FUN_params__. Users can also supply node-wise values in __brw_attr__ for a biased random walk, where larger values will increase transition probabilities to nodes during RWR. The __brw_attr__ argument can also be a character vector of node names. Continuous values for biased random walk will be generated for each node as an inverse function of distance to the given nodes.
#' 
#' Transition matrices can be created through classic column-normalization (__normalize__ = "degree") or through a modified column-normalization (__normalize__ = "modified_degree") which penalizes transitions to nodes as a function of node degree (to combat degree bias). See [transition_matrix()] for more details on the transition matrix construction process. 
#' 
#' The __degree_bias__ argument can further mitigate degree bias by applying a degree bias adjustment method to specific layers in the network (see [bistochastic_scaling()]).
#' 
#' @section Aggregating Layers: 
#' 
#' In multilayer networks, there is the option to aggregate sets of layers after [RWR()] but prior to [solve_MWCS()]. For a given set of layers, the first layer listed is designated as the 'primary' layer. All nodes and edges of the primary layer will be kept in the aggregated network. Then, moving sequentially through the remaining layers in the set, edges (along with their adjacent nodes) are added only if at least one of its adjacent nodes is not already in the aggregated layer. RWR scores of common nodes will be averaged. This can be useful if certain layers only differ by their edge type (e.g., physical vs. functional interactions between proteins) and the user wants the module to only contain edges of a certain type.
#' 
#' @section Additional Parameters: 
#' 
#' __crosstalk_params__ and __seed_weights__ are additional sets of parameters that allow fine control over diffusion dynamics during [RWR()].
#' 
#' __crosstalk_params__ are parameters that control the amount of information shared between layers. They represent the probability of the random walker to jump from the current layer to another through a bipartite edge during RWR.
#' 
#' __seed_weights__ represent the relative weight to give to seeds values of a certain layer or sets of layers. These should be supplied for all sets of categories in hierarchy that share a common parent (i.e., siblings).
#' 
#' @note
#' `module_identification()` uses _node_specific_restart_=TRUE for `RWR()`, with _node_connectivity_score_ coming from `node_connectivity_score(..., inverse=TRUE, mode="core")`.
#' 
#' @inheritParams create_integrated_network
#' @inheritParams RWR
#' @inheritParams transition_matrix
#' @param n Final module size which the algorithm will try to approximate. NULL to get maximum-scoring subnetwork.
#' @param aggregate_layers NULL (default) or a list containing sets of layer names to aggregate after RWR but before node filtering. Set to NULL for no aggregation. The first layer in a set is taken as the primary layer, whose edges will be preferentially used during the node filtering step in AMEND.
#' @param verbose Logical. Whether to print information on algorithm progress.
#' @param eta Numeric scalar or NULL (default). Starting filtering rate for AMEND algorithm. If NULL, this is set automatically to approximate the user-specified final module size _n_.
#' 
#' @return a named list with the following elements:
#' * module: igraph of the final module (i.e., subnetwork)
#' * score: final module score
#' * subnetworks: a list of node names contained in intermediate subnetworks
#' * stats: data.frame of network and algorithm statistics
#' * time: run time
#' * input_params: list of input parameters
#' 
#' @aliases AMEND
#' 
#' @seealso [create_integrated_network()], [create_network_hierarchy], [transition_matrix()], [RWR()], [solve_MWCS()]
#' 
#' @examples
#' # Attach igraph library
#' library(igraph)
#' 
#' # Inspect the multilayer network objects included in the package
#' multilayer_network # List object with 5 igraphs
#' interlayer_links 
#' multilayer_hierarchy # matrix representing a directed edgelist
#' 
#' # Inspect data to be mapped to nodes in multilayer network
#' # Represents p-values generated from runif
#' lapply(multilayer_data_runif, head)
#' 
#' crosstalk_params <- c(protdat = 0.2, phosphdat = 0.4, meta = 0.6)
#' 
#' seed_weights <- list(c(protdat = 0.4, phosphdat = 0.6), 
#'                      c(prot = 0.7, meta = 0.3),
#'                      c(prot2 = 0.7, meta2 = 0.3))
#' 
#' # Identify an active module from multilayer network
#' mod <- module_identification(network_layers = multilayer_network,
#'                             bipartite_networks = interlayer_links,
#'                             network_hierarchy = multilayer_hierarchy,
#'                             n = 50,
#'                             data = multilayer_data_runif,
#'                             FUN = "p_value",
#'                             normalize = "degree",
#'                             crosstalk_params = crosstalk_params,
#'                             seed_weights = seed_weights)
#' # Use alias 'AMEND()'
#' if(FALSE) {
#' mod <- AMEND(network_layers = multilayer_network,
#'              bipartite_networks = interlayer_links,
#'              network_hierarchy = multilayer_hierarchy,
#'              n = 50,
#'              data = multilayer_data_runif,
#'              FUN = "p_value",
#'              normalize = "degree",
#'              crosstalk_params = crosstalk_params,
#'              seed_weights = seed_weights)
#' }
#' 
#' 
#' @export
#' 
module_identification <- function(network_layers, bipartite_networks = NULL, network_hierarchy = NULL, n = NULL, data = NULL, brw_attr = NULL,
                                  FUN = NULL, FUN_params = NULL, directed = FALSE, aggregate_layers = NULL,
                                  normalize = c("degree", "modified_degree"), k = 0.5, degree_bias = NULL,
                                  crosstalk_params = NULL, seed_weights = NULL, verbose = FALSE, 
                                  in_parallel = FALSE, n_cores = NULL, eta = NULL){
  start_time <- Sys.time()
  
  #=== Function Settings ===#
  ma_window <- 5 # Moving average window... set to Inf for a running average
  net_diam_prop <- -1 # Network diameter proportion to use for calculating betweenness and harmonic centrality. Negative value uses entire graph
  avg_rd <- 0.05 # initial average rate difference
  
  normalize <- match.arg(normalize)
  
  # Creating integrated graph
  g_int <- create_integrated_network(network_layers = network_layers, 
                                     bipartite_networks = bipartite_networks,
                                     network_hierarchy = network_hierarchy,
                                     data = data,
                                     FUN = FUN,
                                     FUN_params = FUN_params,
                                     directed = directed,
                                     brw_attr = brw_attr,
                                     lcc = TRUE)
  graph <- g_int$network
  network_hierarchy <- g_int$hierarchy
  
  # Ensuring that 'directed' object accurately reflects the graph object
  if(igraph::is_directed(graph)) {
    directed <- TRUE
  }
  
  # Flag for whether Biased Random Walk value is calculated by inverse distance from set of nodes of interest (NOI)
  BRW_NOI_FLAG <- is.character(brw_attr) && length(brw_attr) > 1
  
  if(is.null(n)) {
    n <- 10
    n_NULL_FLAG <- TRUE
  } else n_NULL_FLAG <- FALSE
  
  #========================#
  # Create aggregated graph
  #========================#
  # User will supply sets of layers, denoted by the corresponding leaf node names in hierarchy, which will be aggregated
  # Accepted input: NULL, list, or list of lists
  if(is.null(aggregate_layers)) {
    agg_graph <- create_aggregated_graph(graph = graph, agg_sets = NULL)
    agg_fun <- NULL
    V(graph)$agg_name <- V(graph)$name
  } else if(is.list(aggregate_layers)) {
    if(!"method" %in% names(aggregate_layers)) aggregate_layers$method <- "mean"
    
    if(length(aggregate_layers$method) > 1) aggregate_layers$method <- aggregate_layers$method[1]
    
    if(!any(aggregate_layers$method %in% c("mean", "median", "sum", "gmean", "hmean"))) {
      aggregate_layers$method <- "mean"
      warning("Unrecognized aggregation function name. Using mean.\nMust be one of 'mean', 'median', 'sum', 'gmean', 'hmean'.")
    }
    
    agg_fun <- aggregate_layers$method
    aggregate_layers <- aggregate_layers[names(aggregate_layers) != "method"]
    
    if(any(!unlist(aggregate_layers) %in% get_leaf_nodes(network_hierarchy, level = 0)[[1]])) stop("Incorrect primary names for layer aggregation. Must be leaves of network_hierarchy.")
    
    aggregate_layers <- lapply(aggregate_layers, unique)
    aggregate_layers <- aggregate_layers[unlist(lapply(aggregate_layers, length)) > 1]
    names(aggregate_layers) <- paste0("agg", seq_along(aggregate_layers))
    
    # Verifying that sets of layers to aggregate have no overlap
    tmp <- NULL
    for(i in seq_along(aggregate_layers)) {
      if(i == 1) {
        tmp <- c(tmp, aggregate_layers[[i]])
      } else {
        if(any(aggregate_layers[[i]] %in% tmp)) stop("For 'aggregate_layers', different sets of layers to aggregate cannot be overlapping.")
        tmp <- c(tmp, aggregate_layers[[i]])
      }
    }
    
    if(length(aggregate_layers) > 0) {
      agg_graph <- create_aggregated_graph(graph = graph, agg_sets = aggregate_layers, agg_fun = agg_fun, directed = directed)
      
      # Create new vertex attr that contains info of aggregated layers
      for(i in seq_along(aggregate_layers)) {
        tmp_nm <- paste(extract_string(V(graph)$name, "\\|", 1), names(aggregate_layers)[i], sep = "|")
        tmp_id <- which(V(graph)$layer %in% aggregate_layers[[i]])
        igraph::vertex_attr(graph, "agg_name", tmp_id) <- tmp_nm[tmp_id]
        tmp_id <- which(is.na(V(graph)$agg_name))
        igraph::vertex_attr(graph, "agg_name", tmp_id) <- V(graph)$name[tmp_id]
      }
    } else {
      agg_graph <- create_aggregated_graph(graph = graph, agg_sets = NULL)
      agg_fun <- NULL
    }
  } else stop("Unrecognized input type for aggregate_layers.")
  
  #======#
  # AMEND
  #======#
  repeat {
    all_scores <- list()
    all_nets <- list()
    
    subg <- vector(mode = "list", length = 2)
    subg[[1]] <- graph
    sub_ag <- vector("list", length = 2)
    sub_ag[[1]] <- agg_graph
    net_hier <- network_hierarchy
    
    iter_num <- 1
    
    if(is.null(eta)) eta <- get_eta0(agg_graph, n, "log_ratio")
    
    if(verbose) message(paste("Starting filtering rate:", round(eta, 4)))
    repeat {
      # if(iter_num == 14) stop("TEST")
      if(verbose) message(paste("Iteration:", iter_num))
      if(iter_num == 1) {
        rate_diff <- avg_rd
        e <- eta
        decay_N <- vcount(agg_graph)
      } else {
        rate_diff <- c(avg_rd, unlist(lapply(all_scores, function(x) x["Observed filtering rate"] - x["Filtering rate"])))
        rate_diff <- rate_diff[ifelse(length(rate_diff) >= ma_window, length(rate_diff) - ma_window + 1, 1):length(rate_diff)]
        e <- all_scores[[iter_num - 1]]["Filtering rate"]
        
        if(iter_num > 2) {
          decay_N <- length(all_nets[[iter_num - 2]])
        }
      }
      # Setting shifting percentile (i.e., filtering rate) for this iteration
      decay <- find_decay(N = decay_N, eta0 = e, n = n, rate_diff = rate_diff)
      l <- exp_filtering_rate(eta0 = e, d = decay)
      ll <- l[ifelse(iter_num == 1, 1, 2)]
      
      # Choosing Restart parameter value through a grid search
      rgs <- restart_grid_search(orig_net = subg[[1]], 
                                 agg_net = sub_ag[[1]], 
                                 network_hierarchy = net_hier, 
                                 normalize = normalize, 
                                 k = k, 
                                 crosstalk_params = crosstalk_params, 
                                 seed_weights = seed_weights, 
                                 degree_bias = degree_bias, 
                                 filtering_rate = ll, 
                                 agg_method = agg_fun, 
                                 in_parallel = in_parallel, 
                                 n_cores = n_cores, 
                                 net_diam_prop = net_diam_prop, 
                                 iteration = iter_num)
      
      subg[[2]] <- expand_graph(ig = subg[[1]], ag = rgs[[1]])
      sub_ag[[2]] <- rgs[[1]]
      
      # Update brw attr for full graph
      subg[[2]] <- set_brw_attr(graph = subg[[2]], brw_attr = brw_attr, BRW_NOI = BRW_NOI_FLAG)
      
      # break out of repeat loop if there is no change in network or if subnet size <= 2
      if(vcount(sub_ag[[2]]) == vcount(sub_ag[[1]]) || vcount(sub_ag[[2]]) <= 2) break
      
      mCS <- rgs[[2]][1] # Mean aggregate connectivity score of aggregated graph
      mZ <- rgs[[2]][2] # Mean standardized experimental scores of aggregated graph
      n_nodes <- rgs[[2]][3] # Size of subnetwork from aggregated graph
      
      sn_score <- mCS * mZ # subnetwork score for aggregated graph
      # wZ <- 0.75
      # sn_score <- (1 - wZ) * mCS + (wZ) * mZ
      
      all_scores[[length(all_scores) + 1]] <- c(rgs[[2]][4], sn_score, mZ, mCS, 
                                                n_nodes, igraph::ecount(sub_ag[[2]]), igraph::edge_density(sub_ag[[2]]), ll, mean(!V(sub_ag[[1]])$name %in% V(sub_ag[[2]])$name), decay)
      names(all_scores[[length(all_scores)]]) <- c("Restart parameter", "Network score", "Avg seed Z-score", "Avg connectivity score", 
                                                   "Node count", "Edge count", "Density", "Filtering rate", "Observed filtering rate", "Decay")
      all_nets[[length(all_nets) + 1]] <- V(sub_ag[[2]])$name
      names(all_nets[[length(all_nets)]]) <- V(sub_ag[[2]])$agg_layer
      
      # Updating network hierarchy object 'net_hier'
      rm_layers <- setdiff(unique(V(subg[[1]])$layer), unique(V(subg[[2]])$layer))
      for(i in seq_along(rm_layers)) {
        orig_classes <- class(net_hier)
        # Get the closest node in hierarchy that has other downstream leaf nodes
        sp <- igraph::shortest_paths(graph = net_hier, from = which.min(V(net_hier)$level), to = which(V(net_hier)$name == rm_layers[i]), mode = "out")
        sp <- as.numeric(sp$vpath[[1]])
        sp <- sp[order(V(net_hier)$level[sp], decreasing = TRUE)]
        for(j in seq_along(sp)) {
          tmp <- setdiff(get_leaf_nodes(net_hier, node = sp[j]), rm_layers[i])
          if(length(tmp) > 0) {
            rm_id <- sp[1:(j-1)] 
            break
          }
        }
        net_hier <- igraph::delete_vertices(net_hier, rm_id)
        class(net_hier) <- orig_classes
      }
      
      iter_num <- iter_num + 1
      subg[[1]] <- subg[[2]]
      sub_ag[[1]] <- sub_ag[[2]]
      
    }
    all_scores <- do.call("rbind", all_scores)
    all_scores <- round(as.data.frame(all_scores), 3)
    
    if(n_NULL_FLAG) {
      size_cond <- rep(TRUE, nrow(all_scores))
    } else {
      # Currently, final module cannot be larger than twice the user-specified 'n'
      size_col <- "Node count"
      size_window <- 2
      size_cond <- all_scores[, size_col] <= size_window * n & all_scores[, size_col] >= (1 / size_window) * n
      if(any(size_cond)) {
        break
      } else {
        if(verbose) message(paste0("Algorithm didn't converge with eta=", round(eta, 4), ". Trying a larger starting filtering rate."))
        if(eta == 1) stop("Algorithm is unable to converge. Try a larger 'n'.")
        eta <- eta * 1.1 # increase by 10%
        if(eta > 1) eta <- 1
      }
    }
  }
  
  # Choosing subnetwork with maximum score and (in the case of ties) closest to final module size 'n'
  # and in case of two networks with same max score and same distance from 'n', choose larger subnetwork
  score_col <- "Network score"
  if(sum(all_scores[size_cond, score_col] == max(all_scores[size_cond, score_col])) > 1 && !n_NULL_FLAG) { # Handling ties
    d <- abs(all_scores[,size_col] - n)
    ids <- which(all_scores[,score_col] == max(all_scores[size_cond, score_col]) & size_cond)
    max_id <- which(all_scores[,score_col] == max(all_scores[size_cond, ]) & d == min(d[ids]) & size_cond)[1]
  } else {
    max_id <- which(all_scores[,score_col] == max(all_scores[size_cond, score_col]) & size_cond)[1]
  }
  
  best_sn <- igraph::induced_subgraph(agg_graph, which(V(agg_graph)$name %in% all_nets[[max_id]]))
  best_score <- all_scores[max_id, score_col]
  
  time <- Sys.time() - start_time
  
  message("*** Converged! ***")
  
  input_params <- list(n = n, normalize = normalize, k = k, FUN = FUN, FUN_params = FUN_params, directed = directed, 
                       crosstalk_params = crosstalk_params, seed_weights = seed_weights, degree_bias = degree_bias, aggregate_layers)
  
  return(list(module = best_sn, score = best_score, subnetworks = all_nets, stats = all_scores, time = time, input_params = input_params))
}


#' @inherit module_identification
#' @export
AMEND <- function(...) {
  module_identification(...)
}


#' @title Get a larger subnetwork from previous iterations of the AMEND algorithm
#'
#' @description
#' This function allows the user to examine intermediate subnetworks that were generated during `module_identification()`, since the algorithm may overshoot the approximate final module size, or the user may be interested in the parent subnetworks of the final module.
#' 
#' @details
#' Since __AMEND__ iteratively filters out nodes to arrive at a final module, there are intermediate subnetworks associated with each iteration. A subnetwork at an earlier iteration contains the subnetworks of all later iterations. This function can be used to investigate these earlier subnetworks, for example, to identify when in the algorithm certain nodes were removed.
#' 
#' @param ig Input network
#' @param amend_object Result from `module_identification()`
#' @param k Iteration associated with the subnetwork you want to retrieve
#'
#' @returns igraph object of the subnetwork associated with the given iteration.
#'
#' @examples
#' # Attach igraph library
#' library(igraph)
#' 
#' # Perform module identification 
#' mod <- module_identification(network_layers = simple_network, n = 25)
#' 
#' # Get the processed igraph of the input network
#' g <- create_integrated_network(network_layers = simple_network)
#' 
#' # Get the subnetwork from the end of the 4th iteration of the AMEND algorithm
#' get_subnetwork(ig = g$network, amend_object = mod, k = 4)
#' 
#' @export
#' 
get_subnetwork <- function(ig, amend_object, k) {
  nodes <- amend_object$subnetworks[[k]]
  igraph::induced_subgraph(ig, which(V(ig)$name %in% nodes))
}


#' @title Calculate core-clustering coefficients
#'
#' @description `core_cc()` calculates the core-clustering coefficients of all nodes in a graph. This network concept was introduced by _Bader & Hogue_ (2003).
#'
#' @details The core-clustering coefficient of a node is the edge density of the maximum k-core of the immediate neighborhood of that node. The k-core of a graph is the maximal subgraph in which every node has degree >= k.
#'
#' @param g igraph. Must be connected
#' @param weight Logical. If true, the edge density is calculated taking into account edge weights
#'
#' @returns Vector of core-clustering coefficients for each node
#'
#' @examples
#' library(igraph)
#' 
#' graph <- sample_pa(n = 100, power = 1.2)
#' core_cc(graph)
#'
#' @export
#'
core_cc <- function(g, weight = TRUE) {
  if(!"weight" %in% igraph::edge_attr_names(g) && weight) {
    weight <- FALSE
  }
  n <- vcount(g)
  w <- numeric(n)
  for(i in seq_len(n)) {
    nborhood <- igraph::make_ego_graph(g, order = 1, nodes = i, mode = "all")[[1]]
    cores <- igraph::coreness(nborhood, mode = "all")
    max_core <- igraph::induced_subgraph(nborhood, which(cores == max(cores)))
    w[i] <- edge_density_weighted(max_core, weight = weight)
  }
  return(w)
}


#' @title Calculate a node-wise connectivity score
#' 
#' @description
#' Calculates a centrality or several centrality metrics and converts to empirical cumulative probabilities. inverse=TRUE multiplies values by -1 before ECDF calculation.
#' 
#' @details
#' For mode="core", only considers the coreness of nodes when calculating ECDF. For mode="mix", it considers a mix of betweenness centrality, harmonic centrality, and core clustering coefficient, aggregated by applying ECDF to each score individually then taking a weighted average of these probabilities.
#'
#' @param graph igraph
#' @param net_diam_prop numeric. Proportion of graph to use when calculating betweenness and harmonic centrality, for mode="mix" only.
#' @param inverse Logical. Whether to return values inversely related to connectivity measure selected in __mode__.
#' @param mode character scalar. One of "core" or "mix"
#' 
#' @returns numeric vector
#' 
#' @examples
#' net <- create_integrated_network(network_layers = simple_network, lcc = TRUE)$network
#' ncs <- node_connectivity_score(graph = net, inverse = FALSE, mode = "core")
#' ncs_inv <- node_connectivity_score(graph = net, inverse = TRUE, mode = "core")
#' node_coreness <- core_cc(g = net)
#' cor(node_coreness, ncs)
#' cor(node_coreness, ncs_inv)
#' 
#' @export
#'
node_connectivity_score <- function(graph, net_diam_prop = -1, inverse = FALSE, mode = c("core", "mix")) {
  
  mode <- match.arg(mode)
  
  # net_diam_prop
  if(net_diam_prop > 1) stop("'net_diam_prop' must be negative or between 0 and 1.")
  
  # Centrality measures
  if(net_diam_prop < 0) {
    kpath <- net_diam_prop
  } else {
    kpath <- net_diam_prop * igraph::diameter(graph = graph, directed = TRUE)
  }
  
  if(mode == "mix") {
    c_val <- data.frame(ccc = core_cc(graph, weight = TRUE),
                        bc = igraph::betweenness(graph = graph, directed = TRUE, normalized = TRUE, cutoff = kpath),
                        hc = igraph::harmonic_centrality(graph = graph, mode = 'out', normalized = TRUE, cutoff = kpath),
                        row.names = V(graph)$name)
  } else if(mode == "core") {
    c_val <- data.frame(core = igraph::coreness(graph = graph, mode = "all"),
                        row.names = V(graph)$name)
  }
  
  
  # Calculate eCDF & Calculate empirical cumulative probabilities
  k <- ifelse(inverse, -1, 1)
  for(j in seq_len(ncol(c_val))) {
    FUN <- stats::ecdf(k * c_val[,j])
    c_val[,j] <- FUN(k * c_val[,j])
  }
  agg_cent <- round(apply(c_val, 1, mean), 4)
  
  return(agg_cent)
}


#' @title Create an aggregated graph
#'
#' @description
#' An aggregated graph, in a multiplex context, is when the layers of a multiplex component are collapsed such that there is only one layer. A user-specified layer is designated as the primary, and edges in this layer take priority over edges of other layers. First, all edges having both of its adjacent nodes in the primary layer are included in this aggregated component. Then, iterating through the remaining layers, nodes and edges are included if the two adjacent nodes of an edge are not both already in the aggregated component. The collapsed multiplex will necessarily contain all of the original nodes, as well as all of the edges in the primary layer.
#'
#' The purpose is to aggregate a multiplex component before the subnetwork identification step (using `solve_MWCS()`) so that the algorithm preferentially considers the edges of the primary layer. This process involves aggregating node attribute values and fusing nodes together such that the fused node is adjacent to all of the edges of its predecessors.
#'
#' @param graph igraph object. The multiplex graph to be aggregated.
#' @param control A named list. The element 'primary' contains the name of the primary layer for a multiplex component to be used during subnetwork identification. The element 'agg.method' contains a character scalar referring to an aggregation function.
#' @param directed logical. TRUE for directed graphs.
#'
#' @return igraph object
#'
#' @seealso [melt_graph()], [expand_graph()], [module_identification()]
#'
#' @noRd
#' 
create_aggregated_graph <- function(graph, agg_sets, agg_fun, directed) {
  
  original_classes <- class(graph)
  
  V(graph)$agg_layer <- V(graph)$layer
  
  for(i in seq_along(agg_sets)) {
    all_nodes <- extract_string(V(graph)$name, "\\|", 1)
    all_layers <- extract_string(V(graph)$name, "\\|", 2)
    
    primary_layer <- agg_sets[[i]][1]
    other_layers <- agg_sets[[i]][-1]
    new_layer_name <- names(agg_sets)[i]
    
    primary_ids <- which(V(graph)$layer == primary_layer)
    
    # Strip layer info from node labels in primary layer
    v_targets <- all_nodes[primary_ids]
    V(graph)$name[primary_ids] <- paste(v_targets, new_layer_name, sep = "|")
    V(graph)$agg_layer[primary_ids] <- new_layer_name
    
    for(j in seq_along(other_layers)) {
      other_flag <- all_layers == other_layers[j]
      # Identify primary and non-primary nodes
      ## primary
      pn <- which(all_nodes %in% v_targets & other_flag)
      v_pn <- V(graph)$name[pn]
      ## non-primary
      npn <- which(!all_nodes %in% v_targets & other_flag)
      v_npn <- V(graph)$name[npn]
      ## Identify and delete primary-primary edges
      # Identify
      el <- igraph::as_edgelist(graph)
      pp_id <- which(apply(el, 1, function(x) all(x %in% v_pn)))
      # Delete
      graph <- igraph::delete_edges(graph = graph, edges = pp_id)
      # Strip layer information from node labels
      V(graph)$name[other_flag] <- paste(extract_string(V(graph)$name[other_flag], "\\|", 1), new_layer_name, sep = "|")
    }
    # Melt graph, i.e., fuse common nodes such that only one of them becomes adjacent to the edges of the others, and then remove duplicate nodes
    graph <- melt_graph(g = graph, agg_fun = agg_fun, directed = directed)
  }
  class(graph) <- original_classes
  return(graph)
}


#' @title Melt a graph
#'
#' @description
#' Melting a graph is finding nodes with the same node name and fusing them together so that there is just one node that is adjacent to all of their edges.
#'
#' @details
#' Fuses any nodes with same node name together, conserving any vertex attributes that may be present.
#'
#' @param g igraph
#' @param agg.method A function or NULL. The aggregation method for vertex attributes. Default NULL takes the first value.
#' @param directed logical. TRUE for directed graphs.
#'
#' @returns igraph with vertex and edge attributes preserved.
#' 
#' @noRd
#'
melt_graph <- function(g, agg_fun = NULL, directed) {
  # Get names of all vertex atrributes in igraph
  v_attr_names <- igraph::vertex_attr_names(g)
  # Collect any vertex attributes from graphs
  v_attrs <- vector("list", length(v_attr_names)); names(v_attrs) <- v_attr_names
  for(j in seq_along(v_attrs)) {
    v_attrs[[j]] <- igraph::vertex_attr(g, v_attr_names[j])
    names(v_attrs[[j]]) <- V(g)$name
  }
  
  # Get names of all edge atrributes in igraphs
  e_attr_names <- igraph::edge_attr_names(g)
  # Collect all edge attributes
  e_attrs <- vector("list", length(e_attr_names)); names(e_attrs) <- e_attr_names
  for(j in seq_along(e_attrs)) {
    e_attrs[[j]] <- igraph::edge_attr(g, names(e_attrs)[j])
  }
  
  
  edge_list <- igraph::as_edgelist(graph = g, names = TRUE)
  # if("weight" %in% igraph::edge_attr_names(g)) {
  #   edge_list <-  cbind(edge_list, igraph::E(g)$weight)
  # } else edge_list <- cbind(edge_list, 1)
  graph <- igraph::graph_from_edgelist(el = edge_list, directed = directed)
  # igraph::E(graph)$weight <- as.numeric(edge_list[,3])
  
  # Assign edge attributes 
  for(i in seq_along(e_attrs)) {
    if(names(e_attrs)[i] == "weight") {
      igraph::edge_attr(graph, names(e_attrs)[i]) <- as.numeric(e_attrs[[i]])
    } else {
      igraph::edge_attr(graph, names(e_attrs)[i]) <- e_attrs[[i]]
    }
  }
  
  # Redistribute the vertex attributes to the integrated graph. Takes the first vertex attribute value for nodes with multiple values.
  if(is.null(agg_fun)) {
    for(i in seq_along(v_attrs)) {
      igraph::vertex_attr(graph, names(v_attrs)[i]) <- unname(v_attrs[[i]][match(V(graph)$name, names(v_attrs[[i]]))])
    }
  } else {
    for(i in seq_along(v_attrs)) {
      if(is.numeric(v_attrs[[i]])) {
        if(any(v_attrs[[i]] < 0) && agg_fun %in% c("gmean", "hmean")) {
          agg_fun_tmp <- mean
        } else agg_fun_tmp <- get_aggregate_method(agg_fun)
        tmp <- stats::aggregate(x = v_attrs[[i]], by = list(names(v_attrs[[i]])), FUN = agg_fun_tmp)
        igraph::vertex_attr(graph, names(v_attrs)[i]) <- tmp$x[match(V(graph)$name, tmp$Group.1)]
      } else { # If non-numeric, take the first value
        igraph::vertex_attr(graph, names(v_attrs)[i]) <- unname(v_attrs[[i]][match(V(graph)$name, names(v_attrs[[i]]))])
      }
    }
  }
  
  # Simplify graph
  graph <- igraph::simplify(graph = graph, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = list(weight = "median", "first"))
  
  return(graph)
}


#' @title Grid search for the restart parameter in RWR
#'
#' @description
#' Performs a grid search for the restart parameter in RWR. For use inside `run_AMEND()`. Can be run in parallel if specified.
#'
#' @details
#' This is a simple grid search for the restart probability parameter in RWR. For a grid of restart probability values, `RWR()` is run, raw scores are shifted by the current filtering rate of that iteration, then `solve_MWCS()` is run to find optimal subnetworks, which are scored. The restart probability resulting in the largest subnetwork score is used.
#'
#' @return a list containing an igraph object, subnetwork score, and restart value
#'
#' @noRd
#'
restart_grid_search <- function(orig_net, agg_net, network_hierarchy, normalize, k, crosstalk_params, seed_weights, degree_bias, filtering_rate, agg_method = NULL, in_parallel = FALSE, n_cores, net_diam_prop, iteration = 1) {
  #=== Function Settings ===#
  FILTERING_RATE_DIFF <- 0.25
  GRID_NET_SIZE <- 30000
  #=========================#
  
  # Normalize adjacency matrix to get transition matrix
  tmat <- transition_matrix(network = orig_net,
                            network_hierarchy = network_hierarchy,
                            normalize = normalize,
                            k = k,
                            crosstalk_params = crosstalk_params,
                            degree_bias = degree_bias,
                            in_parallel = in_parallel,
                            n_cores = n_cores)
  
  # Create seed vector
  seeds <- matrix(V(orig_net)$seeds, ncol = 1, dimnames = list(V(orig_net)$name))
  
  if(!is.null(agg_method)) {
    agg_fun <- get_aggregate_method(agg_method)
    agg_flag <- TRUE
  } else agg_flag <- FALSE
  
  if(iteration == 1) {
    if(igraph::vcount(agg_net) > GRID_NET_SIZE) {
      grid <- seq(0.5, 0.95, by = 0.1)
    } else grid <- seq(0.5, 0.95, by = 0.05)
  } else {
    if(igraph::vcount(agg_net) > GRID_NET_SIZE) {
      grid <- seq(0.1, 0.95, by = 0.1)
    } else grid <- seq(0.1, 0.95, by = 0.05)
  }
  
  ncs <- node_connectivity_score(graph = orig_net, net_diam_prop = net_diam_prop, inverse = TRUE, mode = "core")
  
  if(!in_parallel) {
    score_metrics <- matrix(nrow = length(grid), ncol = 3)
    nets <- vector(mode = "list", length = length(grid))
    rm_id <- NULL
    for(i in seq_along(grid)) {
      # RWR
      rwr_score <- RWR(tmat = tmat, seeds = seeds, restart = grid[i], seed_weights = seed_weights, 
                       network_hierarchy = network_hierarchy, node_specific_restart = TRUE, node_connectivity_scores = ncs)
      rwr_score <- rwr_score[,1]
      names(rwr_score) <- V(orig_net)$name
      
      # Aggregate RWR scores
      if(agg_flag) {
        tmp_val <- stats::aggregate(rwr_score, by = list(V(orig_net)$agg_name), FUN = agg_fun)
        rwr_score <- tmp_val$x
        names(rwr_score) <- tmp_val$Group.1
      }
      
      # Maximum Scoring Subgraph Algorithm
      rwr_score <- rwr_score - quantile(rwr_score, filtering_rate) # if filtering rate close to 0, quantile takes the minimum of rwr_score
      nets[[i]] <- solve_MWCS(agg_net, rwr_score)
      
      nv <- vcount(nets[[i]])
      mZ <- sum(V(nets[[i]])$Z) / nv
      
      # To prevent filtering out too many nodes at one iteration
      # This is the difference between the observed and theoretical filtering rates
      if((vcount(agg_net) - nv) / vcount(agg_net) - filtering_rate >= FILTERING_RATE_DIFF) {
        rm_id <- c(rm_id, i)
        next
      }
      
      score_metrics[i,] <- c(mZ, nv, grid[i])
    }
  }
  if(in_parallel) {
    `%dopar%` <- get("%dopar%", asNamespace("foreach"))
    if(is.null(n_cores)) n_cores <- round(parallel::detectCores() * 0.66)
    cl <- parallel::makeForkCluster(n_cores, outfile = "")
    on.exit(parallel::stopCluster(cl))
    if(!foreach::getDoParRegistered()) doParallel::registerDoParallel(cl)
    res <- foreach::foreach(i = seq_along(grid), .verbose = FALSE, .packages = c("igraph", "Matrix"), .export = NULL, .noexport = NULL) %dopar% {
      # RWR
      rwr_score <- RWR(tmat = tmat, seeds = seeds, restart = grid[i], seed_weights = seed_weights, 
                       network_hierarchy = network_hierarchy, node_specific_restart = TRUE, node_connectivity_scores = ncs)
      rwr_score <- rwr_score[,1]
      names(rwr_score) <- V(orig_net)$name
      
      # Aggregate RWR scores
      if(agg_flag) {
        tmp_val <- stats::aggregate(rwr_score, by = list(V(orig_net)$agg_name), FUN = agg_fun)
        rwr_score <- tmp_val$x
        names(rwr_score) <- tmp_val$Group.1
      }
      
      # Maximum Scoring Subgraph Algorithm
      rwr_score <- rwr_score - quantile(rwr_score, filtering_rate) # if filtering rate close to 0, quantile takes the minimum of rwr_score
      nets <- solve_MWCS(agg_net, rwr_score)
      
      nv <- vcount(nets)
      mZ <- sum(V(nets)$Z) / nv
      
      rm_id <- NULL
      
      # To prevent filtering out too many nodes at one iteration
      # This is the difference between the observed and theoretical filtering rates
      if((vcount(agg_net) - nv) / vcount(agg_net) - filtering_rate >= FILTERING_RATE_DIFF){
        rm_id <- i
        nv <- 0
        mZ <- 0
      }
      
      list(c(mZ, nv, grid[i]), nets, rm_id)
    }
    # parallel::stopCluster(cl)
    score_metrics <- do.call(rbind, lapply(res, function(x) x[[1]]))
    nets <- lapply(res, function(x) x[[2]])
    rm_id <- unlist(lapply(res, function(x) x[[3]]))
  }
  
  rm_id <- unique(c(rm_id, which(score_metrics[,1] %in% c(0, NA))))
  if(length(rm_id) != 0) {
    score_metrics <- score_metrics[-rm_id,]
    nets <- nets[-rm_id]
  }
  
  mCS <- network_connectivity_score(graphs = nets, net_diam_prop = net_diam_prop, in_parallel = in_parallel, n_cores = n_cores)
  score_metrics <- cbind(mCS, score_metrics)
  
  if(isTRUE(nrow(score_metrics) >= 1)) { # if score matrix isn't empty
    max_id <- which.max(score_metrics[, 2] * score_metrics[, 1])
    new_scores <- score_metrics[max_id,]
    new_g <- nets[[max_id]]
  } else {
    new_g <- agg_net
    nv <- vcount(new_g)
    avg_z <- sum(V(new_g)$Z) / nv
    new_scores <- c(network_connectivity_score(graphs = list(new_g), net_diam_prop = net_diam_prop), avg_z, nv, NA)
  }
  return(list(new_g, new_scores))
}


#' @title Aggregated connectivity score for graphs
#'
#' @description
#' Calculate aggregated connectivity scores for a set of two or more graphs based on harmonic centrality, betweenness centrality, and core-clustering coefficient of their constituent nodes. The aggregated score is relative to the input set of graphs.
#'
#' @details
#' Given a list of subgraphs coming from a common parent graph, this function calculates harmonic centrality, betweenness centrality, and core-clustering coefficients for all nodes in each subgraph. The complete set of values for each centrality measure is used to create empirical cumulative distribution functions, which are then used to assign cumulative probabilities to each node. An aggregated score is calculated for each node by taking a weighted mean of their cumulative probabilities of the 3 centrality measures. A subgraph's aggregated connectivity score is the arithmetic mean of the node-wise aggregate scores.
#'
#' @param graphs A list of subgraphs coming from a common parent graph
#' @param net_diam_prop The proportion of a graphs diameter used to calculate the maximum path length (_net_diam_prop_ multiplied by graph diameter) to consider when calculating the betweenness. Can be used to speed up calculations. If zero or negative there is no such limit.
#' @param in_parallel Logical. Whether to run certain operations in parallel, using the _parallel_, _doParallel_ and _foreach_ packages.
#' @param n_cores Numeric scalar or NULL (default). Number of cores to use during parallel processing. When NULL and in_parallel=TRUE, defaults to two-thirds the number of cores detected on machine.
#' 
#' @returns numeric vector of aggregated connectivity scores
#' 
#' @noRd
#'
network_connectivity_score = function(graphs, net_diam_prop = -1, in_parallel = FALSE, n_cores = NULL){
  ## Input checks
  # graphs
  if(length(graphs) == 0) {
    return(numeric(0))
  } 
  
  # net_diam_prop
  if(net_diam_prop > 1) stop("'net_diam_prop' must be negative or between 0 and 1.")
  
  # Calculating centrality measures
  if(in_parallel) {
    `%dopar%` <- get("%dopar%", asNamespace("foreach"))
    if(is.null(n_cores)) n_cores <- round(parallel::detectCores() * 0.66)
    cl <- parallel::makeForkCluster(n_cores, outfile = "")
    on.exit(parallel::stopCluster(cl))
    if(!foreach::getDoParRegistered()) doParallel::registerDoParallel(cl)
    c_val <- foreach::foreach(i = seq_along(graphs), .verbose = FALSE, .packages = c("igraph"), .export = NULL, .noexport = NULL) %dopar% {
      # Centrality measures
      if(net_diam_prop < 0) {
        kpath <- net_diam_prop
      } else {
        kpath <- net_diam_prop * igraph::diameter(graph = graphs[[i]], directed = TRUE)
      }
      ccc <- core_cc(graphs[[i]], weight = TRUE)
      hc <- igraph::harmonic_centrality(graph = graphs[[i]], mode = 'out', normalized = TRUE, cutoff = kpath)
      bc <- igraph::betweenness(graph = graphs[[i]], directed = TRUE, normalized = TRUE, cutoff = kpath)
      data.frame(ccc = ccc, bc = bc, hc = hc, row.names = V(graphs[[i]])$name)
    }
    # parallel::stopCluster(cl)
  } else {
    c_val <- vector("list", length(graphs))
    for(i in seq_along(c_val)) {
      # Centrality measures
      if(net_diam_prop < 0) {
        kpath <- net_diam_prop
      } else {
        kpath <- net_diam_prop * igraph::diameter(graph = graphs[[i]], directed = TRUE)
      }
      ccc <- core_cc(graphs[[i]], weight = TRUE)
      hc <- igraph::harmonic_centrality(graph = graphs[[i]], mode = 'out', normalized = TRUE, cutoff = kpath)
      bc <- igraph::betweenness(graph = graphs[[i]], directed = TRUE, normalized = TRUE, cutoff = kpath)
      c_val[[i]] <- data.frame(ccc = ccc, bc = bc, hc = hc, row.names = V(graphs[[i]])$name)
    }
  }
  
  # Calculate eCDF
  cdf <- vector("list", ncol(c_val[[1]]))
  names(cdf) <- names(c_val[[1]])
  for(i in seq_along(cdf)) {
    vals <- unlist(lapply(c_val, function(x) x[,names(cdf)[i]]))
    cdf[[i]] <- stats::ecdf(vals)
  }
  # Calculate empirical cumulative probabilities
  agg_cent <- numeric(length(c_val))
  for(i in seq_along(agg_cent)) {
    tmp <- vector("list", length(cdf))
    for(j in seq_len(ncol(c_val[[i]]))) {
      tmp[[j]] <- cdf[[colnames(c_val[[i]])[j]]](c_val[[i]][,j])
    }
    tmp <- do.call(data.frame, tmp)
    colnames(tmp) <- colnames(c_val[[i]])
    tmp_agg_cent <- apply(tmp, 1, max)
    agg_cent[i] <- mean(tmp_agg_cent)
  }
  
  return(agg_cent)
}


#' @title Set the starting filtering rate
#'
#' @description
#' Sets the starting filtering rate (denoted as eta, the initial quantile of RWR scores used to shift scores before use in `solve_MWCS()`) in `AMEND()`.
#'
#' eta is linearly increasing with D(N,n), where D is some distance metric between the input network size N (i.e., vcount(g)) and the final subnetwork size n.
#' For method = "difference", D(N,n) = N - n
#' For method = "ratio", D(N,n) = N / n
#' For method = "log_ratio", D(N,n) = log(N / n)
#' Below are the boundaries within which eta changes linearly with D(N,n). Set to capture "normal/expected" input network sizes and user-defined final module sizes.
#' N: (500, 15000), n: (5, 200), eta: (0.3, 0.8)
#'
#' @param g input graph
#' @param n Final module size to approximate
#' @param method One of 'difference', 'ratio', or 'log_ratio'
#' @param eta_range Range for eta within which eta changes linearly with D(N,n)
#' @param N_range Range for N (i.e., vcount(g)) within which eta changes linearly with D(N,n)
#' @param n_range Range for n within which eta changes linearly with D(N,n)
#'
#' @returns Starting filtering rate
#' 
#' @noRd
#'
get_eta0 <- function(g, n, method = c("log_ratio", "difference", "ratio"), eta_range = c(0.2, 0.8), N_range = c(500, 25000), n_range = c(5, 200)) {
  eta.l <- eta_range[1]
  eta.u <- eta_range[2]
  N.l <- min(N_range[1], vcount(g))
  N.u <- max(N_range[2], vcount(g))
  n.l <- min(n_range[1], n)
  n.u <- max(n_range[2], n)
  if(method == "difference") {
    D.l <- N.l - n.u
    D.u <- N.u - n.l
    D <- vcount(g) - n
  } else if(method == "ratio") {
    D.l <- N.l / n.u
    D.u <- N.u / n.l
    D <- vcount(g) / n
  } else if(method == "log_ratio") {
    D.l <- log(N.l / n.u)
    D.u <- log(N.u / n.l)
    D <- log(vcount(g) / n)
  } else stop("'method' should be one of 'difference', 'ratio', or 'log_ratio'")
  # Getting slope and intercept
  a <- (eta.u - eta.l) / (D.u - D.l) # slope
  b <- (eta.u * D.l - eta.l * D.u) / (D.l - D.u) # intercept
  
  eta0 <- max(0, min(1, a * D + b))
  return(eta0)
}


#' @title Exponential decay schedule for filtering rate
#'
#' @description Exponential decay schedule for filtering rate. Used to shift the raw RWR scores
#'
#' @param eta0 Starting filtering rate
#' @param d Decay parameter
#'
#' @returns vector of filtering rates for each iteration of AMEND
#'
#' @noRd
#'
exp_filtering_rate <- function(eta0 = 0.5, d = 0.1) {
  iterations <- 1:100
  return(eta0 * exp(-d * (iterations - 1)))
}


#' @title Sets the decay value
#'
#' @description Finds the max decay value for an exponential filtering rate schedule that will
#' allow the algorithm to arrive at a graph of size n.
#'
#' @param N Size of graph
#' @param eta0 Starting filtering rate
#' @param n Approximate size of final module
#' @param rate_diff Vector of differences between observed filtering rate and theoretical filtering rate
#'
#' @returns decay value, numeric
#'
#' @noRd
#'
find_decay <- function(N, eta0 = 0.5, n = 25, rate_diff = 0.05) {
  rate_diff <- mean(rate_diff)
  
  iters <- 100
  d_grid <- seq(0.01, 1, 0.01)
  decay <- NULL
  k_max <- 0.1
  k_min <- -0.1
  
  for(j in seq_along(d_grid)) {
    under <- FALSE
    perc <- exp_filtering_rate(eta0 = eta0, d = d_grid[j])
    
    d1 <- N
    for(i in seq_len(iters)) {
      # Logic: rate difference up, K up, AFR up, decay up, theoretical filtering rate down, larger networks
      # rate_diff = observed minus theoretical filtering rate, so more nodes filtered out than expected, so need to push on the brakes i.e., increase decay.
      K <- max(min(k_max, (rate_diff * (1 - perc[i])) / perc[i]), k_min) # Will lead to slightly lower adjusted filtering rate (afr) --> lower decay --> larger filtering rates --> smaller networks
      afr <- perc[i] * (1 + K) # Adjusted filtering rate (AFR)
      d2 <- round(d1 * (1 - afr), 0)
      if(d2 <= n) {
        under <- TRUE
        break
      }
      if(d1 == d2) {
        break
      }
      d1 <- d2
    }
    if(under) {
      decay <- c(decay, d_grid[j])
    } else break
  }
  
  if(length(decay) == 0) {
    return(min(d_grid))
  } else return(max(decay))
}


#' @title Calculate the edge density of a weighted graph
#'
#' @description `edge_density_weighted()` calculates the edge density of a graph, taking into account edge weights.
#'
#' @details If weight is true, the edge density equals the sum of the edge weights divided by the maximum number of edges possible for the given graph. If weight is false, the edge density is the number of edges divided by the maximum number of edges possible for the given graph.
#'
#' @param graph Input graph
#' @param weight Logical. If true, the edge density is calculated taking into account edge weights
#'
#' @returns weighted edge density of graph
#'
#' @noRd
#'
edge_density_weighted <- function(graph, weight) {
  if(weight) {
    ew <- igraph::E(graph)$weight
    if(any(ew > 1)) {
      warning("Some edge weights are greater than one. Computing unweighted edge density.")
      return(igraph::edge_density(graph))
    }
    N <- vcount(graph)
    tmp <- ifelse(igraph::is_directed(graph), 1, 2)
    return(sum(ew) / (N * (N - 1) / tmp))
  } else return(igraph::edge_density(graph))
}


#' @title Get a pre-defined aggregation function
#'
#' @description Given a character scalar, this returns a function that will aggregate values of common nodes.
#'
#' @param x character scalar
#'
#' @return aggregation function
#' 
#' @noRd
#'
get_aggregate_method <- function(x) {
  if(x == "mean") {
    agg_fun <- mean
  } else if(x == "median") {
    agg_fun <- stats::median
  } else if(x == "sum") {
    agg_fun <- sum
  } else if(x == "gmean") {
    agg_fun <- function(y) (prod(y)) ^ (1 / length(y))
  } else if(x == "hmean") {
    agg_fun <- function(y) length(y) / (sum(1 / y))
  }
  return(agg_fun)
}


#' @title Expand a graph
#'
#' @description
#' Given an aggregated graph, returns its expanded counterpart, induced from the original, full graph.
#' Takes nodes from aggregated graph, identifies nodes that were in a collapsed multiplex component, duplicates them as necessary, and appends appropriate '|component_layer' information.
#' Returns an induced subgraph of ig that contains only node names matching nodes of aggregated graph.
#'
#' @param ig Input graph
#' @param ag Aggregated graph
#'
#' @return igraph object
#' 
#' @noRd
#'
expand_graph <- function(ig, ag) {
  original_classes <- class(ig)
  ag_names <- V(ag)$name
  ig_names <- V(ig)$agg_name
  ig <- igraph::induced_subgraph(graph = ig, vids = which(ig_names %in% ag_names))
  class(ig) <- original_classes
  return(ig)
}

