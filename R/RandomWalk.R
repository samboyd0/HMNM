# NOTES:
# - Users should use the function create_integrated_network() to correctly format their network and hierarchy prior to transition_matrix

#=================================================================#

#=== IMPORTS ===#
#' @importFrom igraph V 
#' @importFrom foreach %dopar%


#' @title Pipeline for processing a mono-/multi-layer network and running Random Walk with Restart (RWR)
#' 
#' @description `RWR_pipeline()` takes as input a mono- or multi-layer network, constructs an transition matrix, and then performs RWR. Resulting node scores can be returned as a vector or as a list with elements containing scores from each layer.
#' 
#' @param network_layers,bipartite_networks Single network-like object (igraph, adjacency matrix, or edgelist) or a list of these. If a list, it must be named, with names matching category names in network_hierarchy. If multiple layers are contained in a single object, the list name must include these layer names separated by "|". bipartite_networks objects should contain the mappings between different layers. Elements in bipartite_networks list can be set to "common", which will connect all common nodes between the designated layers.
#' @param network_hierarchy A 2-column matrix representing an edgelist for the network hierarchy, where hierarchy vertices represent categories which categorize the network nodes. Or an object of class 'hierarchy' as a result from `create_network_hierarchy()`.
#' @param data Named list of numeric vectors, or a single numeric vector, containing numeric values from which seed values for RWR will be calculated (with _FUN_ and _FUN_params_), a character string, or NULL. Names of list should match layer names. Numeric values must be named with the corresponding node name. If a string, this should be the vertex attribute name (for igraph inputs) containing the data.
#' @param FUN Function, list of functions, or a character string denoting a default function ('binary', 'shift_scale', 'p_value', or 'exp'), to be applied to _data_ to compute seed values for RWR. Names of list must match layer names. NULL (default) applies no transformation of values in _data_. Optional function arguments given in _FUN_params_. 
#' @param FUN_params List or list of lists, containing additional function arguments for functions given in _FUN_. NULL (default) doesn't supply any additional function arguments.
#' @param directed Logical. Whether the input network should be treated as directed.
#' @param brw_attr Similar format as _data_. Contains values to be used in a biased random walk. Should contain non-negative values.
#' @param lcc Logical. Whether to take the largest connected component of the resulting network.
#' @param normalize Adjacency matrix normalization method to construct transition matrix.
#' @param k Penalization factor for normalize="modified_degree". Must be non-negative, with larger values resulting in a greater penalty for node degree, in an effort to mitigate degree bias.
#' @param crosstalk_params A named numeric vector containing the crosstalk parameters for each category in network hierarchy. If NULL (default), a uniform value of 0.5 is set. Hierarchicy categories not given in _crosstalk_params_ will be given this default value of 0.5.
#' @param degree_bias A character vector or list, or NULL (default). The character vector denotes the layers to which the degree bias mitigation method will be applied. The list must contain this character vector of layers (named 'layers') and a numeric scalar (named 'gamma') between 0 and 1 denoting the strength of degree bias mitigation. The default gamma value is 0.2.
#' @param restart restart probability for RWR
#' @param seed_weights A list of named numeric vectors, or NULL (default). List elements should correspond to sibling sets of categories in network hierarchy. Values in each set must sum to one. NULL gives uniform values within each sibling set.
#' @param output "list" for a list of RWR scores separated by network layer. "vector" for a single vector of all RWR scores.
#'
#' @return a numeric vector or a list of numeric vectors
#' 
#' @seealso [create_integrated_graph()], [create_network_hierarchy()], [transition_matrix()], [RWR()]
#' 
#' @export
#' 
RWR_pipeline <- function(network_layers, bipartite_networks = NULL, network_hierarchy = NULL, data = NULL, FUN = NULL, FUN_params = NULL, directed = FALSE, brw_attr = NULL, lcc = FALSE,
                         normalize = c("degree", "modified_degree"), k = 0.5, crosstalk_params = NULL, degree_bias = NULL, restart = 0.5, seed_weights = NULL, output = c("list", "vector")) {
  if(0){ # TEST
    norm = "modified_degree"
    k = 0.5
    crosstalk_params = c(protdat = 0.8, phosphdat = 0.4, meta = 0.2)
    degree_bias = list(layers = V(g$hierarchy)$name[V(g$hierarchy)$level > 0], gamma = 0.1)
    restart = 0.5
    seed_weights = NULL
    output = "list"
  }
  
  output <- match.arg(output)
  
  # Create integrated network
  g <- create_integrated_network(network_layers = network_layers,
                                 bipartite_networks = bipartite_networks,
                                 network_hierarchy = network_hierarchy,
                                 data = data,
                                 FUN = FUN,
                                 FUN_params = FUN_params,
                                 directed = directed,
                                 brw_attr = brw_attr,
                                 lcc = lcc)
  
  # Create transition matrix
  tmat <- transition_matrix(network = g$network,
                            network_hierarchy = g$hierarchy,
                            norm = normalize,
                            k = k,
                            crosstalk_params = crosstalk_params,
                            degree_bias = degree_bias)
  
  # Run RWR
  seed_vals <- V(g$network)$seeds
  names(seed_vals) <- V(g$network)$name
  ncs <- node_connectivity_score(graph = g$network, inverse = TRUE)
  r <- RWR(tmat = tmat, seeds = seed_vals, restart = restart, seed_weights = seed_weights, network_hierarchy = g$hierarchy, node_connectivity_scores = ncs)
  
  if(output == "list") {
    all_layers <- extract_string(names(r), "\\|", 2)
    uniq_layers <- unique(all_layers)
    res <- vector("list", length(uniq_layers)); names(res) <- uniq_layers
    for(i in seq_along(res)) {
      ids <- which(all_layers == uniq_layers[i])
      res[[i]] <- r[ids]
      names(res[[i]]) <- V(g$network)$original_name[ids]
    }
    return(res)
  } else {
    return(r)
  }
}

# TEST
if(0){
  source("/Users/samboyd/Documents/NHNM/R package/TEST/TEST_create_integrated_network.R")
  source("/Users/samboyd/Documents/NHNM/R package/NHNM/R/create_integrated_network.R")
  source("/Users/samboyd/Documents/NHNM/R package/NHNM/R/RandomWalk.R")
  
  library(igraph)
  library(Matrix)
  
  r <- RWR_pipeline(network_layers = network_layers, 
                    bipartite_networks = bipartite_networks,
                    network_hierarchy = network_hierarchy,
                    data = data,
                    FUN = FUN,
                    FUN_params = FUN_params,
                    directed = directed,
                    brw_attr = brw_attr,
                    lcc = lcc,
                    norm = "modified_degree",
                    k = 0.5,
                    crosstalk_params = c(protdat = 0.8, phosphdat = 0.4, meta = 0.2),
                    degree_bias = list(layers = "prot", gamma = 0.1),
                    restart = 0.5,
                    seed_weights = NULL,
                    output = "list")
  
  lapply(r, head)
}


#' @title Perform Random Walk with Restart
#' 
#' @description 
#' Perform Random Walk with Restart (RWR) with a transition matrix, seed set, restart parameter, and optional seed weights.
#' 
#' @details 
#' 
#' The seed values give the probability of the random walker to start at a node and represents the prior importance of each node. Seed weights give the relative importance between seed sets from different layers in a multilayer network.
#' 
#' By default, `RWR()` implements a node-specific restart parameter based on the user-defined _restart_ parameter. For each node _restart_ is increased or decreased as an inverse function of that nodes coreness. Therefore, if a node is located in a dense region of the network (high coreness), then its _restart_ value will decrease to enable the random walker to better explore the neighborhood around this node.
#' 
#' @param tmat transition matrix inheriting class "dgCMatrix" from the "Matrix" package
#' @param seeds 
#' @param seed_weights
#' @param restart
#' @param network_hierarchy
#' @param node_specific_restart
#' @param node_connectivity_scores
#' @param max_iters
#' 
#' @return named numeric vector
#' 
#' @export
#' 
RWR <- function(tmat, seeds = NULL, seed_weights, restart = 0.5, network_hierarchy, node_specific_restart = TRUE, node_connectivity_scores, max_iters = 100) {
  if(0){
    source("/Users/samboyd/Documents/NHNM/R package/TEST/TEST_create_integrated_network.R")
    source("/Users/samboyd/Documents/NHNM/R package/NHNM/R/create_integrated_network.R")
    source("/Users/samboyd/Documents/NHNM/R package/NHNM/R/RandomWalk.R")
    source("/Users/samboyd/Documents/NHNM/R package/NHNM/R/AMEND.R")
    
    library(igraph)
    library(Matrix)
    
    s.t = Sys.time()
    g = create_integrated_network(network_layers = network_layers, 
                                  bipartite_networks = bipartite_networks,
                                  network_hierarchy = network_hierarchy,
                                  data = data,
                                  FUN = FUN,
                                  FUN_params = FUN_params,
                                  directed = directed,
                                  brw_attr = brw_attr,
                                  lcc = lcc)

    tmat = transition_matrix(network = g$network,
                             network_hierarchy = g$hierarchy,
                             norm = "modified_degree",
                             k = 0.5,
                             crosstalk_params = c(protdat = 0.8, phosphdat = 0.4, meta = 0.2),
                             degree_bias = list(layers = V(g$hierarchy)$name[V(g$hierarchy)$level > 0], gamma = 0.1))


    seeds = NULL
    restart = 0.05
    seed_weights = NULL
    network_hierarchy = g$hierarchy
    max_iters = 100
    node_specific_restart <- TRUE
    
    node_connectivity_scores = node_connectivity_score(graph = g$network, net_diam_prop = -1, inverse =TRUE)
    
  }
  
  #=== Function settings ===#
  DEFAULT_RESTART <- 0.5
  convergence_tol <- 1e-08
  zero_tol <- 1e-10
  max_iter_radj <- 400
  
  #=== restart checks ===#
  if(is.null(restart) || is.na(restart) || restart < 0 || restart > 100) {
    restart <- DEFAULT_RESTART
  } else if(restart > 1 && restart < 100) {
    restart <- restart / 100
  } 
  
  #=== network_hierarchy checks ===#
  if(!inherits(network_hierarchy, "hierarchy")) stop("network_hierarchy must be an object of class 'hierarchy'.")
  
  #=== seed_weights checks ===#
  # list of named numeric vectors or NULL
  if(is.null(seed_weights)) {
    sw <- seed_weight_sets(network_hierarchy)
    names(sw) <- unlist(lapply(sw, function(x) paste(x, collapse = "|")))
    
    seed_weights <- vector("list", length(sw))
    for(i in seq_along(seed_weights)) {
      seed_weights[[i]] <- rep(1/length(sw[[i]]), length(sw[[i]]))
      names(seed_weights[[i]]) <- sw[[i]]
    }
    names(seed_weights) <- names(sw)
  } else if(is.list(seed_weights)) {
    sw <- seed_weight_sets(network_hierarchy)
    names(sw) <- unlist(lapply(sw, function(x) paste(x, collapse = "|")))
    
    for(i in seq_along(seed_weights)) {
      if(any(seed_weights[[i]] < 0) || all(seed_weights[[i]] == 0)) stop("Numeric vectors in seed_weights must be non-negative with at least one non-zero value.")
      if(sum(seed_weights[[i]]) != 1){
        seed_weights[[i]] = seed_weights[[i]] / sum(seed_weights[[i]])
      } 
    }
    names(seed_weights) <- unlist(lapply(seed_weights, function(x) paste(names(x), collapse = "|")))
    
    all_sw <- unlist(lapply(sw, function(x) paste(x, collapse = "|")))
    missing_sw <- setdiff(all_sw, names(seed_weights))
    if(length(missing_sw) > 0) {
      for(i in seq_along(missing_sw)) {
        seed_weights[[length(seed_weights) + 1]] <- rep(1/length(sw[[missing_sw[i]]]), length(sw[[missing_sw[i]]]))
        names(seed_weights[[length(seed_weights)]]) <- extract_string(missing_sw[i], "\\|", 0)
      }
      names(seed_weights)[(length(seed_weights) - length(missing_sw) + 1):length(seed_weights)] <- missing_sw
    }
  } else stop("Unrecognized input for seed_weights. Must be a list of named numeric vectors, or NULL.")
  
  
  #=== seeds check ===#
  if(is.null(seeds)) {
    seeds <- matrix(1, ncol = 1, nrow = nrow(tmat), dimnames = list(rownames(tmat)))
  } else if(is.data.frame(seeds)) {
    seeds <- as.matrix(seeds)[,1, drop=FALSE]
  } else if(is.numeric(seeds)) {
    seeds <- as.matrix(seeds, ncol = 1)
  } else if(is.matrix(seeds)) {
    seeds <- seeds[,1, drop=FALSE]
  } else stop("Unrecognized input for 'seeds'. Must be a matrix, data.frame, numeric vector, or NULL.")
  
  if(is.null(rownames(seeds))) {
    stop("The function requires the rownames (if matrix/data.frame) or names (if numeric) of 'seeds'.")
  } 
  if(any(is.na(rownames(data)))) {
    warning("seeds with NA as rownames will be removed")
    seeds <- seeds[!is.na(rownames(seeds)), , drop=FALSE]
  }
  # Add missing seeds
  if(any(!rownames(tmat) %in% rownames(seeds))) {
    new_names <- setdiff(rownames(tmat), rownames(seeds))
    zero_seeds <- matrix(0, ncol = 1, nrow = length(new_names), dimnames = list(new_names))
    seeds <- rbind(seeds, zero_seeds)
  }
  # Remove seeds not in tmat
  seeds <- seeds[rownames(seeds) %in% rownames(tmat),, drop=FALSE]
  # Get seeds in same order as tmat
  id <- match(rownames(tmat), rownames(seeds))
  seeds <- seeds[id,, drop=FALSE]
  
  
  #================#
  # Normalize seeds
  #================#
  seeds <- normalize_seeds(seeds = seeds, layers = extract_string(rownames(seeds), "\\|", 2), network_hierarchy = network_hierarchy, seed_weights = seed_weights)
  
  
  #====#
  # RWR
  #====#
  # By Iterative Matrix Multiplication
  p_i <- seeds
  if(restart != 1) {
    if(node_specific_restart) {
      d0 <- node_connectivity_scores
      l <- (sum(d0 * p_i)) / (sum(p_i ^ 2))
      d <- d0 - l * p_i
      for(i in seq_len(max_iter_radj)) {
        if(all(d >= -restart) && all(d <= 1-restart)) break
        d[d < -restart] <- -restart
        d[d > 1-restart] <- 1 - restart
        l <- (sum(d * p_i)) / (sum(p_i ^ 2))
        d <- d - l * p_i
      }
    } else {
      d <- 0
    }
    
    for(i in seq_len(max_iters)) {
      p_new <- (1 - restart) * tmat %*% p_i + (restart + d) * seeds
      if(sum(abs(p_new - p_i)) < convergence_tol) break  # Convergence check
      p_i <- p_new
    }
    p_i[p_i < zero_tol] <- 0
    p_i <- sum2one(p_i, mode = "dense")
  }
  p_i <- p_i[,1]
  
  if(0){ # diagnostics
    plot(p_i0, p_i1); abline(a=0,b=1)
    delta = ecdf(p_i1)(p_i1) - ecdf(p_i0)(p_i0)
    cn = coreness(g$network)
    dg = degree(g$network)
    
    m = c("pearson", "spearman", "kendall")[2]
    
    par(mfrow=c(1,1))
    plot(cn, delta); abline(h=0)
    hist(delta[cn < 20])
    cor(delta, cn, method = m)
    cor(delta, dg, method = m)
    
    par(mfrow=c(1,1))
    plot(d0, d); abline(a=0,b=1)
    summary(d)
    
    cor(dg, d, method = m)
    cor(cn, d, method = m)
    
    cor(dg, p_i0, method = m)
    cor(dg, p_i1, method = m)
    cor(cn, p_i0, method = m)
    cor(cn, p_i1, method = m)
    
    par(mfrow=c(1,2))
    plot(dg, p_i0, ylim=c(0, 0.0056))
    plot(dg, p_i1, ylim=c(0, 0.0056))
    
    plot(dg, p_i0, ylim=c(0, 0.001), xlim = c(0,200))
    plot(dg, p_i1, ylim=c(0, 0.001), xlim = c(0,200))
    
    par(mfrow=c(1,1))
    plot(rank(p_i0), rank(p_i1)); abline(a=0,b=1, col='red')
  }
  
  return(p_i)
}

if(0){ # TEST
  source("/Users/samboyd/Documents/NHNM/R package/TEST/TEST_create_integrated_network.R")
  source("/Users/samboyd/Documents/NHNM/R package/NHNM/R/create_integrated_network.R")
  source("/Users/samboyd/Documents/NHNM/R package/NHNM/R/RandomWalk.R")
  
  library(igraph)
  library(Matrix)
  
  s.t = Sys.time()
  g = create_integrated_network(network_layers = network_layers, 
                                bipartite_networks = bipartite_networks,
                                network_hierarchy = network_hierarchy,
                                data = data,
                                FUN = FUN,
                                FUN_params = FUN_params,
                                directed = directed,
                                brw_attr = brw_attr,
                                lcc = lcc)
  Sys.time() - s.t
  
  s.t = Sys.time()
  tmat = transition_matrix(network = g$network,
                           network_hierarchy = g$hierarchy,
                           norm = "modified_degree",
                           k = 0.5,
                           crosstalk_params = c(protdat = 0.8, phosphdat = 0.4, meta = 0.2),
                           degree_bias = list(layers = V(g$hierarchy)$name[V(g$hierarchy)$level > 0], gamma = 0.1))
  Sys.time() - s.t

  s.t = Sys.time()
  r = RWR(tmat = tmat, seeds = NULL, restart = 0.05, seed_weights = NULL, network_hierarchy = g$hierarchy)
  Sys.time() - s.t

}


#' @title Create a transition matrix
#' 
#' @description Create a transition matrix
#' 
#' @param t t
#' 
#' @return t
#' 
#' @export
#' 
transition_matrix <- function(network, network_hierarchy, norm = c("degree", "modified_degree"), k = 0.5,
                              crosstalk_params = NULL, degree_bias = NULL, in_parallel = FALSE, n_cores = NULL) {
  if(0){ # TEST
    network = g$network
    network_hierarchy = g$hierarchy
    norm = "modified_degree"
    k = 0.5
    crosstalk_params = c(protdat = 0.8, phosphdat = 0.4, meta = 0.2)
    degree_bias = list(layers = V(g$hierarchy)$name[V(g$hierarchy)$level > 0], gamma = 0.25)
    in_parallel = FALSE
    n_cores = NULL
    #===========#
    network = orig_net
    network_hierarchy = network_hierarchy
    norm = norm
    k = k
    crosstalk_params = crosstalk_params
    degree_bias = degree_bias
    in_parallel = in_parallel
    n_cores = n_cores
  }
  
  #=== Function Settings ===#
  DEFAULT_CROSSTALK <- 0.5
  DEFAULT_GAMMA <- 0.2
  
  norm <- match.arg(norm)
  
  #=== network_hierarchy checks ===#
  if(!inherits(network_hierarchy, "hierarchy")) stop("network_hierarchy must be an object of class 'hierarchy'.")
  
  #=== network checks ===#
  if(!inherits(network, "NHNMgraph")) stop("network must be an object of class 'NHNMgraph'.")
  
  # Convert to adjacency matrix with class dgCMatrix from Matrix package
  adj_mat <- igraph::as_adjacency_matrix(graph = network, attr = "weight", sparse = TRUE)
  
  #=== crosstalk_params Checks ===#
  # crosstalk_params is a named numeric vector, or NULL. If non-NULL, names correspond to nodes in hierarchy.
  # Values must be between zero and one.
  if(is.null(crosstalk_params)) {
    crosstalk_params <- rep(DEFAULT_CROSSTALK, igraph::vcount(network_hierarchy) - 1)
    names(crosstalk_params) <- V(network_hierarchy)$name[V(network_hierarchy)$level > 0]
  }else if(is.numeric(crosstalk_params)) {
    if(is.null(names(crosstalk_params))) stop("crosstalk_params must be a named numeric vector or NULL.")
    if(any(crosstalk_params > 1) || any(crosstalk_params < 0)) stop("crosstalk_params values must be between 0 and 1.")
    
    add_ct_params <- rep(DEFAULT_CROSSTALK, sum(!V(network_hierarchy)$name %in% names(crosstalk_params)) - 1)
    names(add_ct_params) <- V(network_hierarchy)$name[!V(network_hierarchy)$name %in% names(crosstalk_params) & V(network_hierarchy)$level != 0]
    
    crosstalk_params <- c(crosstalk_params, add_ct_params)
  }else stop("Unrecognized input for crosstalk_params. Must be a named numeric vector or NULL.")
  
  #=== degree_bias Checks ===#
  # Must be NULL, a character vector, or a named list
  if(!is.null(degree_bias)) {
    if(is.character(degree_bias)) {
      if(any(!degree_bias %in% V(network_hierarchy)$name)) {
        degree_bias <- degree_bias[degree_bias %in% V(network_hierarchy)$name]
      } 
      degree_bias <- list(layers = degree_bias, gamma = DEFAULT_GAMMA)
    } 
    if(is.list(degree_bias)) {
      if(is.null(names(degree_bias))) stop("List is not named. If a list, degree_bias must have names 'layers' and 'gamma'.")
      if(any(!c("layers", "gamma") %in% names(degree_bias))) stop("If a list, degree_bias must have names 'layers' and 'gamma'.")
      if(any(!degree_bias$layers %in% V(network_hierarchy)$name)) {
        degree_bias$layers <- degree_bias$layers[degree_bias$layers %in% V(network_hierarchy)$name]
      } 
      if(!is.numeric(degree_bias$gamma) || length(gamma) != 1) stop("degree_bias$gamma must be a numeric scalar.")
    }else stop("Unrecognized input for degree_bias. Must be a character vector, a named list, or NULL.")
  }
  
  #===============================#
  # Normalize the adjacency matrix
  #===============================#
  # Layers of each node in adjacency matrix
  all_layers <- extract_string(rownames(adj_mat), "\\|", 2)
  
  # BOTTLENECK
  tmat <- normalize_adjmat(adj_mat = adj_mat,
                           norm = norm,
                           k = k,
                           layers = all_layers,
                           brw_attr = V(network)$brw_attr,
                           network_hierarchy = network_hierarchy,
                           crosstalk_params = crosstalk_params,
                           degree_bias = degree_bias,
                           in_parallel = in_parallel,
                           n_cores = n_cores)
  
  return(tmat)
}

# TEST
if(0){
  source("/Users/samboyd/Documents/NHNM/R package/TEST/TEST_create_integrated_network.R")
  source("/Users/samboyd/Documents/NHNM/R package/NHNM/R/create_integrated_network.R")
  source("/Users/samboyd/Documents/NHNM/R package/NHNM/R/RandomWalk.R")
  
  g = create_integrated_network(network_layers = network_layers, 
                                bipartite_networks = bipartite_networks,
                                network_hierarchy = network_hierarchy,
                                data = data,
                                FUN = FUN,
                                FUN_params = FUN_params,
                                directed = directed,
                                brw_attr = brw_attr,
                                lcc = lcc)
  net = g
  
  vcount(g$network)
  s.t = Sys.time()
  tmat = transition_matrix(network = g$network,
                           network_hierarchy = g$hierarchy,
                           norm = "degree",
                           crosstalk_params = c(protdat = 0.8, phosphdat = 0.4, meta = 0.2),
                           degree_bias = list(layers = V(g$hierarchy)$name[V(g$hierarchy)$level > 0], gamma = 0.25))
  Sys.time() - s.t
  
  
  # Investigate
  if(0){
    tmat[1:5, 1:5]
    nnzero(tmat)
    min(tmat@x)
    summary(tmat@x)
    d = drop0(tmat1 - tmat)
    summary(d@x)
    median(d@x)
    # Intra-layer
    leaves = V(net$hierarchy)$name[V(net$hierarchy)$level == 3]
    median_diff = numeric(length(leaves)); names(median_diff) = leaves
    for(i in seq_along(median_diff)) {
      ids <- which(V(net$network)$layer %in% leaves[i])
      tmp_mat <- tmat[ids, ids]
      tmp_mat1 <- tmat1[ids, ids]
      d = drop0(tmp_mat1 - tmp_mat)
      median_diff[i] = median(d@x) 
    }
    median_diff
    # Inter-layer k
    k = 1
    sibs <- get_siblings(net$hierarchy, level = k)
    median_diff = NULL
    for(i in seq_along(sibs)) {
      if(length(sibs[[i]]) == 1) next
      leaf_nodes <- get_leaf_nodes(net$hierarchy, node = sibs[[i]])
      for(j in 1:(length(sibs[[i]])-1)) {
        r_ids <- which(V(net$network)$layer %in% leaf_nodes[j])
        for(l in (j+1):length(sibs[[i]])) {
          c_ids <- which(V(net$network)$layer %in% leaf_nodes[l])
          tmp_mat <- tmat[r_ids, c_ids]
          tmp_mat1 <- tmat1[r_ids, c_ids]
          d = drop0(tmp_mat1 - tmp_mat)
          median_diff = c(median_diff, median(d@x))
        }
      }
    }
    median_diff
    # The new transition matrix indeed has... 
    # 1) larger intra-layer transition probabilities,
    # 2) larger inter-category transition probabilities for categories closer together in hierarchy,
    # 3) smaller inter-category transition probabilities for categories further apart in hierarchy.
    
  }
  
  
  
}


#' @title Extract the diagonal elements of a matrix
#'
#' @description Extract the diagonal elements of a square matrix.
#'
#' @param X dgCMatrix sparse matrix
#'
#' @returns vector of diagonal elements of X.
#'
#' @export
#' 
get_diagonal <- function(X) {
  d <- numeric(ncol(X))
  for(i in seq_len(ncol(X))) {
    id_tmp <- (X@p[i] + 1):X@p[i+1] # ids of X@x that are non-zero and in col i
    row_ids <- X@i[id_tmp] + 1 # row ids of non-zero elements in col i
    if(i %in% row_ids) { # if diagonal is non-zero
      d[i] <- X@x[id_tmp[row_ids == i]]
    }
  }
  return(d)
}


#' @title Iterative Proportional Fitting
#'
#' @description Scales input matrix X to have row and column sums approximately equal to 1, thereby transforming X to a bistochastic (doubly stochastic) matrix.
#'
#' @param X Matrix
#' @param gamma
#' @param diag_weight Edge weight to add to zero-value diagonal elements of X to aid on convergence
#'
#' @returns A named list
#'  B: scaled matrix from the matrix product PXQ
#'  p: diagonal elements of diagonal matrix P
#'  q: diagonal elements of diagonal matrix Q
#'
#'  @noRd
#'  
ipf <- function(X, gamma = 1, diag_weight = 1e-8) {
  # For starting transition matrix X, obtain B = PXQ, where B is bistochastic, P,Q are diagonal matrices, and B and X are as similar as possible (minimize relative entropy between B & X)
  # Set target row sums depending on how much to homogenize degree influence. gamma = 1 for complete homogenization, gamma = 0 for no correction.
  
  #=== Function Settings ===#
  stop_delta <- 1e-10
  stop_step <- 200
  #=========================#
  
  if(gamma < 0 || gamma > 1) stop("gamma must be between 0 and 1.")
  
  get_rowsum_targets <- function(rowsum, gamma, target) {
    rowsum - gamma * (rowsum - target)
  } 
  target_row_sums <- get_rowsum_targets(Matrix::rowSums(X), gamma, 1)
  target_col_sums <- 1
  
  x_dimnames <- dimnames(X)
  
  # Adding self-loops to aid in convergence
  d <- get_diagonal(X)
  d[d == 0] <- diag_weight
  X <- X + Matrix::Diagonal(n = nrow(X), x = d)
  
  q1 <- rep(1, nrow(X)) # initialize the diagonal elements of Q
  p1 <- rep(1, nrow(X)) # initialize the diagonal elements of P
  for(i in seq_len(stop_step)) {
    # set p, given q (p are the diagonal elements of P)
    p2 <- target_row_sums / Matrix::rowSums(X %*% Matrix::.sparseDiagonal(x = q1))
    # set q, given p found above (q are the diagonal elements of Q)
    q2 <- target_col_sums / Matrix::colSums(Matrix::.sparseDiagonal(x = p2) %*% X)
    
    delta_p <- sum(abs(p2 - p1)) <= stop_delta
    delta_q <- sum(abs(q2 - q1)) <= stop_delta
    if(delta_p && delta_q) break
    q1 <- q2
    p1 <- p2
  }
  P <- Matrix::.sparseDiagonal(x = p2)
  Q <- Matrix::.sparseDiagonal(x = q2)
  B <- P %*% X %*% Q
  dimnames(B) <- x_dimnames
  return(B)
}


#' @title Bistochastic Scaling
#'
#' @description 
#' Directly modify the transition matrix to attenuate the influence of degree on diffusion scores (as evidenced by an increased entropy of stationary distribution associated with the modified transition matrix).
#' This is done by scaling the transition matrix to be approximately bistochastic (all row & column sums equal 1).
#'
#' @param tmat transition matrix
#' @param gamma 
#'
#' @returns A modified transition matrix that is approximately bistochastic
#'
#' @export
#' 
bistochastic_scaling <- function(tmat, gamma = 1) {
  if(gamma == 0) return(tmat)
  B <- ipf(X = tmat, gamma = gamma)
  b <- get_diagonal(B)
  tmp_res <- numeric(Matrix::nnzero(B))
  for(i in seq_len(ncol(B))){
    if(b[i] == 0) next
    id_tmp <- (B@p[i] + 1):B@p[i+1] # ids of B@x that are non-zero and in col i
    row_ids <- B@i[id_tmp] + 1 # row ids of non-zero elements in col i
    diag_id <- id_tmp[row_ids == i]
    off_id <- id_tmp[row_ids != i]
    # Evenly distribute to neighbors of node i. Preserves column sums
    tmp_res[diag_id] <- 0
    tmp_res[off_id] <- B@x[off_id] + b[i] / (B@p[i+1] - B@p[i] - 1) # Denominator: number of non-zero elements in col i i.e., degree of node i. Minus 1 b/c of self-loops
  }
  B@x <- tmp_res
  return(Matrix::drop0(B))
}

# TEST
if(0){
  res = bistochastic_scaling(tmat, 1)
  res[1:5,1:5]
  mean(abs(res[res != 0] - tmat[tmat != 0]))
  
  res = bistochastic_scaling(tmat, 0.9)
  res[1:5,1:5]
  mean(abs(res[res != 0] - tmat[tmat != 0]))
  
  res = bistochastic_scaling(tmat, 0.8)
  res[1:5,1:5]
  mean(abs(res[res != 0] - tmat[tmat != 0]))
  
  res = bistochastic_scaling(tmat, 0.25)
  res[1:5,1:5]
  mean(abs(res[res != 0] - tmat[tmat != 0]))
  
  res = bistochastic_scaling(tmat, 0.1)
  res[1:5,1:5]
  mean(abs(res[res != 0] - tmat[tmat != 0]))
  
  tmat[1:5,1:5]
  
}


#' @title Scale columns of a matrix to sum to one
#'
#' @description Scale a matrix such that all columns sum to one. Relies on the 'Matrix' package.
#'
#' @param X A sparse matrix of class "dgCMatrix" from "Matrix" package
#' @param mode 
#' @param tol tolerance for numerical precision
#'
#' @returns A matrix, with all columns summing to one.
#'
#' @export
#'
sum2one <- function(X, mode = "sparse", tol = 1e-12) {
  if(mode == "sparse") {
    col_sums <- Matrix::colSums(X)
    col_sums[col_sums < tol] <- 1
    col_sums_inv <- 1 / col_sums
    X@x <- X@x * rep(col_sums_inv, diff(X@p))
    X@x[abs(X@x) < tol] <- 0
    return(Matrix::drop0(X))
  } else if(mode == "dense") {
    if(ncol(X) == 1) {
      X <- X / colSums(X)
    } else X <- X / rep(colSums(X), each = nrow(X))
    return(X)
  }else stop("Unrecognized string for 'mode'.")
}

# TEST
if(0){
  M=100000000
  X=matrix(rexp(M),ncol=1)
  time=matrix(nrow=10,ncol=2)
  for(i in seq_len(nrow(time))){
    s.t = Sys.time()
    x = sum2one(X,mode="dense")
    time[i,1] = Sys.time() - s.t
    s.t = Sys.time()
    x = sum2one(X[,1],mode="vector")
    time[i,2] = Sys.time() - s.t
  }
  apply(time, 2, summary)
}


#' @title Required seed_weight parameters
#' 
#' @description Get sets of siblings in hierarchy for which seed weight parameters are necessary. These will be any sibling sets with two or more nodes.
#' 
#' @param hierarchy igraph of network hierarchy
#' 
#' @return list. list element will be an empty list if no seed weight parameters are required.
#' 
#' @export
#' 
seed_weight_sets <- function(hierarchy) {
  sibs <- get_siblings(hierarchy)
  res <- lapply(sibs, function(x) x[sapply(x, length) > 1] )
  res <- unlist(res, recursive = FALSE)
  return(res)
}


#' @title Normalize an adjacency matrix in a hierarchical fashion
#'
#' @description `normalize_adjmat()` efficiently constructs a transition matrix from an adjacency matrix whose nodes are categorized by a hierarchy. 
#'
#' @details This is a recursive function.
#'
#' @param adj_mat an adjacency matrix. Should of class dgCMatrix.
#' @param layers vector of all layer names in same order as rows of adj_mat.
#' @param network_hierarchy 
#' @param crosstalk_params
#' @param level
#' @param degree_bias
#' @param in_parallel
#' @param n_cores
#'
#' @returns A column-normalized matrix containing pairwise transition probabilities
#' 
#' @noRd
#'
normalize_adjmat <- function(adj_mat, norm, k, layers, brw_attr, network_hierarchy, crosstalk_params, level = 1, degree_bias, ...) {
  if(0){ # TEST
    layers = all_layers
    brw_attr = V(network)$brw_attr
    #==============#
    adj_mat = adj_mat
    norm = norm
    k = k
    layers = all_layers
    brw_attr = V(network)$brw_attr
    network_hierarchy = network_hierarchy
    crosstalk_params = crosstalk_params
    degree_bias = degree_bias
    in_parallel = in_parallel
    n_cores = n_cores
  }
  
  if(level == max(V(network_hierarchy)$level)) { 
    #=== BASE CASE ===#
    # For each sibling set at bottom level:
    #   normalize intra-layer adj mat
    #   normalize inter-layer adj mat
    #   apply crosstalk parameters to sibling sets
    
    sibs <- get_siblings(network_hierarchy, level)
    
    for(s in seq_along(sibs)) {
      # Normalize intra-layer and inter-layer adj mat
      for(i in seq_along(sibs[[s]])) { # Row indices
        r_ids <- which(layers == sibs[[s]][i])
        for(j in seq_along(sibs[[s]])) { # Column indices
          c_ids <- which(layers == sibs[[s]][j])
          if(norm == "degree") {
            adj_mat[r_ids, c_ids] <- sum2one(scale_rows(adj_mat[r_ids, c_ids, drop=FALSE], x = brw_attr[r_ids]))
          } else if(norm == "modified_degree") {
            adj_mat[r_ids, c_ids] <- sum2one(scale_rows(adjM = adj_mat[r_ids, c_ids, drop=FALSE], x = brw_attr[r_ids], k = k))
          } else stop("Unrecognized string for norm argument. Must be one of 'degree' or 'modified_degree'.")
          # Apply degree bias adjustment method
          if(!is.null(degree_bias)) {
            if(i == j) {
              if(sibs[[s]][i] %in% get_leaf_nodes(network_hierarchy, node = degree_bias$layers)) {
                adj_mat[r_ids, c_ids] <- bistochastic_scaling(adj_mat[r_ids, c_ids, drop=FALSE], gamma = degree_bias$gamma)
              }
            }
          }
        }
      }
      
      if(length(sibs[[s]]) > 1) {
        # Apply crosstalk parameters
        ids <- which(layers %in% sibs[[s]])
        tmp_mat <- adj_mat[ids, ids, drop=FALSE]
        tmp_layers <- extract_string(rownames(tmp_mat), "\\|", 2)
        parents <- get_parents(g = network_hierarchy, h_nodes = tmp_layers)
        loop_ids <- which(Matrix::colSums(tmp_mat) != 0)
        for(i in loop_ids) {
          x_ids <- (tmp_mat@p[i] + 1):tmp_mat@p[i+1] # IDs of @x for column i
          current <- tmp_layers[i] # layer of node associated with column i
          others <- tmp_layers[tmp_mat@i[x_ids] + 1] # Other layers that this node is connected to
          if(all(others != current)) {
            tmp_mat@x[x_ids] <- tmp_mat@x[x_ids] * (1 / length(unique(others)))
            next
          }
          jump_prob <- (1 - crosstalk_params[parents[i]]) * prod(crosstalk_params[tmp_layers[i]])
          n_interlayer_links <- sum(unique(others) != current)
          if(n_interlayer_links == 0) {
            next
          } 
          ctp <- ifelse(others == current, 1 - crosstalk_params[current], jump_prob / n_interlayer_links)
          tmp_mat@x[x_ids] <- tmp_mat@x[x_ids] * ctp
        }
        
        adj_mat[ids, ids] <- sum2one(tmp_mat)
      }
    }
    
    return(adj_mat)
  }
  
  #=== RECURSIVE CASE ===#
  adj_mat <- normalize_adjmat(adj_mat = adj_mat, norm = norm, k = k, layers = layers, brw_attr = brw_attr, network_hierarchy = network_hierarchy, crosstalk_params = crosstalk_params, level = level + 1, degree_bias = degree_bias, in_parallel = in_parallel, n_cores = n_cores)
  
  # For each sibling set at current level:
  #   normalize inter-unit sections of norm_adj_mat
  #   apply crosstalk parameters to sibling sets in norm_adj_mat
  
  sibs <- get_siblings(network_hierarchy, level)
  
  for(s in seq_along(sibs)) {
    # Normalize inter-layer adj mat
    for(i in seq_along(sibs[[s]])) { # Row indices
      r_ids <- which(layers %in% get_leaf_nodes(network_hierarchy, node = sibs[[s]][i]))
      for(j in seq_along(sibs[[s]])) { # Column indices
        if(i == j) next
        c_ids <- which(layers %in% get_leaf_nodes(network_hierarchy, node = sibs[[s]][j]))
        if(norm == "degree") {
          adj_mat[r_ids, c_ids] <- sum2one(scale_rows(adj_mat[r_ids, c_ids, drop=FALSE], x = brw_attr[r_ids]))
        } else if(norm == "modified_degree") {
          adj_mat[r_ids, c_ids] <- sum2one(scale_rows(adjM = adj_mat[r_ids, c_ids, drop=FALSE], x = brw_attr[r_ids], k = k))
        } else stop("Unrecognized string for norm argument. Must be one of 'degree' or 'modified_degree'.")
      }
    }
    
    if(length(sibs[[s]]) > 1) {
      # Apply crosstalk parameters
      leaf_nodes <- get_leaf_nodes(network_hierarchy, node = sibs[[s]])
      ids <- which(layers %in% leaf_nodes)
      tmp_mat <- adj_mat[ids, ids, drop=FALSE]
      tmp_layers <- extract_string(rownames(tmp_mat), "\\|", 2)
      tmp_cats <- names(leaf_nodes)[match(tmp_layers, leaf_nodes)]
      parents <- get_parents(g = network_hierarchy, h_nodes = tmp_cats)
      offspring <- get_descendants(g = network_hierarchy, h_nodes = tmp_cats, leaves = tmp_layers)
      loop_ids <- which(Matrix::colSums(tmp_mat) != 0)
      for(i in loop_ids) {
        x_ids <- (tmp_mat@p[i] + 1):tmp_mat@p[i+1] # IDs of @x for column i
        current <- tmp_cats[i] # layer of node associated with column i
        others <- tmp_cats[tmp_mat@i[x_ids] + 1] # Other layers that this node is connected to
        if(all(others != current)) {
          tmp_mat@x[x_ids] <- tmp_mat@x[x_ids] * (1 / length(unique(others)))
          next
        }
        ctp_tmp <- ifelse(level == 1, 0, crosstalk_params[parents[i]])
        jump_prob <- (1 - ctp_tmp) * prod(crosstalk_params[offspring[[i]]])
        n_interlayer_links <- sum(unique(others) != current)
        if(n_interlayer_links == 0) { # equivalent to all(others == current)
          next
        } 
        ctp <- ifelse(others == current, 1, jump_prob / n_interlayer_links)
        tmp_mat@x[x_ids] <- tmp_mat@x[x_ids] * ctp
      }
      
      adj_mat[ids, ids] <- sum2one(tmp_mat)
    }
  }
  
  return(adj_mat)
}


#' @title Normalize a seed vector in a hierarchical fashion
#'
#' @description `normalize_seeds()` efficiently normalizes a seed vector or 1-column matrix, where nodes are categorized by a hierarchy. 
#'
#' @details This is a recursive function.
#'
#' @param seeds an adjacency matrix. Should of class dgCMatrix.
#' @param layers vector of all layer names in same order as rows of adj_mat.
#' @param network_hierarchy 
#' @param seed_weights
#' @param level
#'
#' @returns A named numeric vector containing initial probabilities for RWR
#' 
#' @noRd
#'
normalize_seeds <- function(seeds, layers, network_hierarchy, seed_weights, level = 1) {
  if(0){ # TEST
    layers = extract_string(rownames(seeds), "\\|", 2)
    
  }
  
  if(level == max(V(network_hierarchy)$level)) { 
    #=== BASE CASE ===#
    # For each sibling set at bottom level:
    #   normalize seeds
    #   apply seed weight parameters to sibling sets
    
    sibs <- get_siblings(network_hierarchy, level)
    for(s in seq_along(sibs)) {
      sw <- seed_weights[[paste(sibs[[s]], collapse = "|")]]
      # Normalize seeds of each layer
      for(i in seq_along(sibs[[s]])) { # Row indices
        ids <- which(layers == sibs[[s]][i])
        seeds[ids,] <- sum2one(seeds[ids,, drop=FALSE], mode = "dense")
        if(length(sibs[[s]]) > 1) {
          # Apply seed_weights
          seeds[ids,] <- sw[sibs[[s]][i]] * seeds[ids,]
        }
      }
    }
    return(seeds)
  }
  
  #=== RECURSIVE CASE ===#
  seeds <- normalize_seeds(seeds = seeds, layers = layers, network_hierarchy = network_hierarchy, seed_weights = seed_weights, level = level + 1)
  
  # For each sibling set at current level:
  #   apply seed weight parameters to sibling sets
  
  sibs <- get_siblings(network_hierarchy, level)
  for(s in seq_along(sibs)) {
    if(length(sibs[[s]]) > 1) {
      # Apply seed_weights
      sw <- seed_weights[[paste(sibs[[s]], collapse = "|")]]
      for(i in seq_along(sibs[[s]])) { 
        ids <- which(layers %in% get_leaf_nodes(network_hierarchy, node = sibs[[s]][i]))
        seeds[ids,] <- sw[sibs[[s]][i]] * seeds[ids,]
      }
    }
  }
  if(level == 1) {
    seeds <- sum2one(seeds, mode = "dense")
  }
  return(seeds)
}


#' @title 
#' 
#' @description 
#' 
#' @param 
#' 
#' @return 
#' 
#' @noRd
#'
get_parents <- function(g, h_nodes) {
  # For given hierarchical nodes, find parents
  uniq_h_nodes <- unique(h_nodes)
  ids <- match(uniq_h_nodes, V(g)$name)
  nei <- unlist(lapply(igraph::neighborhood(graph = g, order = 1, nodes = ids, mode = 'in'), function(x) x[!x %in% ids]))
  h_node_parents <- names(nei)[match(h_nodes, uniq_h_nodes)]
  return(h_node_parents)
}

# TEST
if(0){ 
  table(get_parents(net$hierarchy, V(net$network)$layer))
  table(get_parents(net$hierarchy, get_parents(net$hierarchy, V(net$network)$layer)))
  table(get_parents(net$hierarchy, get_parents(net$hierarchy, get_parents(net$hierarchy, V(net$network)$layer))))
}


#' @title 
#' 
#' @description 
#' 
#' @param 
#' 
#' @return 
#' 
#' @noRd
#'
get_descendants <- function(g, h_nodes, leaves) {
  comb <- paste(h_nodes, leaves, sep = "|")
  uniq_comb <- unique(comb)
  tmp <- vector("list", length(uniq_comb)); names(tmp) <- uniq_comb
  for(i in seq_along(tmp)) {
    l <- extract_string(names(tmp)[i], "\\|", 0)
    node_id <- match(l[1], V(g)$name)
    leaf_id <- match(l[2], V(g)$name)
    tmp[[i]] <- as.numeric(igraph::shortest_paths(graph = g, from = node_id, to = leaf_id, mode = 'out')$vpath[[1]])
    tmp[[i]] <- V(g)$name[tmp[[i]]]
  }
  res <- unname(tmp[match(comb, uniq_comb)])
  return(res)
}

# TEST
if(0){
  leaves = V(net$network)$layer
  head(get_descendants(g = net$hierarchy, h_nodes = leaves, leaves = leaves))
  layers = get_parents(net$hierarchy, leaves)
  head(get_descendants(g = net$hierarchy, h_nodes = layers, leaves = leaves))
  layers = get_parents(net$hierarchy, layers)
  head(get_descendants(g = net$hierarchy, h_nodes = layers, leaves = leaves))
}


#' @title Get all sets of siblings from network hierarchy
#' 
#' @description Find all sets of child nodes that have a common parent in network hierarchy, i.e., siblings
#' 
#' @param g igraph of network hierarchy
#' @param level An integer giving the level to query, or "all" (default) to query all levels
#' 
#' @return list of sibling sets at each level, given by node name
#' 
#' @noRd
#'
get_siblings <- function(g, level = "all") {
  if(level == 'all') {
    sibs <- vector('list', max(V(g)$level))
    for(j in seq_along(sibs)) {
      # Step 1: For nodes of a given level, find their parent nodes 
      ids <- which(V(g)$level == j)
      nei <- unlist(lapply(igraph::neighborhood(graph = g, order = 1, nodes = ids, mode = 'in'), function(x) x[!x %in% ids]))
      
      # Step 2: Create sets of child nodes that have a common parent, i.e., siblings 
      uniq_parents <- unique(nei)
      sibs_tmp <- vector('list', length(uniq_parents))
      for(i in seq_along(sibs_tmp)) {
        sibs_tmp[[i]] <- V(g)$name[ids[nei == uniq_parents[i]]]
      }
      sibs[[j]] <- sibs_tmp
    }
  } else {
    # Step 1: For nodes of a given level, find their parent nodes 
    ids <- which(V(g)$level == level)
    nei <- unlist(lapply(igraph::neighborhood(graph = g, order = 1, nodes = ids, mode = 'in'), function(x) x[!x %in% ids]))
    
    # Step 2: Create sets of child nodes that have a common parent, i.e., siblings 
    uniq_parents <- unique(nei)
    sibs <- vector('list', length(uniq_parents))
    for(i in seq_along(sibs)) {
      sibs[[i]] <- V(g)$name[ids[nei == uniq_parents[i]]]
    }
  }
  return(sibs)
}


#' @title Scale rows of adjacency matrix
#'
#' @description Scale rows of adjacency matrix by `x`, or by row sums raised to the power `-k` (if `x` is NULL).
#'
#' @param adjM A sparse matrix of class dgCMatrix
#' @param x numeric vector in same order as rows of `adjM`
#' @param k non-negative number. Row sums will be raised to the power `-k` before row scaling.
#'
#' @returns A matrix
#' 
#' @noRd
#'
scale_rows <- function(adjM, x = NULL, k = NULL) {
  if(!is.null(x) && !is.null(k)) {
    row_scaling_values <- x * Matrix::rowSums(adjM)^(-k)
  } else if(!is.null(x) && is.null(k)) {
    row_scaling_values <- x
  } else if(is.null(x) && !is.null(k)){
    row_scaling_values <- Matrix::rowSums(adjM)^(-k)
  } else stop("At least one of 'x' or 'k' must be non-null.")
  # Ensure no NaN or Inf values (to handle division by zero)
  row_scaling_values[!is.finite(row_scaling_values)] <- 0
  # Scale rows
  adjM@x <- adjM@x * row_scaling_values[adjM@i + 1]
  return(Matrix::drop0(adjM))
}



#=============#
# EXPERIMENTAL
#=============#
#' @title Inflation-Normalization Procedure
#'
#' @description Directly modify the transition matrix to attenuate the influence of degree on diffusion scores (as evidenced by an increased entropy of stationary distribution associated with the modified transition matrix).
#' This is done by raising the values in each row of a left-stochastic transition matrix by an exponent greater than one which is a positive linear function of the stationary probability of the node associated with that row. This is followed by column normalization. This procedure displaces incoming transition probabilities to a node to the other outgoing transition probabilities of its neighbors as a function of degree.
#'
#' @param nadjM a column-normalized adjacency matrix (i.e., transition matrix)
#'
#' @returns A modified transition matrix that has been adjusted for degree bias.
#' 
#' @noRd
#'
inflate_normalize <- function(tmat) {
  # Getting stationary distribution of tmat
  stat_distr1 <- stationary_distr(tmat)
  e0 <- entropy(stat_distr1)
  # Performing Inflation-normalization on transition matrix
  kf <- c(1, 10, seq(50, 2000, 50))
  entropy_res <- numeric(length(kf))
  for(j in seq_along(kf)) {
    inflation <- 1 + kf[j] * stat_distr1
    tmat_tmp <- methods::as(tmat, "RsparseMatrix")
    tmp <- numeric(Matrix::nnzero(tmat_tmp))
    for(i in which(Matrix::rowSums(tmat_tmp) != 0)) { # i corresponds to rows of tmat
      tmp[(tmat_tmp@p[i] + 1):tmat_tmp@p[i+1]] = rep(inflation[i], tmat_tmp@p[i+1] - tmat_tmp@p[i])
    }
    tmat_tmp@x <- tmat_tmp@x ^ tmp
    tmat_tmp <- methods::as(tmat_tmp, "CsparseMatrix")
    tmat_tmp <- sum2one(tmat_tmp)
    sd_tmp <- stationary_distr(tmat_tmp)
    if(any(sd_tmp < 0 | sd_tmp > 1)) {
      if(j == 1) {
        return(tmat)
      } else {
        j <- j - 1
        break
      }
    }
    entropy_res[j] <- entropy(sd_tmp)
    if(j == 1) {
      if(entropy_res[j] < e0) {
        return(tmat)
      }
    } else if(entropy_res[j] <= entropy_res[j-1]) {
      break
    }
  }
  entropy_res <- entropy_res[1:j]
  j <- which.max(entropy_res)
  inflation <- 1 + kf[j] * stat_distr1
  tmat <- methods::as(tmat, "RsparseMatrix")
  tmp <- numeric(Matrix::nnzero(tmat))
  for(i in which(Matrix::rowSums(tmat) != 0)) { # i corresponds to rows of tmat
    tmp[(tmat@p[i] + 1):tmat@p[i+1]] = rep(inflation[i], tmat@p[i+1] - tmat@p[i])
  }
  tmat@x <- tmat@x ^ tmp
  tmat <- methods::as(tmat, "CsparseMatrix")
  tmat <- sum2one(tmat)
  return(tmat)
}

# TEST
if(0) {
  tmat_tmp = bistochastic_scaling(tmat)
  nnzero(tmat)
  nnzero(tmat_tmp)
  
  tmat_tmp = inflate_normalize(tmat)
  nnzero(tmat)
  nnzero(tmat_tmp)
}


#' @title Compute entropy of a discrete probability distribution
#'
#' @description
#' Computes entropy (with natural log) of a discrete probability distribution vector (i.e., sums to 1, values between 0 and 1).
#'
#' @param x Probability vector
#'
#' @returns Numeric value representing entropy
#' 
#' @noRd
#'
entropy <- function(x) {
  x <- x[x > 0]
  -sum(x * log(x))
}


#' @title Compute stationary distribution of a transition matrix
#'
#' @description `stationary_distr()` relies on iterative matrix multiplication to compute the stationary distribution of a transition matrix. This is equivalent to the eigenvector associated with the absolute largest eigenvalue of the transition matrix (eigenvector normalized to sum to one). 
#'
#' @param P Left-stochastic transition matrix
#'
#' @returns Probability vector of the stationary distribution.
#' 
#' @noRd
#'
stationary_distr <- function(P, tol = 1e-8, max_iter = 1000) {
  n <- nrow(P)
  p_i <- rep(1 / n, n)  # Initial uniform distribution
  
  for(i in seq_len(max_iter)) {
    p_new <- P %*% p_i  # Matrix-vector multiplication
    if(sum(abs(p_new - p_i)) < tol) break  # Convergence check
    p_i <- p_new
  }
  return(as.numeric(p_i))
}


#=========#
# DETRITUS
#=========#
if(0){
  normalize_adjmat <- function(adj_mat, norm, k, layers, brw_attr, network_hierarchy, crosstalk_params, level = 1, degree_bias, ...) {
    if(0){ # TEST
      layers = all_layers
      brw_attr = V(network)$brw_attr
      #==============#
      adj_mat = adj_mat
      norm = norm
      k = k
      layers = all_layers
      brw_attr = V(network)$brw_attr
      network_hierarchy = network_hierarchy
      crosstalk_params = crosstalk_params
      degree_bias = degree_bias
      in_parallel = in_parallel
      n_cores = n_cores
    }
    
    if(level == max(V(network_hierarchy)$level)) { 
      #=== BASE CASE ===#
      # For each sibling set at bottom level:
      #   normalize intra-layer adj mat
      #   normalize inter-layer adj mat
      #   apply crosstalk parameters to sibling sets
      
      sibs <- get_siblings(network_hierarchy, level)
      
      for(s in seq_along(sibs)) {
        # Normalize intra-layer and inter-layer adj mat
        for(i in seq_along(sibs[[s]])) { # Row indices
          r_ids <- which(layers == sibs[[s]][i])
          for(j in seq_along(sibs[[s]])) { # Column indices
            c_ids <- which(layers == sibs[[s]][j])
            if(norm == "degree") {
              adj_mat[r_ids, c_ids] <- sum2one(scale_rows(adj_mat[r_ids, c_ids], x = brw_attr[r_ids]))
            } else if(norm == "modified_degree") {
              adj_mat[r_ids, c_ids] <- sum2one(scale_rows(adjM = adj_mat[r_ids, c_ids], x = brw_attr[r_ids], k = k))
            } else stop("Unrecognized string for norm argument. Must be one of 'degree' or 'modified_degree'.")
            # Apply degree bias adjustment method
            if(!is.null(degree_bias)) {
              if(i == j) {
                if(sibs[[s]][i] %in% get_leaf_nodes(network_hierarchy, node = degree_bias$layers)) {
                  adj_mat[r_ids, c_ids] <- bistochastic_scaling(adj_mat[r_ids, c_ids], gamma = degree_bias$gamma)
                }
              }
            }
          }
        }
        
        if(length(sibs[[s]]) > 1) {
          # Apply crosstalk parameters [!!!!]
          ids <- which(layers %in% sibs[[s]])
          tmp_mat <- adj_mat[ids, ids]
          tmp_layers <- extract_string(rownames(tmp_mat), "\\|", 2)
          loop_ids <- which(Matrix::colSums(tmp_mat) != 0)
          for(i in loop_ids) {
            x_ids <- (tmp_mat@p[i] + 1):tmp_mat@p[i+1] # IDs of @x for column i
            current_layer <- tmp_layers[i] # layer of node associated with column i
            other_layers <- tmp_layers[tmp_mat@i[x_ids] + 1] # Other layers that this node is connected to
            if(all(other_layers != current_layer)) {
              jump_prob <- 1
            } else {
              jump_prob <- crosstalk_params[current_layer]
            }
            n_interlayer_links <- sum(unique(other_layers) != current_layer)
            if(n_interlayer_links == 0) {
              ctp <- 1
            } else {
              ctp <- ifelse(other_layers == current_layer, 1 - crosstalk_params[current_layer], jump_prob / n_interlayer_links)
            }
            tmp_mat@x[x_ids] <- tmp_mat@x[x_ids] * ctp
          }
          
          adj_mat[ids, ids] <- sum2one(tmp_mat)
        }
      }
      
      return(adj_mat)
    }
    
    #=== RECURSIVE CASE ===#
    adj_mat <- normalize_adjmat(adj_mat = adj_mat, norm = norm, k = k, layers = layers, brw_attr = brw_attr, network_hierarchy = network_hierarchy, crosstalk_params = crosstalk_params, level = level + 1, degree_bias = degree_bias, in_parallel = in_parallel, n_cores = n_cores)
    
    # For each sibling set at current level:
    #   normalize inter-unit sections of norm_adj_mat
    #   apply crosstalk parameters to sibling sets in norm_adj_mat
    
    sibs <- get_siblings(network_hierarchy, level)
    
    for(s in seq_along(sibs)) {
      # Normalize inter-layer adj mat
      for(i in seq_along(sibs[[s]])) { # Row indices
        r_ids <- which(layers %in% get_leaf_nodes(network_hierarchy, node = sibs[[s]][i]))
        for(j in seq_along(sibs[[s]])) { # Column indices
          if(i == j) next
          c_ids <- which(layers %in% get_leaf_nodes(network_hierarchy, node = sibs[[s]][j]))
          if(norm == "degree") {
            adj_mat[r_ids, c_ids] <- sum2one(scale_rows(adj_mat[r_ids, c_ids], x = brw_attr[r_ids]))
          } else if(norm == "modified_degree") {
            adj_mat[r_ids, c_ids] <- sum2one(scale_rows(adjM = adj_mat[r_ids, c_ids], x = brw_attr[r_ids], k = k))
          } else stop("Unrecognized string for norm argument. Must be one of 'degree' or 'modified_degree'.")
        }
      }
      
      if(length(sibs[[s]]) > 1) {
        # Apply crosstalk parameters [!!!!]
        leaf_nodes <- get_leaf_nodes(network_hierarchy, node = sibs[[s]])
        ids <- which(layers %in% leaf_nodes)
        tmp_mat <- adj_mat[ids, ids]
        tmp_layers <- names(leaf_nodes)[match(extract_string(rownames(tmp_mat), "\\|", 2), leaf_nodes)]
        loop_ids <- which(Matrix::colSums(tmp_mat) != 0)
        for(i in loop_ids) {
          x_ids <- (tmp_mat@p[i] + 1):tmp_mat@p[i+1] # IDs of @x for column i
          current_layer <- tmp_layers[i] # layer of node associated with column i
          other_layers <- tmp_layers[tmp_mat@i[x_ids] + 1] # Other layers that this node is connected to
          if(all(other_layers != current_layer)) {
            jump_prob <- 1
          } else {
            jump_prob <- crosstalk_params[current_layer]
          }
          n_interlayer_links <- sum(unique(other_layers) != current_layer)
          if(n_interlayer_links == 0) {
            ctp <- 1
          } else {
            ctp <- ifelse(other_layers == current_layer, 1 - crosstalk_params[current_layer], jump_prob / n_interlayer_links)
          }
          tmp_mat@x[x_ids] <- tmp_mat@x[x_ids] * ctp
        }
        
        adj_mat[ids, ids] <- sum2one(tmp_mat)
      }
    }
    
    return(adj_mat)
  }
}
