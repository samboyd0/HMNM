# NOTES:
# - Users should use the function create_integrated_network() to correctly format their network and hierarchy prior to transition_matrix

#=================================================================#

#=== IMPORTS ===#
#' @importFrom igraph V 


#' @title Pipeline for running Random Walk with Restart (RWR)
#' 
#' @description `RWR_pipeline()` takes as input a mono- or multi-layer network, constructs a transition matrix, and then performs RWR. Resulting node scores can be returned as a vector or a list with elements containing scores from each layer.
#' 
#' @inheritParams create_integrated_network
#' @inheritParams RWR
#' @inheritParams transition_matrix
#' @param node_specific_restart Logical. Whether to use node-specific restart probabilities based around __restart__ and inversely proportional to node coreness.
#' @param output "list" for a list of RWR scores separated by network layer. "vector" for a single vector of all RWR scores.
#'
#' @return a numeric vector or a list of numeric vectors
#'
#' @examples
#' # Attach igraph package
#' library(igraph)
#' 
#' # Unique node names in 'simple_network'
#' uniq_names <- unique(c(simple_network[,1], simple_network[,2]))
#' 
#' # Generate random data simulating Fold Changes
#' fold_changes <- rexp(length(uniq_names))
#' names(fold_changes) <- uniq_names
#' 
#' # Run the RWR pipeline
#' res <- RWR_pipeline(network_layers = simple_network, 
#'                     data = fold_changes, 
#'                     FUN = function(x) abs(log(x)),
#'                     restart = 0.75,
#'                     output = "vector")
#' head(res)
#' 
#' 
#' @seealso [create_integrated_network()], [create_network_hierarchy()], [transition_matrix()], [RWR()]
#' 
#' @export
#' 
RWR_pipeline <- function(network_layers, bipartite_networks = NULL, network_hierarchy = NULL, data = NULL, FUN = NULL, FUN_params = NULL, directed = FALSE, brw_attr = NULL, lcc = FALSE,
                         normalize = c("degree", "modified_degree"), k = 0.5, crosstalk_params = NULL, degree_bias = NULL, restart = 0.5, seed_weights = NULL, node_specific_restart = FALSE, output = c("list", "vector")) {
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
                            normalize = normalize,
                            k = k,
                            crosstalk_params = crosstalk_params,
                            degree_bias = degree_bias)
  
  # Run RWR
  seed_vals <- V(g$network)$seeds
  names(seed_vals) <- V(g$network)$name
  ncs <- node_connectivity_score(graph = g$network, inverse = TRUE, mode = "core")
  r <- RWR(tmat = tmat, seeds = seed_vals, restart = restart, seed_weights = seed_weights, network_hierarchy = g$hierarchy, node_specific_restart = node_specific_restart, node_connectivity_scores = ncs)
  r <- r[,1]
  
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
    if(length(unique(V(g$network)$layer)) == 1) {
      names(r) <- V(g$network)$original_name
    }
    return(r)
  }
}


#' @title Perform Random Walk with Restart
#' 
#' @description 
#' Perform Random Walk with Restart (RWR) with a transition matrix, seed set(s), restart parameter, and optional seed weights.
#' 
#' @details 
#' 
#' The seed values give the probability of the random walker to start at a node and represents the prior importance of each node. This can be a single set of seed values, or _k_ sets of seed values in a _N x k_ matrix (_N_=nrow(tmat)), in which case RWR scores will be obtained for all seed sets and will return a _N x k_ matrix.
#' 
#' Seed weights give the relative importance between seed sets from different layers in a multilayer network. Seed weights can be given for each set of siblings in the network hierarchy (i.e., categories sharing the same parent). For upper-level categories of the hierarchy, the seeds weights are applied to the layers of their downstream leaf categories. This allows fine control over the relative importance of seed values coming from different categories at all levels of the hierarchy.
#' 
#' `RWR()` can implement a node-specific restart parameter based on the user-defined __restart__ parameter. For each node, __restart__ is increased or decreased as an inverse function of that nodes coreness. Therefore, if a node is located in a dense region of the network (high coreness), then its __restart__ value will decrease to enable the random walker to better explore the neighborhood around this node.
#' 
#' `RWR()` uses the power iteration method to calculate steady-state probability distributions.
#' 
#' @param tmat transition matrix inheriting class "dgCMatrix" from the "Matrix" package.
#' @param seeds matrix, data.frame, numeric vector, or NULL (default). If NULL, uniform seed values are given to all nodes.
#' @param seed_weights A list of named numeric vectors, or NULL (default). List elements should correspond to sibling sets of categories in network hierarchy. Values in each set must sum to one. NULL gives uniform values within each sibling set.
#' @param restart Restart probability. Larger values give more weight to seed values in RWR.
#' @param network_hierarchy an object of class 'hierarchy' as a result from `create_network_hierarchy()`.
#' @param node_specific_restart Logical. Whether to use node-specific restart probabilities based around __restart__ proportional to __node_connectivity_scores__.
#' @param node_connectivity_scores numeric vector from a call to `node_connectivity_score(..., inverse=TRUE)`
#' @param max_iters Maximum number of iterations for RWR.
#' 
#' @return matrix
#' 
#' @seealso [node_connectivity_score()], [RWR_pipeline()], [transition_matrix()], [module_identification()]
#' 
#' @examples
#' # Attach igraph package
#' library(igraph)
#' 
#' # Unique node names in 'simple_network'
#' uniq_names <- unique(c(simple_network[,1], simple_network[,2]))
#' 
#' # Generate random data simulating Fold Changes
#' fold_changes <- rexp(length(uniq_names))
#' names(fold_changes) <- uniq_names
#' 
#' # Create integrated network
#' g <- create_integrated_network(network_layers = simple_network, 
#'                                data = fold_changes, 
#'                                FUN = function(x) abs(log(x)))
#' 
#' # Create transition matrix
#' tmat <- transition_matrix(network = g$network,
#'                           network_hierarchy = g$hierarchy,
#'                           normalize = "modified_degree",
#'                           k = 0.5)
#'                           
#' # Run RWR
#' seed_vals <- V(g$network)$seeds
#' names(seed_vals) <- V(g$network)$name
#' r <- RWR(tmat = tmat, 
#'          seeds = seed_vals, 
#'          restart = 0.3, 
#'          network_hierarchy = g$hierarchy, 
#'          node_specific_restart = FALSE)
#' 
#' @export
#' 
RWR <- function(tmat, seeds = NULL, seed_weights = NULL, restart = 0.5, network_hierarchy, node_specific_restart = FALSE, node_connectivity_scores, max_iters = 100) {
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
    seeds <- as.matrix(seeds)
  } else if(is.numeric(seeds)) {
    seeds <- as.matrix(seeds, ncol = 1)
  } else stop("Unrecognized input for 'seeds'. Must be a matrix, data.frame, numeric vector, or NULL.")
  
  if(is.null(rownames(seeds))) {
    stop("The function requires the rownames (if matrix/data.frame) or names (if numeric) of 'seeds'.")
  } 
  if(any(is.na(rownames(seeds)))) {
    warning("seeds with NA as rownames will be removed")
    seeds <- seeds[!is.na(rownames(seeds)), , drop=FALSE]
  }
  # Add missing seeds
  if(any(!rownames(tmat) %in% rownames(seeds))) {
    new_names <- setdiff(rownames(tmat), rownames(seeds))
    zero_seeds <- matrix(0, ncol = ncol(seeds), nrow = length(new_names), dimnames = list(new_names))
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
      l <- (colSums(d0 * p_i)) / (colSums(p_i ^ 2))
      d <- d0 - l * p_i
      for(i in seq_len(max_iter_radj)) {
        if(all(d >= -restart) && all(d <= 1-restart)) break
        d[d < -restart] <- -restart
        d[d > 1-restart] <- 1 - restart
        l <- (colSums(d * p_i)) / (colSums(p_i ^ 2))
        d <- d - p_i %*% diag(l, ncol = ncol(seeds))
      }
    } else {
      d <- 0
    }
    
    for(i in seq_len(max_iters)) {
      p_new <- (1 - restart) * tmat %*% p_i + (restart + d) * seeds
      if(all(apply(p_new - p_i, 2, function(x) sum(abs(x)) < convergence_tol))) break # Convergence check
      p_i <- p_new
    }
    p_i[p_i < zero_tol] <- 0
    p_i <- sum2one(as.matrix(p_i), mode = "dense")
  }
  
  return(p_i)
}

#' @title Construct a transition matrix
#' 
#' @description Constructs a transition matrix from a given network and network hierarchy. There are various options for normalization methods and degree bias adjustment.
#' 
#' @details
#' `transition_matrix()` first column-normalizes intra-layer and inter-layer adjacency matrices individually. Then values in __crosstalk_params__ are applied to inter-layer adjacency matrices to control the amount of information sharing or 'crosstalk' between layers.
#' 
#' Values in __crosstalk_params__ represent the probability of the random walker to jump from a layer of the specified category to layers of other categories in the same sibling set of the hierarchy. The network hierarchy is used such that it is easier for information to spread between layers that lie closer together in the hierarchy.
#' 
#' For normalize="modified_degree", the adjacency matrix is first transformed by `D^(-k) %*% A`, where `A` is an adjacency matrix, `D` is a diagonal matrix of columns sums of `A`, and `k` is a penalization factor. This transformed matrix is then column-normalized to get a transition matrix. This is equivalent to a biased random walk which penalizes transitions to nodes as a function of node degree, with penalization factor _k_ controlling the strength of this inverse relationship.
#' 
#' The __degree_bias__ argument can mitigate degree bias by applying a degree bias adjustment method to specific layers in the network (see [bistochastic_scaling()]).
#' 
#' @inheritParams create_integrated_network
#' @param network igraph object of class 'HMNMgraph' as a result from `create_integrated_network()`.
#' @param network_hierarchy igraph object of class 'hierarchy' as a result from `create_network_hierarchy()`.
#' @param normalize Adjacency matrix normalization method to construct transition matrix.
#' @param k Penalization factor for normalize="modified_degree". Must be non-negative, with larger values resulting in a greater penalty for node degree, in an effort to mitigate degree bias.
#' @param crosstalk_params A named numeric vector containing the crosstalk parameters for each category in network hierarchy. If NULL (default), a uniform value of 0.5 is set. Hierarchicy categories not given in _crosstalk_params_ will be given this default value of 0.5.
#' @param degree_bias A character vector or list, or NULL (default). The character vector denotes the layers to which the degree bias mitigation method will be applied. The list must contain this character vector of layers (named 'layers') and a numeric scalar (named 'gamma') between 0 and 1 denoting the strength of degree bias mitigation. The default gamma value is 0.2. Set to NULL for no degree bias adjustment.
#' 
#' @return sparse matrix of class dgCMatrix
#' 
#' @seealso [create_integrated_network()], [module_identification()]
#' 
#' @examples
#' # Attach igraph package
#' library(igraph)
#' 
#' # Create integrated network
#' g <- create_integrated_network(network_layers = simple_network)
#' 
#' # Create transition matrix
#' tmat <- transition_matrix(network = g$network,
#'                           network_hierarchy = g$hierarchy,
#'                           normalize = "modified_degree",
#'                           k = 0.5)
#' 
#' @export
#' 
transition_matrix <- function(network, network_hierarchy, normalize = c("degree", "modified_degree"), k = 0.5,
                              crosstalk_params = NULL, degree_bias = NULL, in_parallel = FALSE, n_cores = NULL) {
  #=== Function Settings ===#
  DEFAULT_CROSSTALK <- 0.5
  DEFAULT_GAMMA <- 0.2
  
  normalize <- match.arg(normalize)
  
  #=== network_hierarchy checks ===#
  if(!inherits(network_hierarchy, "hierarchy")) stop("network_hierarchy must be an object of class 'hierarchy'.")
  
  #=== network checks ===#
  if(!inherits(network, "HMNMgraph")) stop("network must be an object of class 'HMNMgraph'.")
  
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
                           norm = normalize,
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


#' @title Bistochastic Scaling
#'
#' @description 
#' Directly modify the transition matrix to attenuate the influence of degree on diffusion scores (as evidenced by an increased entropy of stationary distribution associated with the modified transition matrix).
#' This is done by scaling the transition matrix to be approximately bistochastic (all row & column sums equal 1).
#'
#' @details
#' `bistochastic_scaling()` uses the Iterative Proportional Fitting (IPF) algorithm to scale the input matrix to have target row and column sums. The taret column sums are set to 1, while the target row sums are set to `rowsum - gamma * (rowsum - 1)`, where `rowsum` are the original row sums of __tmat__ and `gamma` is the parameter controlling the extent of adjustment. 
#'
#' @param tmat transition matrix
#' @param gamma factor controlling to what extent the transition matrix will be adjusted to mitigate degree influence on diffusion processes. Between 0 (no adjustment) and 1 (maximum adjustment).
#'
#' @returns A modified transition matrix that is approximately bistochastic
#' 
#' @examples
#' # Attach igraph package
#' library(igraph)
#' 
#' # Create integrated network
#' g <- create_integrated_network(network_layers = simple_network)
#' 
#' # Create transition matrix
#' tmat <- transition_matrix(network = g$network,
#'                           network_hierarchy = g$hierarchy,
#'                           normalize = "modified_degree",
#'                           k = 0.5)
#'                           
#' tmat2 <- bistochastic_scaling(tmat = tmat, gamma = 1)
#' mean(abs(tmat2[tmat2 != 0] - tmat[tmat != 0]))
#' 
#' tmat2 <- bistochastic_scaling(tmat = tmat, gamma = 0.1)
#' mean(abs(tmat2[tmat2 != 0] - tmat[tmat != 0]))
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


#' @title Possible seed_weight parameters
#' 
#' @description Get sets of siblings in hierarchy that can have seed weight parameters. These will be any sibling sets with two or more nodes.
#' 
#' @details This function may be useful to determine which sets of categories in your network hierarchy should be given seed weights for use in `RWR()`.
#' 
#' @param hierarchy igraph of network hierarchy
#' 
#' @return list. list element will be an empty list if no seed weight parameters are required.
#' 
#' @seealso [RWR()], [RWR_pipeline()], [module_identification()]
#' 
#' @examples
#' # Attach igraph package
#' library(igraph)
#' 
#' # Inspect 'multilayer_hierarchy' provided in HMNM package
#' multilayer_hierarchy
#' 
#' # Create network hierarchy as an igraph object with additional class 'hierarchy'
#' net_hier <- create_network_hierarchy(multilayer_hierarchy)
#' 
#' # Get the sibling sets of hierarchy
#' seed_weight_sets(net_hier)
#' 
#' @export
#' 
seed_weight_sets <- function(hierarchy) {
  sibs <- get_siblings(hierarchy)
  res <- lapply(sibs, function(x) x[sapply(x, length) > 1] )
  res <- unlist(res, recursive = FALSE)
  return(res)
}


#' @title Iterative Proportional Fitting
#'
#' @description Scales input matrix X to have row and column sums approximately equal to 1, thereby transforming X to a bistochastic (doubly stochastic) matrix.
#'
#' @param X Matrix
#' @param gamma factor controlling to what extent the transition matrix will be adjusted to mitigate degree influence on diffusion processes. Between 0 (no adjustment) and 1 (maximum adjustment).
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


#' @title Extract the diagonal elements of a matrix
#'
#' @description Extract the diagonal elements of a square matrix.
#'
#' @param X dgCMatrix sparse matrix
#'
#' @returns vector of diagonal elements of X.
#'
#' @noRd
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
#' @noRd
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
    X <- X / rep(colSums(X), each = nrow(X))
    return(X)
  }else stop("Unrecognized string for 'mode'.")
}


#' @title Normalize an adjacency matrix in a hierarchical fashion
#'
#' @description `normalize_adjmat()` efficiently constructs a transition matrix from an adjacency matrix whose nodes are categorized by a hierarchy. 
#'
#' @details This is a recursive function.
#'
#' @param adj_mat an adjacency matrix. Should of class dgCMatrix.
#' @param layers vector of all layer names in same order as rows of adj_mat.
#'
#' @returns A column-normalized matrix containing pairwise transition probabilities
#' 
#' @noRd
#'
normalize_adjmat <- function(adj_mat, norm, k, layers, brw_attr, network_hierarchy, crosstalk_params, level = 1, degree_bias, ...) {
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
  adj_mat <- normalize_adjmat(adj_mat = adj_mat, norm = norm, k = k, layers = layers, brw_attr = brw_attr, network_hierarchy = network_hierarchy, crosstalk_params = crosstalk_params, level = level + 1, degree_bias = degree_bias, ...)
  
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
#'
#' @returns A named numeric vector containing initial probabilities for RWR
#' 
#' @noRd
#'
normalize_seeds <- function(seeds, layers, network_hierarchy, seed_weights, level = 1) {
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


# Get the parents of categories (h_nodes) in a network hierarchy (g)
get_parents <- function(g, h_nodes) {
  # For given hierarchical nodes, find parents
  uniq_h_nodes <- unique(h_nodes)
  ids <- match(uniq_h_nodes, V(g)$name)
  nei <- unlist(lapply(igraph::neighborhood(graph = g, order = 1, nodes = ids, mode = 'in'), function(x) x[!x %in% ids]))
  h_node_parents <- names(nei)[match(h_nodes, uniq_h_nodes)]
  return(h_node_parents)
}


# Get all descendants, including the queried items, between categories (h_nodes) and leaf nodes (leaves) in a network hierarchy (g).
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


