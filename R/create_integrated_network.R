
#' @importFrom igraph V V<- E E<- 


#' @title Create an integrated network
#' 
#' @description Create an integrated network from either a list of network objects (igraph, adjacency matrix, or edgelist) or a single network object. Multilayer networks must have the __network_hierarchy__ argument specified.
#' 
#' @details 
#' `create_integrated_network()` constructs a single- or multi-layer network object from one or more distinct network objects (igraph, adjacency matrix, or edgelist). For multilayer networks, a network hierarchy must be supplied. Networks can be (un)weighted and/or (un)directed.
#' 
#' The hierarchy can have an arbitrary number of levels, where each level is associated with a factor (e.g., tissue, molecule, data type). Each leaf in the hierarchy (i.e., category with no child categories) corresponds to a layer in the multilayer network, and the nodes in a layer are defined by the categories of its ancestors in the hierarchy.
#' 
#' Layers are connected with mappings given in __bipartite_networks__. When any input network object has nodes from more than one layer, additional information is required to associate each node in that network object with a specific layer. These requirements differ depending on the network type. For igraphs, there must be a 'layer' vertex attribute containing the layer names of each node. For adjacency matrices, there must be a 'layer' attribute in the same order as rows/columns. For edgelists, there must be two additional columns, "layer1" and "layer2", containing the layer names of nodes in column 1 and 2, respectively.
#' 
#' Seed values for RWR are computed from __data__, __FUN__, and __FUN_params__. Users can also supply node-wise values in __brw_attr__ for a biased random walk, where larger values will increase transition probabilities to nodes during RWR.
#' 
#' @param network_layers,bipartite_networks Single network-like object (igraph, adjacency matrix, or edgelist) or a list of these. If a list, it must be named, with names matching category names in __network_hierarchy__. If multiple layers are contained in a single object, the list name must include these layer names separated by "|". __bipartite_networks__ should contain the mappings between different layers. Elements in __bipartite_networks__ list can be set to "common", which will connect all common nodes between the designated layers. Third column of edgelists are assumed to be edge weights unless colname "weight" is present.
#' @param network_hierarchy A 2-column matrix representing an edgelist for the network hierarchy, where hierarchy vertices represent categories which categorize the network nodes. Or an object of class 'hierarchy' as a result from `create_network_hierarchy()`.
#' @param data Named list of numeric vectors, a single numeric vector, a character string, or NULL (default). Used to calculate seeds values for RWR (with __FUN__ and __FUN_params__). Names of list should match layer names. Numeric values must be named with the corresponding node name. If a string, this should be the vertex attribute name (for igraph inputs) containing the data. NULL gives uniform seed values within each layer.
#' @param FUN Function, list of functions, or a character string denoting a default function ('binary', 'shift_scale', 'p_value', or 'exp'), to be applied to __data__ to compute seed values for RWR. Names of list must match layer names. NULL (default) applies no transformation of values in __data__. Optional function arguments given in __FUN_params__. 
#' @param FUN_params List or list of lists, containing additional function arguments for functions given in __FUN__. NULL (default) doesn't supply any additional function arguments.
#' @param directed Logical. Whether the input network should be treated as directed.
#' @param brw_attr Similar format as __data__. Contains values to be used in a biased random walk. Should contain non-negative values.
#' @param lcc Logical. Whether to take the largest connected component of the resulting network.
#' @param in_parallel Logical. Whether to run certain operations in parallel, using the _parallel_, _doParallel_ and _foreach_ packages.
#' @param n_cores Numeric scalar or NULL (default). Number of cores to use during parallel processing. If NULL and in_parallel=TRUE, defaults to two-thirds the number of cores detected on machine.
#' 
#' @return A named list:
#' * network: the integrated network of class 'igraph' and 'HMNMgraph'
#' * hierarchy: the network hierarchy of class 'igraph' and 'hierarchy'
#' 
#' @seealso [create_network_hierarchy()]
#' 
#' @examples
#' # Attach igraph package
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
#' # Construct an integrated multilayer network from 
#' #   network objects ('network_layers' arg) and bipartite mappings ('bipartite_networks' arg)
#' net <- create_integrated_network(network_layers = multilayer_network,
#'                                  bipartite_networks = interlayer_links,
#'                                  network_hierarchy = multilayer_hierarchy,
#'                                  data = multilayer_data_runif,
#'                                  FUN = "p_value")
#' class(net$network)
#' class(net$hierarchy)
#' 
#' @export
#' 
create_integrated_network <- function(network_layers, bipartite_networks = NULL, network_hierarchy = NULL, data = NULL, FUN = NULL, 
                                      FUN_params = NULL, directed = FALSE, brw_attr = NULL, lcc = FALSE, in_parallel = FALSE, n_cores = NULL) {
  #== START INPUT CHECK: network_hierarchy 
  if(is.null(network_hierarchy)) {
    if(!is.null(bipartite_networks) || (is.list(network_layers) && length(network_layers) > 1 && !igraph::is_igraph(network_layers))) stop("network_hierarchy must be supplied if bipartite_networks is non-null or network_layers is a list of several network objects, i.e., if the network inputs form a multilayer network.")
    
    # Create hierarchy matrix with default node name "layer1"
    network_hierarchy <- matrix(c("root", "layer1"), ncol = 2, byrow = TRUE)
    network_hierarchy <- create_network_hierarchy(network_hierarchy)
    
    # Format single network object as named list, with default name "layer1"
    network_layers <- list(layer1 = network_layers)
    
  } else if(!inherits(network_hierarchy, "hierarchy")) {
    network_hierarchy <- create_network_hierarchy(network_hierarchy)
  }
  
  root_node_id <- which(igraph::degree(graph = network_hierarchy, mode = "in") == 0)
  n_levels <- max(V(network_hierarchy)$level)
  #== END INPUT CHECK: network_hierarchy 
  
  
  #== START INPUT CHECK: network_layers
  if(inherits(network_layers, "HMNMgraph")) {
    if(any(!V(network_layers)$layer %in% V(network_hierarchy)$name)) stop("If network_layers is an HMNMgraph object, values of 'layer' vertex attribute must match names in network_hierarchy.")
    return(list(network = network_layers, hierarchy = network_hierarchy))
  }
  
  # Convert non-list objects to a list, with default name "layer1"
  if(!is.list(network_layers) || igraph::is_igraph(network_layers)) {
    network_layers <- list(layer1 = network_layers)
  }
  
  # Check that names of network_layers match names in network_hierarchy
  if(is.null(names(network_layers))) stop("If network_layers is input as a list of network objects, it must be named.")
  if(any(!extract_string(names(network_layers), "\\|", 0) %in% V(network_hierarchy)$name)) stop("List names of network_layers must match names in network_hierarchy.")
  
  # Get object types 
  layer_types <- sapply(network_layers, get_network_type)
  
  # Determine if each object is multilayer or single-layer
  multi_layer <- sapply(names(network_layers), function(obj_name) {
    grepl("\\|", obj_name) || V(network_hierarchy)$level[V(network_hierarchy)$name == obj_name] < n_levels
  })
  
  # Checking that network objects have necessary layer information 
  # Also, checking correct formatting of igraph, adj mat, and edgelist
  # Then, converting all non-igraph objects to igraphs with "layer", "name", and "original_ID" vertex attributes
  if(in_parallel) {
    # .Platform$OS.type # one of "unix" or "windows"
    `%dopar%` <- get("%dopar%", asNamespace("foreach"))
    nl_names <- names(network_layers)
    if(is.null(n_cores)) n_cores <- round(parallel::detectCores() * 0.66)
    cl <- parallel::makeForkCluster(n_cores, outfile = "")
    on.exit(parallel::stopCluster(cl))
    # if(!foreach::getDoParRegistered()) doParallel::registerDoParallel(cl)
    if(foreach::getDoParRegistered()) {
      doParallel::stopImplicitCluster()
      foreach::registerDoSEQ()  # Reset to sequential backend
    }
    doParallel::registerDoParallel(cl)
    network_layers <- foreach::foreach(i = seq_along(network_layers), .packages = c("igraph", "Matrix"), .verbose = FALSE, .export = NULL, .noexport = NULL) %dopar% {
      process_network_layers(obj = network_layers[[i]], 
                             obj_name = names(network_layers)[i], 
                             obj_type = layer_types[i], 
                             multi = multi_layer[i],
                             directed = directed,
                             network_hierarchy = network_hierarchy)
    }
    names(network_layers) <- nl_names
  } else {
    for(i in seq_along(network_layers)) {
      network_layers[[i]] <- process_network_layers(obj = network_layers[[i]], 
                                                    obj_name = names(network_layers)[i], 
                                                    obj_type = layer_types[i], 
                                                    multi = multi_layer[i],
                                                    directed = directed,
                                                    network_hierarchy = network_hierarchy)
    }
  }
  
  
  
  # Now, every object is an igraph (with vertex attr "layer", "name", and "original_name", and edge attr "weight"), 
  # where "name" v.attr has format "ID|layer", and "layer" attr matches a name of hierarchical leaf node
  #== END INPUT CHECK: network_layers
  
  
  #== START INPUT CHECK: bipartite_networks
  if(!is.null(bipartite_networks)) {
    if(!is.list(bipartite_networks) || is.null(names(bipartite_networks))) stop("bipartite_networks must be a named list, or null.")
    
    # Check that names of network_layers match names in network_hierarchy
    if(any(!extract_string(names(bipartite_networks), "\\|", 0) %in% V(network_hierarchy)$name)) stop("List names of bipartite_networks must match names in network_hierarchy.")
    
    # Get object types 
    bp_types <- sapply(bipartite_networks, get_bp_network_type)
    
    # Determine if each object is multilayer or single-layer
    multi_bp <- sapply(names(bipartite_networks), function(obj_name) {
      grepl("\\|", obj_name) || V(network_hierarchy)$level[V(network_hierarchy)$name == obj_name] < n_levels
    })
    
    if(any(!multi_bp)) stop("Each bipartite network must include nodes from two or more layers.")
    
    if(in_parallel) {
      bp_names <- names(bipartite_networks)
      bipartite_networks <- foreach::foreach(i = seq_along(bipartite_networks), .packages = c("igraph", "Matrix"), .verbose = FALSE, .export = NULL, .noexport = NULL) %dopar% {
        process_bipartite_networks(obj = bipartite_networks[[i]], 
                                   obj_name = names(bipartite_networks)[i], 
                                   obj_type = bp_types[i],
                                   network_layers = network_layers, 
                                   network_hierarchy = network_hierarchy)
      }
      # parallel::stopCluster(cl)
      names(bipartite_networks) <- bp_names
    } else {
      for(i in seq_along(bipartite_networks)) {
        bipartite_networks[[i]] <- process_bipartite_networks(obj = bipartite_networks[[i]], 
                                                              obj_name = names(bipartite_networks)[i], 
                                                              obj_type = bp_types[i],
                                                              network_layers = network_layers, 
                                                              network_hierarchy = network_hierarchy,
                                                              directed = directed)
      }
    }
    # At this point, all bipartite networks are in edgelist form with a third column for weights
    # All bipartite mappings have been applied between the specific network layers
    # They are ready to be added to the network layers through rbind to get the final integrated network
    
    # Combine layers and bipartite networks
    all_networks <- c(network_layers, bipartite_networks)
    
  } else {
    # Combine layers
    all_networks <- c(network_layers)
  }
  #== END INPUT CHECK: bipartite_networks
  
  
  #== Merge the separate network objects into one igraph
  #   All elements in all_networks are either igraph or edgelist
  #   All nodes follow "ID|layer" format
  #   All igraphs have vattrs: name, original_name, layer... and eattrs: weight
  #   All edgelists have 3rd weight column
  
  # Get names of all vertex atrributes in igraphs
  uniq.v.attrs <- unique(unlist(lapply(all_networks, function(x) {
    if(igraph::is_igraph(x)) {
      igraph::vertex_attr_names(x)
    } else NULL
  })))
  
  # Collect any vertex attributes from graphs
  v.attrs <- vector("list", length(uniq.v.attrs)); names(v.attrs) <- uniq.v.attrs
  for(j in seq_along(v.attrs)) {
    v.attrs[[j]] <- unlist(lapply(unname(all_networks), function(x) {
      if(igraph::is_igraph(x)) {
        if(uniq.v.attrs[j] %in% igraph::vertex_attr_names(x)) {
          res <- igraph::vertex_attr(x, uniq.v.attrs[j])
          names(res) <- igraph::vertex_attr(x, "name")
          res
        } else NULL
      } else NULL
    }))
  }
  
  # Get names of all edge atrributes in igraphs
  uniq.e.attrs <- unique(unlist(lapply(all_networks, function(x) {
    if(igraph::is_igraph(x)) {
      igraph::edge_attr_names(x)
    } else NULL
  })))
  

  # Convert all network objects to edgelists, collecting any edge attributes
  e.attrs <- vector("list", length(uniq.e.attrs)); names(e.attrs) <- uniq.e.attrs
  edge_list <- vector("list", length(all_networks))
  for(i in seq_along(edge_list)) {
    if(igraph::is_igraph(all_networks[[i]])) {
      edge_list[[i]] <- igraph::as_edgelist(graph = all_networks[[i]], names = TRUE)
      for(j in seq_along(e.attrs)) {
        if(names(e.attrs)[j] %in% igraph::edge_attr_names(all_networks[[i]])) {
          e.attrs[[j]] <- c(e.attrs[[j]], igraph::edge_attr(all_networks[[i]], names(e.attrs)[j]))
        } else e.attrs[[j]] <- c(e.attrs[[j]], rep(NA, igraph::ecount(all_networks[[i]])))
      }
    } else {
      edge_list[[i]] <- all_networks[[i]][,1:2]
      for(j in seq_along(e.attrs)) {
        if(names(e.attrs)[j] == "weight") {
          e.attrs[[j]] <- c(e.attrs[[j]], all_networks[[i]][,3])
        } else e.attrs[[j]] <- c(e.attrs[[j]], rep(NA, nrow(all_networks[[i]])))
      }
    }
  }
  # Merge to a single edgelist and convert to igraph
  edge_list <- do.call(rbind, edge_list)
  graph <- igraph::graph_from_edgelist(el = edge_list, directed = directed) 
  
  # Assign edge attributes 
  for(i in seq_along(e.attrs)) {
    if(names(e.attrs)[i] == "weight") {
      igraph::edge_attr(graph, names(e.attrs)[i]) <- as.numeric(e.attrs[[i]])
    } else {
      igraph::edge_attr(graph, names(e.attrs)[i]) <- e.attrs[[i]]
    }
  }
  # Assign vertex attributes
  for(i in seq_along(v.attrs)) {
    id <- match(V(graph)$name, names(v.attrs[[i]]))
    igraph::vertex_attr(graph, names(v.attrs)[i], which(!is.na(id))) <- unname(v.attrs[[i]][id[!is.na(id)]])
  }
  # Simplify graph
  graph <- igraph::simplify(graph = graph, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = list(weight = "median", "first"))
  
  
  #=====#
  # DATA 
  #=====#
  if(is.list(data) || is.numeric(data)) {
    if(is.null(names(data))) stop("If data is a list, it must be named.")
    
    if(is.numeric(data)) {
      data <- list(data)
      names(data) <- V(network_hierarchy)$name[root_node_id]
    } 
    
    if(any(!names(data) %in% V(network_hierarchy)$name)) stop("For data, names of list elements must match nodes in hierarchy.")
    
    for(i in seq_along(data)) {
      leaf_nodes <- get_leaf_nodes(network_hierarchy, node = names(data)[i])
      for(j in seq_along(leaf_nodes)) {
        key_names <- paste(gsub(pattern = "\\|", replacement = "_", x = names(data[[i]])), leaf_nodes[j], sep = "|")
        id <- match(V(graph)$name, key_names)
        igraph::vertex_attr(graph, "data", which(!is.na(id))) <- unname(data[[i]][id[!is.na(id)]])
      }
    }
  } else if(is.character(data) && length(data) == 1) {
    if(!data %in% igraph::vertex_attr_names(graph)) stop("For data: No vertex attribute '", data,"' found.")
    V(graph)$data <- igraph::vertex_attr(graph, data)
  } else if(is.null(data)) {
    V(graph)$data <- 1
  } else stop("Unrecognized data argument. Must be a named list, numeric vector, character string, or NULL.")
  
  
  #=========#
  # brw_attr 
  #=========#
  if(is.list(brw_attr) || is.numeric(brw_attr)) {
    if(is.null(names(brw_attr))) stop("If brw_attr is a list, it must be named.")
    
    if(is.numeric(brw_attr)) {
      brw_attr <- list(brw_attr)
      names(brw_attr) <- V(network_hierarchy)$name[root_node_id]
    } 
    
    if(any(!names(brw_attr) %in% V(network_hierarchy)$name)) stop("For brw_attr, names of list elements must match nodes in hierarchy.")
    
    for(i in seq_along(brw_attr)) {
      leaf_nodes <- get_leaf_nodes(network_hierarchy, node = names(brw_attr)[i])
      for(j in seq_along(leaf_nodes)) {
        key_names <- paste(gsub(pattern = "\\|", replacement = "_", x = names(brw_attr[[i]])), leaf_nodes[j], sep = "|")
        id <- match(V(graph)$name, key_names)
        igraph::vertex_attr(graph, "brw_attr", which(!is.na(id))) <- unname(brw_attr[[i]][id[!is.na(id)]])
      }
    }
  } else if(is.character(brw_attr) && length(brw_attr) == 1) {
    if(!brw_attr %in% igraph::vertex_attr_names(graph)) stop("For brw_attr: No vertex attribute '", brw_attr,"' found.")
    V(graph)$brw_atrr <- igraph::vertex_attr(graph, brw_attr)
  } else if(is.character(brw_attr) && length(brw_attr) > 1) {
    graph <- set_brw_attr(graph = graph, brw_attr = brw_attr, BRW_NOI = TRUE)
  } else if(is.null(brw_attr)) {
    V(graph)$brw_attr <- 1
  } else stop("Unrecognized brw_attr argument. Must be a named list, numeric vector, character string, or NULL.")
  
  # Coercing negative BRW values to 1
  if(any(igraph::vertex_attr(graph, "brw_attr") < 0)) {
    warning("'brw_attr' values must be non-negative. Coercing negative 'brw_attr' values to 1.")
    igraph::vertex_attr(graph, "brw_attr", which(igraph::vertex_attr(graph, "brw_attr") < 0)) <- 1
  }
  
  
  #===========#
  # FUN Checks
  #===========#
  if(is.null(data) || is.null(FUN)) FUN <- function(x, ...) x
  
  if(is.function(FUN) || (is.character(FUN) && length(FUN) == 1)) {
    FUN <- list(FUN)
    names(FUN) <- V(network_hierarchy)$name[root_node_id]
  }
  
  if(is.list(FUN)) {
    if(is.null(names(FUN))) stop("If FUN is a list, it must be named.")
    if(any(!names(FUN) %in% V(network_hierarchy)$name)) stop("For FUN, names of list elements must match nodes in hierarchy.")
  }else stop("Unrecognized FUN argument. Must be a function, a named list of functions/character strings, a character string, or NULL.")
  
  
  #==================#
  # FUN_params Checks
  #==================#
  if(is.null(data) || is.null(FUN_params)) {
    FUN_params <- list(list())
    names(FUN_params) <- V(network_hierarchy)$name[root_node_id]
  } 
  
  if(is.list(FUN_params)) {
    if(all(unlist(lapply(FUN_params, function(x) is.list(x))))) {
      if(any(!names(FUN_params) %in% V(network_hierarchy)$name)) stop("For FUN_params, names of list elements must match nodes in hierarchy when FUN_params contains lists of arguments to be passed to FUN.")
    } else {
      FUN_params <- list(FUN_params)
      names(FUN_params) <- V(network_hierarchy)$name[root_node_id]
    }
  } else stop("Unrecognized FUN_params argument. Must be a list, or NULL.")
  
  #====================================#
  # Applying FUN and FUN_params to data
  #====================================#
  # FUN and FUN_params are named lists, with names corresponding to nodes in network_hierarchy
  for(i in seq_along(FUN)) {
    if(is.character(FUN[[i]])) {
      f_tmp <- get_default_function(FUN[[i]])
    } else if(is.function(FUN[[i]])) {
      f_tmp <- FUN[[i]]
    } else stop("Unrecognized FUN argument. Must be a function, a named list of functions/character strings, a character string, or NULL.")
    
    leaf_nodes <- get_leaf_nodes(network_hierarchy, node = names(FUN)[i])
    for(j in seq_along(leaf_nodes)) {
      # Query FUN_params to get list of function parameters for leaf_nodes[j]
      for(l in seq_along(FUN_params)) {
        fp_leaf_nodes <- get_leaf_nodes(network_hierarchy, node = names(FUN_params)[l])
        if(leaf_nodes[j] %in% fp_leaf_nodes) {
          fp_tmp <- FUN_params[[l]]
          break
        }
      }
      
      ids <- which(V(graph)$layer == leaf_nodes[j])
      
      seeds_tmp <- do.call(f_tmp, c(list(V(graph)$data[ids]), fp_tmp))
      seeds_tmp[is.na(seeds_tmp)] <- 0
      
      if(stats::sd(seeds_tmp) == 0) {
        Z_tmp <- rep(1, length(seeds_tmp))
      } else {
        Z_tmp <- (seeds_tmp - mean(seeds_tmp)) / stats::sd(seeds_tmp) # Z-score
        # Z_tmp <- (seeds_tmp - stats::median(seeds_tmp)) / stats::mad(seeds_tmp, constant = 1) # 'robust' Z-score
        # Z_tmp <- stats::ecdf(seeds_tmp)(seeds_tmp) 
      }
      
      
      if(!all_good(seeds_tmp)) stop("For layer ", leaf_nodes[j],": Applying 'FUN' to 'data' to compute seed values resulted in either negative or all zero values.")
      
      igraph::vertex_attr(graph, "seeds", ids) <- seeds_tmp
      igraph::vertex_attr(graph, "Z", ids) <- Z_tmp
    }
  }
  
  # Take largest connected component if specified by 'lcc' argument
  if(!igraph::is_connected(graph)) {
    if(lcc) {
      graph <- largest_connected_component(graph)
    } else warning("Integrated network is disconnected.")
  }
  
  # Assigning additional class to graph for processing in other package functions
  class(graph) <- c("HMNMgraph", "igraph")
  
  return(list(network = graph, hierarchy = network_hierarchy))
}


#' @title Create a network hierarchy
#' 
#' @description Creates an igraph object with additional class "hierarchy" representing the network hierarchy.
#' 
#' @details 
#' The network hierarchy contains levels representing factors (e.g., tissue, molecule, data type), with each factor containing categories (e.g., tissue=\{lung, blood\}, molecule=\{protein, metabolite\}, data type=\{proteomic, metabolomic\}). All categories must have a unique name, even if they represent the same category at the same level, but different location (e.g., 'proteomic' and 'proteomic2' when there are multiple layers with proteomic data). The leaves in the hierarchy (i.e., categories with no children categories, bottom level) correspond to layers in the multilayer network. The nodes in each layer are defined by the categories of its ancestors in the hierarchy.
#' 
#' The purpose of the hierarchy is to control how information is spread through the multilayer network during Random Walk with Restart (RWR). Information is more readily shared between layers located closer together in the hierarchy. It is important to carefully consider the hierarchy structure in light of the research domain or experimental design.
#' 
#' For a hierarchy with M levels, the top level is level 1 and the bottom level is level M.
#' 
#' @param el a 2-column matrix representing a directed network hierarchy (col 1 = parent, col 2 = child).
#' 
#' @return an igraph object with class "hierarchy"
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
#' class(net_hier)
#' if(interactive()) {
#'  plot(net_hier)
#' }
#' 
#' @export
#' 
create_network_hierarchy <- function(el) {
  if(!is.matrix(el) || ncol(el) != 2) stop("network_hierarchy must be a 2-column matrix representing a directed edgelist (col1 = source, col2 = target).")
  
  # Convert hierarchy matrix to igraph
  network_hierarchy <- igraph::graph_from_edgelist(el = el, directed = TRUE)
  
  # Verifying that hierarchy nodes do not contain "|"
  if(any(grepl("\\|", V(network_hierarchy)$name))) stop("Hierarchy node names cannot contain '|'. This is used as a delimiter for internal processing.")
  
  # The hierarchy network must be a tree
  if(!igraph::is_tree(network_hierarchy, mode = "out")) {
    # Add a dummy 'root' node
    new_root_name <- "root"
    repeat {
      if(!new_root_name %in% V(network_hierarchy)$name) break
      new_root_name <- paste(new_root_name, sample(0:9,1), sep = "")
    }
    root_node_ids <- which(igraph::degree(graph = network_hierarchy, mode = "in") == 0)
    root_names <- V(network_hierarchy)$name[root_node_ids]
    add_to_el <- c(rbind(rep(new_root_name, length(root_names)), root_names))
    el <- rbind(el, matrix(add_to_el, ncol = 2, byrow = TRUE))
    # Convert new hierarchy matrix to igraph
    network_hierarchy <- igraph::graph_from_edgelist(el = el, directed = TRUE)
    
    if(!igraph::is_tree(network_hierarchy, mode = "out")) stop("Network hierarchy must be a tree.")
  } 
  
  root_node_id <- which(igraph::degree(graph = network_hierarchy, mode = "in") == 0)
  leaf_node_ids <- which(igraph::degree(graph = network_hierarchy, mode = "out") == 0)
  
  d <- igraph::distances(graph = network_hierarchy, v = root_node_id, weights = NA, mode = "out")
  igraph::vertex_attr(graph = network_hierarchy, name = "level") <- apply(d, 2, function(x) x[x != Inf])
  
  # Checking that the hierarchy is fully specified at all levels
  d_root2leaf <- igraph::distances(graph = network_hierarchy, v = root_node_id, to = leaf_node_ids, weights = NA, mode = "out")
  flag <- any(apply(d_root2leaf, 1, function(x){
    x <- x[x != Inf]
    all(x != max(V(network_hierarchy)$level))
  }))
  if(flag) stop("All path lengths from root node to a leaf node must be equal to the number of levels in hierarchy.\nAlong each root-to-leaf path, there must be a node specified for each level.")
  
  class(network_hierarchy) <- c("hierarchy", "igraph")
  
  return(network_hierarchy)
}


#' @title Get a built-in function from keyword
#'
#' @description
#' Return a default function used to compute seed values for random walk with restart.
#'
#' @details
#' For the input _f_, NULL means no transformation is done. 
#'
#' 'binary' coerces to numeric and then sets all positive values to 1, and to 0 otherwise. 
#'
#' 'shift_scale' is for data types with a range centered about zero and takes two parameters: _DOI_ (direction of interest: 1 for positive, -1 for negative, and 0 for both) and _w_ (numeric value between 0 and 1). It takes the absolute values of the data and then down-weights nodes that weren't in DOI by _w_ (if _DOI_=0, _w_ is coerced to 1). 
#'
#' 'p_value' assumes the data are p-values and calculates _-log10(x)_ . 
#'
#' 'exp' exponentiates the values and has a _DOI_ argument. For _DOI_=-1 or 1, it is exp(_DOI_ * _x_). For _DOI_=0, it is exp(abs( _x_)).
#'
#' @param f character string: 'binary', 'shift_scale', 'p_value', 'exp', or NULL
#'
#' @returns a function
#'
#' @examples
#' get_default_function("p_value")
#' get_default_function(NULL)
#' 
#' @export
#'
get_default_function <- function(f) {
  if(is.null(f)) {
    res <- function(x, ...) x
  } else if(f == "binary") {
    res <- function(x, ...) ifelse(as.numeric(x) > 0, 1, 0)
  } else if(f == "shift_scale") {
    res <- function(x, DOI = 1, w = 0.5, ...) {
      if(DOI == -1) {
        ifelse(x < 0, abs(x), w * x)
      } else if(DOI == 1) {
        ifelse(x > 0, x, w * abs(x))
      } else if(DOI == 0) {
        abs(x)
      } else stop("DOI must be -1, 0, or 1.")
    }
  } else if(f == "p_value") {
    res <- function(x, ...) -log(x + 1e-10, base = 10)
  } else if(f == "exp") {
    res <- function(x, DOI = 1, ...) {
      if(DOI == 0) {
        exp(abs(x))
      } else if(DOI %in% c(-1,1)) {
        exp(x * DOI)
      } else stop("DOI must be -1, 0, or 1.")
    }
  } else stop(paste0("Unrecognized character string '", f,"' for FUN."))
  return(res)
}


#' @title Process network-like objects representing layers
#' 
#' @description Check valid formatting for network objects and transform to igraph object if necessary.
#' 
#' @param obj network object (igraph, adjacency matrix, or edgelist)
#' @param obj_name name of network
#' @param obj_type type of network ('igraph', 'adjmat', or 'edgelist')
#' @param multi Logical. Whether network contains multiple layers
#' @param directed Logical. Whether the network should be treated as directed.
#' @param network_hierarchy igraph object representing network hierarchy
#' 
#' @return list of processed networks as igraph objects
#' 
#' @noRd
#'
process_network_layers <- function(obj, obj_name, obj_type, multi, directed, network_hierarchy) {
  if(multi) { # Multilayer object
    # Get leaf nodes in hierarchy downstream of the node that matches object name 
    tmp_leaf_nodes <- get_leaf_nodes(g = network_hierarchy, node = extract_string(obj_name, "\\|", 0))
    
    if(obj_type == "igraph") {
      
      if(!"name" %in% igraph::vertex_attr_names(obj)) stop("For ", obj_name, ": igraphs must have 'name' vertex attribute.")
      if(!"layer" %in% igraph::vertex_attr_names(obj)) stop("For ", obj_name, ": Multilayer igraphs must have 'layer' vertex attribute.")
      
      if(any(!tmp_leaf_nodes %in% V(obj)$layer) || any(!unique(V(obj)$layer) %in% tmp_leaf_nodes)) stop("For ", obj_name, ": Values of 'layer' vertex attribute must match hierarchy leaf nodes corresponding to the multilayer input.")
      
    } else if(obj_type == "adjmat") {
      
      # Is it a square matrix?
      if(nrow(obj) != ncol(obj)) stop("For ", obj_name, ": Adjacency matrix must be square.")
      # Row or colnames?
      if(is.null(rownames(obj))) {
        if(is.null(colnames(obj))) {
          stop("For ", obj_name, ": Adjacency matrices must have rownames and/or colnames.")
        } else rownames(obj) <- colnames(obj)
      }
      
      if(!"layer" %in% names(attributes(obj))) stop("For ", obj_name, ": Multilayer adjacency matrices must have 'layer' attribute.")
      
      if(any(!tmp_leaf_nodes %in% attr(obj, "layer")) || any(!unique(attr(obj, "layer")) %in% tmp_leaf_nodes)) stop("For ", obj_name, ": Values of 'layer' attribute must match hierarchy leaf nodes corresponding to the multilayer input.")
      
    } else if(obj_type == "edgelist") {
      if(ncol(obj) > 2 && !"weight" %in% colnames(obj)) colnames(obj)[3] <- "weight"
      
      if(ncol(obj) < 4 || any(!c("layer1", "layer2") %in% colnames(obj))) stop("For ", obj_name, ": Multilayer edgelists must have additional columns named 'layer1' and 'layer2', corresponding to the layers of nodes in the first and second columns, respectively.")
      
      if(any(!tmp_leaf_nodes %in% c(obj[,"layer1"], obj[,"layer2"])) || any(!unique(c(obj[,"layer1"], obj[,"layer2"])) %in% tmp_leaf_nodes)) stop("For ", obj_name, ": Values of 'layer1' & 'layer2' columns must match hierarchy leaf nodes corresponding to the multilayer input.")
    }
  } else { # Single-layer object
    
    if(obj_type == "igraph") {
      if(!"name" %in% igraph::vertex_attr_names(obj)) stop(paste0("For ", obj_name, ": igraphs must have 'name' vertex attribute."))
      V(obj)$layer <- obj_name
    } else if(obj_type == "adjmat") {
      # Is it a square matrix?
      if(nrow(obj) != ncol(obj)) stop(paste0("For ", obj_name, ": Adjacency matrix must be square."))
      # Row or colnames?
      if(is.null(rownames(obj))) {
        if(is.null(colnames(obj))) {
          stop(paste0("For ", obj_name, ": Adjacency matrices must have rownames and/or colnames."))
        } else rownames(obj) <- colnames(obj)
      }
      attr(obj, "layer") <- rep(obj_name, nrow(obj))
    } else if(obj_type == "edgelist") {
      if(ncol(obj) > 2 && !"weight" %in% colnames(obj)) colnames(obj)[3] <- "weight"
      obj <- cbind(obj, "layer1" = rep(obj_name, nrow(obj)), "layer2" = rep(obj_name, nrow(obj)))
    }
  }
  
  # Converting to igraph
  if(obj_type == "igraph") {
    V(obj)$original_name <- V(obj)$name
    V(obj)$name <- gsub(pattern = "\\|", replacement = "_", x = V(obj)$name)
    V(obj)$name <- paste(V(obj)$name, V(obj)$layer, sep = "|")
    if(!"weight" %in% igraph::edge_attr_names(obj)) igraph::edge_attr(obj, "weight") <- 1
  } else if(obj_type == "adjmat") {
    orig_names <- rownames(obj)
    rownames(obj) <- gsub(pattern = "\\|", replacement = "_", x = rownames(obj))
    rownames(obj) <- paste(rownames(obj), attr(obj, "layer"), sep = "|")
    obj <- igraph::graph_from_adjacency_matrix(adjmatrix = obj, mode = ifelse(directed, "directed", "undirected"), weighted = TRUE, diag = FALSE, add.rownames = TRUE)
    V(obj)$original_name <- orig_names
    V(obj)$layer <- extract_string(V(obj)$name, "\\|", 2)
  } else if(obj_type == "edgelist") {
    
    if("weight" %in% colnames(obj)) {
      obj_weights <- as.numeric(obj[,"weight"])
    } else obj_weights <- 1
    
    if(any(grepl("\\|", obj))) {
      orig_el <- obj[,1:2]
      random_delim <- "XxXA5t6V1b9xXx" # random delimiter to ensure preservation of original node names
      max_iterations <- 100 # Safety limit for attempts
      attempts <- 0
      repeat {
        if(!any(grepl(random_delim, orig_el))) break
        random_delim <- paste0(random_delim, paste(sample(c(0:9, letters, LETTERS), 5), collapse = ""))
        attempts <- attempts + 1
        if(attempts > max_iterations) stop("Failed to generate a unique delimiter after ", max_iterations, " attempts. Ensure node names do not contain '|'.")
      }
      orig_el[,1] <- paste(orig_el[,1], orig_el[,"layer1"], sep = random_delim)
      orig_el[,2] <- paste(orig_el[,2], orig_el[,"layer2"], sep = random_delim)
      
      obj[,1] <- gsub(pattern = "\\|", replacement = "_", x = obj[,1])
      obj[,2] <- gsub(pattern = "\\|", replacement = "_", x = obj[,2])
      
      obj[,1] <- paste(obj[,1], obj[,"layer1"], sep = "|")
      obj[,2] <- paste(obj[,2], obj[,"layer2"], sep = "|")
      
      obj <- igraph::graph_from_edgelist(el = obj[,1:2], directed = directed)
      V(obj)$original_name <- unique(as.character(t(orig_el[,1:2])))
    } else {
      obj[,1] <- paste(obj[,1], obj[,"layer1"], sep = "|")
      obj[,2] <- paste(obj[,2], obj[,"layer2"], sep = "|")
      
      obj <- igraph::graph_from_edgelist(el = obj[,1:2], directed = directed)
      V(obj)$original_name <- extract_string(V(obj)$name, "\\|", 1)
    }
    
    igraph::edge_attr(obj, "weight") <- obj_weights
    V(obj)$layer <- extract_string(V(obj)$name, "\\|", 2)
  }
  return(obj)
}


#' @title Process network-like objects containing inter-layer edges
#' 
#' @description Check valid formatting of network objects containing inter-layer (i.e., bipartite) edges. Transform to edgelist if necessary.
#' 
#' @param obj network object (igraph, adjacency matrix, or edgelist)
#' @param obj_name name of network
#' @param obj_type type of network ('igraph', 'adjmat', or 'edgelist')
#' @param network_layers list of network layers (as igraphs)
#' @param network_hierarchy igraph object representing network hierarchy
#' 
#' @return List of processed bipartite networks as edgelists
#' 
#' @noRd
#' 
process_bipartite_networks <- function(obj, obj_name, obj_type, network_layers, network_hierarchy, directed) {
  # Create a flag to denote if uses a keyword (e.g., "common")
  LINKING_COMMON_NODES <- obj_type == "keyword"
  
  # Get leaf nodes in hierarchy downstream of the node that matches object name 
  tmp_leaf_nodes <- get_leaf_nodes(g = network_hierarchy, node = extract_string(obj_name, "\\|", 0))
  
  if(obj_type == "igraph") {
    
    if(!"name" %in% igraph::vertex_attr_names(obj)) stop("For ", obj_name, ": igraphs must have 'name' vertex attribute.")
    if(!"layer" %in% igraph::vertex_attr_names(obj)) stop("For ", obj_name, ": Bipartite igraphs must have 'layer' vertex attribute.")
    
    obj_leaf_nodes <- get_leaf_nodes(g = network_hierarchy, node = unique(V(obj)$layer))
    
    if(any(!tmp_leaf_nodes %in% obj_leaf_nodes) || any(!obj_leaf_nodes %in% tmp_leaf_nodes)) stop("For ", obj_name, ": Values of 'layer' vertex attribute must match hierarchy nodes corresponding to the bipartite network input.")
    
    if(any(!names(tmp_leaf_nodes) %in% names(obj_leaf_nodes)) || any(!names(obj_leaf_nodes) %in% names(tmp_leaf_nodes))) stop("For ", obj_name, ": Values of 'layer' vertex attribute must match names of the bipartite network input.")
    
  } else if(obj_type == "adjmat") {
    
    # Is it a square matrix?
    if(nrow(obj) != ncol(obj)) stop("For ", obj_name, ": Adjacency matrix must be square.")
    # Row or colnames?
    if(is.null(rownames(obj))) {
      if(is.null(colnames(obj))) {
        stop("For ", obj_name, ": Adjacency matrices must have rownames and/or colnames.")
      } else rownames(obj) <- colnames(obj)
    }
    
    if(!"layer" %in% names(attributes(obj))) stop("For ", obj_name, ": Bipartite adjacency matrices must have 'layer' attribute.")
    
    obj_leaf_nodes <- get_leaf_nodes(g = network_hierarchy, node = unique(attr(obj, "layer")))
    
    if(any(!tmp_leaf_nodes %in% obj_leaf_nodes) || any(!obj_leaf_nodes %in% tmp_leaf_nodes)) stop("For ", obj_name, ": Values of 'layer' attribute must match hierarchicy nodes corresponding to the bipartite network input.")
    
  } else if(obj_type == "edgelist") {
    
    if(ncol(obj) > 2 && !"weight" %in% colnames(obj)) colnames(obj)[3] <- "weight"
    
    if(ncol(obj) < 4 || any(!c("layer1", "layer2") %in% colnames(obj))) stop("For ", obj_name, ": Bipartite edgelists must have additional columns named 'layer1' and 'layer2', corresponding to the layers of nodes in the first and second columns, respectively.")
    
    obj_leaf_nodes <- get_leaf_nodes(g = network_hierarchy, node = unique(c(obj[,"layer1"], obj[,"layer2"])))
    
    if(any(!tmp_leaf_nodes %in% obj_leaf_nodes) || any(!obj_leaf_nodes %in% tmp_leaf_nodes)) stop("For ", obj_name, ": Values of 'layer1' & 'layer2' columns must match hierarchicy nodes corresponding to the bipartite network input.")
  } else if(obj_type == "keyword") {
    obj <- connect_layers(bp_name = obj_name, network_layers = network_layers, network_hierarchy = network_hierarchy)
    obj_type <- "edgelist"
  }
  
  # Convert to edgelist
  if(obj_type == "igraph") {
    V(obj)$name <- gsub(pattern = "\\|", replacement = "_", x = V(obj)$name)
    V(obj)$name <- paste(V(obj)$name, V(obj)$layer, sep = "|")
    
    if("weight" %in% igraph::edge_attr_names(obj)) {
      obj_weights <- E(obj)$weight
    } else obj_weights <- 1
    
    obj <- igraph::as_edgelist(obj, names = TRUE)
    obj <- cbind(obj, obj_weights)
  } else if(obj_type == "adjmat") {
    rownames(obj) <- gsub(pattern = "\\|", replacement = "_", x = rownames(obj))
    rownames(obj) <- paste(rownames(obj), attr(obj, "layer"), sep = "|")
    obj <- igraph::graph_from_adjacency_matrix(adjmatrix = obj, mode = ifelse(directed, "directed", "undirected"), weighted = TRUE, diag = FALSE, add.rownames = TRUE)
    tmp_obj <- igraph::as_edgelist(obj, names = TRUE)
    tmp_obj <- cbind(tmp_obj, E(obj)$weight)
    obj <- tmp_obj
  } else if(obj_type == "edgelist") {
    if(!LINKING_COMMON_NODES) {
      obj[,1] <- gsub(pattern = "\\|", replacement = "_", x = obj[,1])
      obj[,2] <- gsub(pattern = "\\|", replacement = "_", x = obj[,2])
      obj[,1] <- paste(obj[,1], obj[,"layer1"], sep = "|")
      obj[,2] <- paste(obj[,2], obj[,"layer2"], sep = "|")
    }
    if("weight" %in% colnames(obj)) {
      obj <- obj[,c(1,2,which(colnames(obj) == "weight"))]
    } else obj <- cbind(obj[,1:2], 1)
  }
  
  # At this point, all objects will be an edgelist with a third column of weights
  # Create final bipartite mappings
  if(!LINKING_COMMON_NODES) { # BOTTLENECK
    obj <- connect_layers(obj = obj, bp_name = obj_name, network_layers = network_layers, network_hierarchy = network_hierarchy)
  }
  
  return(obj)
}


#' @title Create an edgelist of inter-layer edges
#' 
#' @description Create an edgelist of inter-layer (i.e., bipartite) edges between specific layers. `connect_layers()` makes calls to either `connect_common_nodes()` or `apply_bp_mapping()`.
#' 
#' @param obj NULL or bipartite edgelist
#' @param bp_name Name of bipartite edgelist containing the relevant edges
#' @param network_layers list of network layers (as igraphs)
#' @param network_hierarchy igraph object representing network hierarchy
#' 
#' @return an edgelist
#' 
#' @noRd
#' 
connect_layers <- function(obj = NULL, bp_name, network_layers, network_hierarchy) {
  layer_names <- extract_string(bp_name, "\\|", 0)
  leaves <- get_leaf_nodes(network_hierarchy, node = layer_names)
  
  if(!is.null(obj)) {
    c1_ids <- extract_string(obj[,1], "\\|", 1)
    c2_ids <- extract_string(obj[,2], "\\|", 1)
    c1_layers <- extract_string(obj[,1], "\\|", 2)
    c2_layers <- extract_string(obj[,2], "\\|", 2)
  }
  
  if(grepl("\\|", bp_name)) {
    if(length(extract_string(bp_name, "\\|", 0)) != 2) stop("At most two names, separated by '|', can be given in bipartite network name.")
    
    np <- expand.grid(leaves[names(leaves) == layer_names[1]], leaves[names(leaves) == layer_names[2]], stringsAsFactors = FALSE)
    names(np) <- layer_names
    
    el <- vector("list", nrow(np))
    for(i in seq_along(el)) {
      # Mapping between np[i,1] and np[i,2]
      if(is.null(obj)) {
        el[[i]] <- connect_common_nodes(layers = c(np[i,1], np[i,2]), 
                                        network_layers = network_layers, 
                                        network_hierarchy = network_hierarchy)
      } else {
        el[[i]] <- apply_bp_mapping(layers = c(np[i,1], np[i,2]), 
                                    c1_ids = c1_ids, 
                                    c2_ids = c2_ids, 
                                    c1_layers = c1_layers,
                                    c2_layers = c2_layers,
                                    weights = obj[,3], 
                                    network_layers = network_layers, 
                                    network_hierarchy = network_hierarchy)
      }
    }
  } else {
    if(length(leaves) < 2) stop("Name of bipartite mapping must correspond to at least two layers.")
    
    el <- vector("list", choose(length(leaves), 2))
    id <- 1
    for(i in 1:(length(leaves)-1)) {
      for(j in (i+1):length(leaves)) {
        # Mapping between leaves[i] and leaves[j]
        if(is.null(obj)) {
          el[[id]] <- connect_common_nodes(layers = c(leaves[i], leaves[j]), 
                                           network_layers = network_layers, 
                                           network_hierarchy = network_hierarchy)
        } else {
          el[[id]] <- apply_bp_mapping(layers = c(leaves[i], leaves[j]), 
                                       c1_ids = c1_ids, 
                                       c2_ids = c2_ids, 
                                       c1_layers = c1_layers,
                                       c2_layers = c2_layers, 
                                       weights = obj[,3], 
                                       network_layers = network_layers, 
                                       network_hierarchy = network_hierarchy)
        }
        id <- id + 1
      }
    }
  }
  el <- do.call(rbind, el)
  return(el)
}


#' @title Create an edgelist of inter-layer edges between common nodes
#' 
#' @description Create an edgelist of inter-layer (i.e., bipartite) edges between specific layers by connecting nodes with same name.
#' 
#' @param layers The layers (or categories in hierarchy) between which this mapping will be applied.
#' @param network_layers list of network layers (as igraphs)
#' @param network_hierarchy igraph object representing network hierarchy
#' 
#' @return an edgelist
#' 
#' @noRd
#' 
connect_common_nodes <- function(layers, network_layers, network_hierarchy) {
  layer_net <- layer_node <- vector("list", length(layers))
  for(i in seq_along(layer_net)) {
    id <- which(unlist(sapply(names(network_layers), function(x) layers[i] %in% get_leaf_nodes(network_hierarchy, node = extract_string(x, "\\|", 0)))))
    layer_net[[i]] <- network_layers[[id]]
    layer_node[[i]] <- unique(extract_string(V(layer_net[[i]])$name[V(layer_net[[i]])$layer == layers[i]], "\\|", 1))
  }
  
  common_nodes <- intersect(layer_node[[1]], layer_node[[2]])
  
  # if(length(common_nodes) == 0) stop("Cannot connect common nodes between non-overlapping layers: ", paste(layers, collapse = ", "), ".")
  if(length(common_nodes) == 0) {
    warning("Cannot connect common nodes between non-overlapping layers: ", paste(layers, collapse = ", "), ".")
    return(NULL)
  }
  
  if(igraph::is_directed(layer_net[[1]])) {
    el <- matrix(c(paste(common_nodes, layers[1], sep = "|"),
                   paste(common_nodes, layers[2], sep = "|"),
                   paste(common_nodes, layers[2], sep = "|"),
                   paste(common_nodes, layers[1], sep = "|")), ncol = 2)
  } else {
    el <- matrix(c(paste(common_nodes, layers[1], sep = "|"),
                   paste(common_nodes, layers[2], sep = "|")), ncol = 2)
  }
  return(el)
}


#' @title Create an edgelist of inter-layer edges from a bipartite network
#' 
#' @description Create an edgelist of inter-layer (i.e., bipartite) edges between specific layers from a pre-specified bipartite network.
#' 
#' @param layers The layers (or categories in hierarchy) between which this mapping will be applied.
#' @param c1_ids Names of nodes in column one of a bipartite edgelist.
#' @param c2_ids Names of nodes in column two of a bipartite edgelist.
#' @param c1_layers Layer names of nodes in column one of a bipartite edgelist.
#' @param c2_layers Layer names of nodes in column two of a bipartite edgelist.
#' @param weight Weights of edges in bipartite edgelist.
#' @param network_layers list of network layers (as igraphs)
#' @param network_hierarchy igraph object representing network hierarchy
#' 
#' @return an edgelist
#' 
#' @noRd
#' 
apply_bp_mapping = function(layers, c1_ids, c2_ids, c1_layers, c2_layers, weights, network_layers, network_hierarchy) {
  # Unique nodes in each of the non-BP network objects associated with layers
  layer_node <- vector("list", length(layers)); names(layer_node) <- layers
  for(i in seq_along(layer_node)) {
    id <- which(unlist(sapply(names(network_layers), function(x) layers[i] %in% get_leaf_nodes(network_hierarchy, node = extract_string(x, "\\|", 0)))))
    layer_node[[i]] <- unique(extract_string(V(network_layers[[id]])$name[V(network_layers[[id]])$layer == layers[i]], "\\|", 1))
  }
  
  c1_leafnodes <- layers[c1_layers]
  c2_leafnodes <- layers[c2_layers]
  
  valid_keys <- unlist(lapply(names(layer_node), function(x) {
    paste(layer_node[[x]], x, sep = "|")
  }))
  keys_1 <- paste(c1_ids, c1_leafnodes, sep = "|")
  keys_2 <- paste(c2_ids, c2_leafnodes, sep = "|")
  
  # Check if both keys are in `valid_keys`
  matches <- keys_1 %in% valid_keys & keys_2 %in% valid_keys
  
  # Extract results based on matches
  new_el <- matrix(c(keys_1[matches],
                     keys_2[matches],
                     weights[matches]), ncol = 3)
  
  return(new_el)
}


#' @title Get the type of a network-like object
#' 
#' @description Helper function to determine the type of a network object
#' 
#' @param obj network object (igraph, adjacency matrix, or edgelist)
#' 
#' @return character string. Either "igraph", "adjmat", or "edgelist"
#' 
#' @noRd
#' 
get_network_type <- function(obj) {
  if (inherits(obj, "igraph")) {
    type <- "igraph"
  } else if((inherits(obj, "Matrix") || is.matrix(obj)) && nrow(obj) == ncol(obj)) {
    type <- "adjmat"
  } else if((is.matrix(obj) || is.data.frame(obj)) && (ncol(obj) >= 2)) {
    type <- "edgelist"
  } else {
    stop("Invalid input type for network_layers: The input or list elements must be one of igraph, matrix, or data.frame (for edgelist).")
  }
  return(type)
}


#' @title Get the type of a network-like object for bipartite networks
#' 
#' @description Helper function to determine the type of a network object
#' 
#' @param obj network object (igraph, adjacency matrix, or edgelist)
#' 
#' @return character string. Either "igraph", "adjmat", "edgelist", or "keyword"
#' 
#' @noRd
#' 
get_bp_network_type <- function(obj) {
  keywords = c("common")
  if (inherits(obj, "igraph")) {
    type <- "igraph"
  } else if((inherits(obj, "Matrix") || is.matrix(obj)) && nrow(obj) == ncol(obj)) {
    type <- "adjmat"
  } else if((is.matrix(obj) || is.data.frame(obj)) && (ncol(obj) >= 2)) {
    type <- "edgelist"
  } else if(is.character(obj) && length(obj) == 1 && obj %in% keywords) {
    type <- "keyword"
  } else {
    stop("Invalid input type for bipartite_networks: The input or list elements must be one of igraph, matrix, data.frame (for edgelist), or character string 'common'.")
  }
  return(type)
}


#' @title Get descendant leaf nodes in hierarchy
#' 
#' @description For a given level in the hierarchy, find the descendant leaf nodes of all nodes at that level
#' 
#' @param g Hierarchy as igraph
#' @param level Numeric. Level of hierarchy to query
#' @param node Character string. Node of hierarchy to query
#' 
#' @return list of character vectors of names of leaf nodes in hierarchy for each node at the specified level.
#' 
#' @noRd
#' 
get_leaf_nodes <- function(g, level = NULL, node = NULL) {
  if(!is.null(level)) {
    nodes <- which(V(g)$level == level)
    d <- igraph::distances(graph = g, v = nodes, weights = NA, mode = 'out')
    res <- apply(d, 1, function(x){
      x <- x[x != Inf]
      names(x)[x == max(x)]
    }, simplify = FALSE)
  } else if(!is.null(node)) {
    if(is.character(node)) {
      nodes <- which(V(g)$name %in% node)
    } else {
      nodes <- node
    }
    d <- igraph::distances(graph = g, v = nodes, weights = NA, mode = 'out')
    tmp <- apply(d, 1, function(x) {
      x <- x[x != Inf]
      names(x)[x == max(x)]
    }, simplify = FALSE)
    res <- unlist(tmp, use.names = FALSE)
    names(res) <- rep(names(tmp), times = sapply(tmp, length))
  } else stop("One of level or node must be non-null.")
  return(res)
}


#' @title Extract a substring 
#' 
#' @description Extract a substring from a specified position (_pos_) after splitting on a given delimiter (_k_).
#' 
#' @param x character vector
#' @param k delimiter to split on
#' @param pos position of split to extract. 0 to extract all splits
#' 
#' @return character vector 
#' 
#' @examples
#' x <- paste(letters, LETTERS, sep = "|")
#' extract_string(x, "\\|", 2)
#' extract_string(x, "\\|", 0)
#' 
#' @export
#' 
extract_string <- function(x, k, pos) {
  if(pos == 0) {
    res <- unlist(lapply(strsplit(x, k), function(y) y[1:length(y)]))
  } else {
    res <- unlist(lapply(strsplit(x, k), function(y) y[min(length(y), pos)]))
  }
  return(res)
} 


#' @title Get largest connected component of a graph
#'
#' @description Get largest connected component of a graph. For directed graphs, returns the largest weakly connected component.
#'
#' @param g igraph object
#'
#' @returns igraph object of largest (weakly) connected component.
#' 
#' @noRd
#'
largest_connected_component <- function(g) {
  if(!igraph::is_connected(g, mode = "weak")) {
    comps <- igraph::components(g, mode = "weak")
    largest_comp_id <- which.max(comps$csize)
    g <- igraph::induced_subgraph(g, which(comps$membership == largest_comp_id))
  }
  return(g)
}


#' @title Check if values of a vector are non-negative with at least one positive value.
#'
#' @description Check if values of a vector are non-negative with at least one positive value.
#'
#' @param x numeric vector
#'
#' @returns logical: TRUE if all values are non-negative with at least one positive value.
#' 
#' @noRd
#'
all_good <- function(x) all(x >= 0) && any(x > 0)


#' @title Set the Biased Random Walk attribute
#'
#' @description Calculates values for a biased random walk. This value is a decreasing function of minimum distance in the full graph to a node of interest (NOI): f(d,k) = exp(-k*d), where k >= 0, d=distance.
#'
#' @param graph Input graph
#' @param brw_attr Character vector. Denotes the nodes of interest (NOI) used to artificially seed the active module.
#'
#' @return igraph object
#' 
#' @noRd
#'
set_brw_attr <- function(graph, brw_attr, BRW_NOI) {
  if(BRW_NOI) {
    noi <- brw_attr
    k <- 1
    # Get IDs of NOIs in graph
    node_names <- extract_string(V(graph)$name, "\\|", 1)
    noi_id <- which(node_names %in% noi)
    if(length(noi_id) > 0) {
      # Get distance to closest NOI for each node in graph
      d_mat <- igraph::distances(graph = graph, v = noi_id, mode = "all")
      node_dists <- apply(d_mat, 2, function(x) mean(x[is.finite(x)]))
      node_dists[is.na(node_dists)] <- max(node_dists, na.rm = TRUE)
      # Calculate values for biased random walk
      alpha <- exp(-node_dists * k)
      while(any(alpha == 0 & !node_dists %in% c(0, Inf))) {
        k <- k * 0.9
        alpha <- exp(-node_dists * k)
      }
      # Assign values to 'brw_attr' vertex attribute in graph
      igraph::vertex_attr(graph, "brw_attr") <- alpha[match(V(graph)$name, names(alpha))]
    } else {
      igraph::vertex_attr(graph, "brw_attr") <- 1
    }
  } 
  return(graph)
}

