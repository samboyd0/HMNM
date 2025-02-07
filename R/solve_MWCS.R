#' @importFrom igraph vcount V V<- E E<- vertex_attr vertex_attr<-

#' @title Heuristic solution of Maximum-weight Connected Subgraph problem
#'
#' @description Given a graph and a named vector of node scores, `solve_MWCS()` heuristically finds a solution to the maximum-weight connected subgraph problem.
#'
#' @details Details are taken from the _dnet_ package documentation for the `dnetfind()` function:
#'
#'  1) transform the input graph into a new graph by collapsing connected positive nodes into a meta-node. As such, meta-nodes are isolated to each other but are linked via negative nodes (single-nodes). Clearly, meta-nodes have positive scores, and negative scores for the single-nodes.
#'
#'  2) append the weight attribute to the edges in the transformed graph. There are two types of edges: i) the single-single edge with two single-nodes as two ends, and ii) single-meta edge with a single-node as one end and a meta-node as the other end. The weight for a single-single edge is the absolute sum of the scores in its two-end single-nodes but normalised by their degrees. The weight for a single-meta edge is simply the absolute score in its single-node end normalised by the degree. As such, weights are all non-negative.
#'
#'  3) find minimum spanning tree (MST) in the weighted transformed graph. A spanning tree of the weighted graph is a subgraph that is tree and connects all the node together. The MST is a spanning tree with the sum of its edge weights minimised amongst all possible spanning trees.
#'
#'  4) find all shortest paths between any pair of meta-nodes in the MST. Within the weighted transformed graph in 2), a subgraph is induced containing nodes (only occuring in these shortest paths) and all edges between them.
#'
#'  5) within the induced subgraph, identify single-nodes that are direct neighbors of meta-nodes. For each of these single-nodes, also make sure it has the absolute scores no more than the sum of scores in its neighboring meta-nodes. These single-nodes meeting both criteria are called "linkers".
#'
#'  6) still within the induced subgraph from 5), find the linker graph(s) that contains only linkers and edges between them. This will most likely be disconnected, resulting in multiple components. Similarly to 3), find MST of the linker graph(s), called 'linker MST'. Notably, these linker MST(s) serves as the scaffold, which only contains linkers but has meta-nodes being direcly attached to.
#'
#'  7) In each linker MST plus its attached meta-nodes, find the optimal path that has the sum of scores of its nodes and attached meta-nodes maximized amongest all possible paths. Nodes along these optimal paths plus their attached meta-nodes are called 'subgraph nodes'.
#'
#'  8) For each set of 'subgraph nodes', extract a subgraph (called 'subgraph') from the input graph that only contains subgraph nodes and edges betwen them. 
#'  
#'  9) Starting with the highest scoring subgraph, and in descending order based on sum of scores for each 'subgraph', merge two subgraphs. If the resulting subgraph sum of scores is greater than previous score, keep this merged subgraph. Else, move to the next subgraph.
#'  
#'  10) Take the largest connected component of this (possibly merged) subgraph. This subgraph is the maximum scoring subgraph containing the positive nodes as many as possible, but the negative nodes as few as possible.
#'
#' @note Adapted from the _dnet_ package to be much faster, with solutions closer to optimality. Based on method from the _BioNet_ package.
#'
#' @param ig igraph
#' @param scores Named vector of scores for nodes in ig
#' @param min_cluster_size Minimum size of disconnected components of linker graph to consider when searching for maximum-weight subgraph. If NULL, only considers largest connected component of linker graph.
#'
#' @return a subnetwork as an igraph object
#'
#' @seealso [module_identification()], [RWR()]
#'
#' @examples
#' library(igraph)
#'
#' graph <- sample_pa(n = 100, power = 1.2)
#' V(graph)$name <- 1:100
#' g_scores <- rnorm(100)
#' names(g_scores) <- 1:100
#'
#' new_g <- solve_MWCS(ig = graph, scores = g_scores)
#'
#' @export
#' 
solve_MWCS <- function(ig, scores, min_cluster_size = 2) {
  if(!igraph::is_igraph(ig)) {
    stop("The function must apply to an 'igraph' object.\n")
  }
  original_classes <- class(ig)
  if(is.null(V(ig)$name)) {
    V(ig)$name <- as.character(V(ig))
  }
  if(is.null(names(scores))) {
    stop("The function must require the names of the input scores.")
  } else if(any(is.na(names(scores)))) {
    warning("Those scores with NA as names will be removed")
    scores <- scores[!is.na(names(scores))]
  }
  # Avoiding possible vertex attribute naming conflicts
  van <- "MWCS_score"
  while(van %in% igraph::vertex_attr_names(ig)) {
    van <- paste0(van, 1)
  }
  igraph::vertex_attr(ig, van) <- scores[V(ig)$name]
  pos_nodes <- which(igraph::vertex_attr(ig, van) > 0)
  if(length(pos_nodes) == 0) {
    warning("No positive nodes")
    subgraph <- igraph::make_empty_graph(n = 0, directed = FALSE)
    return(subgraph)
  } else if(length(pos_nodes) == 1) {
    subgraph <- igraph::induced_subgraph(ig, pos_nodes)
    V(subgraph)$type <- "desired"
    return(subgraph)
  } else {
    # Get the subgraphs consisting of only positive nodes
    pos_subgraph <- igraph::induced_subgraph(ig, pos_nodes)
    # Returns list of connected components of positive subgraph above. These will be the meta-nodes
    conn_comp_graph <- igraph::decompose(pos_subgraph, mode = "weak")
    # Returns the sum of the scores for each meta-node
    score_comp <- unlist(lapply(lapply(conn_comp_graph, function(x) as.numeric(igraph::vertex_attr(x, van))), sum))
    ind_order <- order(score_comp, decreasing = TRUE)
    conn_comp_graph <- conn_comp_graph[ind_order]
    score_comp <- score_comp[ind_order]
    for(i in seq_along(conn_comp_graph)) {
      conn_comp_graph[[i]]$score <- score_comp[i]
    }
    v_id <- seq(1, vcount(ig))
    names(v_id) <- V(ig)$name
    edgelist <- igraph::as_edgelist(ig, names = FALSE)
    edgelist1 <- edgelist[, 1]
    edgelist2 <- edgelist[, 2]
    # This for loop is changing the edgelist to treat the nodes in each meta-node as one node
    ig_size <- vcount(ig)
    conn_comp_gt1_ids <- which(unlist(lapply(conn_comp_graph, vcount)) > 1)
    conn_comp_eq1_ids <- which(unlist(lapply(conn_comp_graph, vcount)) == 1)
    names(conn_comp_eq1_ids) <- unlist(lapply(conn_comp_graph[conn_comp_eq1_ids], function(x) V(x)$name))
    for(i in seq_along(conn_comp_gt1_ids)) { # for each meta-node i
      # change the source, target ids of all edges containing nodes in meta-node i to be the same
      v_id_tmp <- v_id[V(conn_comp_graph[[conn_comp_gt1_ids[i]]])$name]
      edgelist1[edgelist1 %in% v_id_tmp] <- ig_size + i # new id
      edgelist2[edgelist2 %in% v_id_tmp] <- ig_size + i # new id
    }
    # Order of new_ids and new_names need to correspond to order of conn_comp_graph
    # IDs of nodes in graph that are meta nodes connected only to negative nodes.
    v_meta_id <- v_id[match(names(conn_comp_eq1_ids), names(v_id))] # v_meta_id is in same order as conn_comp_eq1_ids
    new_ids <- numeric(length(conn_comp_graph))
    new_ids[c(conn_comp_eq1_ids, conn_comp_gt1_ids)] <- c(unname(v_meta_id), seq(vcount(ig) + 1, vcount(ig) + length(conn_comp_gt1_ids)))
    # For the names cluster'x', the 'x' must correspond to the xth element of conn_comp_graph
    new_names <- paste("cluster", 1:length(conn_comp_graph), sep = "")
    names(new_ids) <- new_names
    v_id <- c(v_id[!v_id %in% v_meta_id], new_ids)
    v_name <- names(v_id)
    names(v_name) <- v_id
    # This is extracting from v_name only those nodes whose IDs appear in edgelist.
    new_edgelist <- cbind(v_name[as.character(edgelist1)],
                          v_name[as.character(edgelist2)])
    mig <- igraph::graph_from_edgelist(new_edgelist, directed = FALSE)
    mig <- igraph::simplify(mig, remove.loops = TRUE, remove.multiple = TRUE)
    node_score <- scores[V(mig)$name]
    names(node_score) <- V(mig)$name
    node_score_cluster <- sapply(conn_comp_graph, igraph::get.graph.attribute, "score")
    names(node_score_cluster) <- new_names
    ind_cluster <- grep("cluster", names(node_score))
    node_score[ind_cluster] <- node_score_cluster[names(node_score[ind_cluster])]
    # Appending edge weight attribute: absolute value of sum of degree-normalized node scores (only including negative nodes in this sum). Larger score = less desirable
    igraph::vertex_attr(mig, van) <- node_score
    score_degree <- 1/(igraph::degree(mig, mode = "all") + 1)
    tmp_score <- igraph::vertex_attr(mig, van)
    tmp_score[tmp_score > 0] <- 0
    V(mig)$score_degree <- score_degree * tmp_score
    E(mig)$weight <- rep(0, length(E(mig)))
    tmp_n1 <- igraph::as_edgelist(mig, names = FALSE)[, 1]
    tmp_n2 <- igraph::as_edgelist(mig, names = FALSE)[, 2]
    E(mig)$weight <- -(V(mig)[tmp_n1]$score_degree + V(mig)[tmp_n2]$score_degree)
    # This is to ensure the transformed graph is connected. If input graph is connected, this will always be connected.
    # If not connected, choose component with largest sum of meta-node scores
    if(!igraph::is_connected(mig, mode = "weak")) {
      decomp_graphs <- igraph::decompose(mig, mode = "weak")
      sum_pos <- lapply(decomp_graphs, function(x) {
        sum(node_score[names(which(node_score[V(x)$name] > 0))])
      })
      mig <- decomp_graphs[[which.max(sum_pos)]]
    }
    # Minimum spanning tree of transformed graph
    mst <- igraph::mst(mig, weights = E(mig)$weight)
    mst_cluster_id <- grep("cluster", V(mst)$name)
    names(mst_cluster_id) <- V(mst)[mst_cluster_id]$name
    tmp <- unlist(strsplit(names(mst_cluster_id), "cluster"))
    ttmp <- as.numeric(matrix(tmp, nrow = 2)[2, ])
    mst_cluster_id <- mst_cluster_id[order(ttmp)]
    if(length(mst_cluster_id) == 1) {
      neg_node_ids2 <- NULL
    } else {
      # Iteratively exclude all negative nodes of degree 1 from MST
      loop <- TRUE
      while(loop) {
        neg_deg1 <- which(igraph::degree(mst, mode = "all") == 1 & igraph::vertex_attr(mst, van) < 0)
        if(length(neg_deg1) == 0) {
          loop <- FALSE
        } else {
          mst <- igraph::induced_subgraph(mst, -neg_deg1)
        }
      }
      # Get induced subgraph from transformed graph with these negative leaf nodes removed
      sub_mig <- igraph::induced_subgraph(mig, which(V(mig)$name %in% V(mst)$name))
      if(!igraph::is_simple(sub_mig)) sub_mig <- igraph::simplify(sub_mig)
      # Identify linker nodes: negative nodes that have meta-node neighbors and an absolute score less than sum of meta-node neighbor scores (Step v)
      neg_node_ids <- which(igraph::vertex_attr(sub_mig, van) < 0) # negative nodes in induced subgraph of step iv
      for(i in neg_node_ids) {
        tmp_nei <- igraph::neighbors(sub_mig, v = i, mode = "all")
        tmp_nei_meta <- grep("cluster", V(sub_mig)[tmp_nei]$name) # If "cluster" is in the name, this means it's a meta node
        V(sub_mig)[i]$clusters <- list(tmp_nei[tmp_nei_meta])
      }
      score_neg_nodes <- NULL
      for(i in neg_node_ids) {
        if(!is.na(V(sub_mig)[i]$clusters[1])) { # If neg_node i has at least one meta-node neighbor...
          borders <- c(i, V(sub_mig)[i]$clusters) # vector of id of neg_node and ids of its meta-node neighbors
          borders <- unlist(borders)
          score_neg_nodes <- c(score_neg_nodes, sum(igraph::vertex_attr(sub_mig, van, borders))) # sum of scores of neg.node and its meta-node neighbors
        } else {
          score_neg_nodes <- c(score_neg_nodes, igraph::vertex_attr(sub_mig, van, i))
        }
      }
      neg_node_ids2 <- neg_node_ids[score_neg_nodes >= 0] 
      # If score_neg_nodes is negative for a neg-node, this necessarily means it wasn't connected to any meta nodes and will be removed.
    }
    # Getting MST of linker graph (containing linker nodes only) (step vi)
    if(length(neg_node_ids2) == 0) {
      tmp <- unlist(strsplit(names(node_score_cluster)[which.max(node_score_cluster)], "cluster"))
      ttmp <- as.numeric(matrix(tmp, nrow = 2)[2, ])
      tmp_nodes <- unlist(lapply(conn_comp_graph, igraph::vertex_attr, "name")[ttmp])
      subgraph <- igraph::induced_subgraph(ig, which(V(ig)$name %in% tmp_nodes))
      subgraph <- largest_connected_component(subgraph)
      V(subgraph)$type <- "desired"
      subgraph <- igraph::delete_vertex_attr(subgraph, van)
      class(subgraph) <- original_classes
      return(subgraph)
    } else {
      linker_graph <- igraph::induced_subgraph(sub_mig, neg_node_ids2) # 'linker' graph(s)
      if(!igraph::is_connected(linker_graph, mode = "weak")) {
        clust <- igraph::components(linker_graph, mode = "weak")
        max_cluster_size <- max(clust$csize)
        min_cluster_size <- min(min_cluster_size, max_cluster_size)
        clust_ids <- which(clust$csize >= min_cluster_size)
        getPathScore <- function(path, g1_score, g1_clust, g2_score) {
          s1 <- g1_score[path] # scores from a path in the MST linker graph
          tmp <- unique(unlist(g1_clust[path])) # meta nodes that linker nodes along 'path' are connected to
          s2 <- g2_score[tmp] # scores of meta nodes that linker nodes along 'path' are connected to
          sum(s1, s2)
        }
        subgraph_list <- vector("list", length(clust_ids))
        for(j in seq_along(clust_ids)) {
          linker_graph_tmp <- igraph::induced_subgraph(linker_graph, which(clust$membership == clust_ids[j]))
          if(!igraph::is_simple(linker_graph_tmp)) linker_graph_tmp <- igraph::simplify(linker_graph_tmp)
          mst_linker_graph_tmp <- igraph::mst(linker_graph_tmp, E(linker_graph_tmp)$weight) # MST of linker graph
          # Find highest scoring path (by vertex weights) in the linker MST plus attached meta-nodes (Step vii in function documentation)
          max_score <- 0
          mst_linker_graph_score_tmp <- igraph::vertex_attr(mst_linker_graph_tmp, van)
          mst_linker_graph_clust_tmp <- V(mst_linker_graph_tmp)$clusters
          sub_mig_score_tmp <- igraph::vertex_attr(sub_mig, van)
          for(i in seq_len(vcount(mst_linker_graph_tmp))) {
            path <- igraph::all_shortest_paths(mst_linker_graph_tmp, from = V(mst_linker_graph_tmp)[i], mode = "all")
            path_score <- unlist(lapply(path$res, getPathScore, g1_score = mst_linker_graph_score_tmp, g1_clust = mst_linker_graph_clust_tmp, g2_score = sub_mig_score_tmp))
            best_pos <- which.max(path_score)
            if(path_score[[best_pos]] > max_score) {
              best_path <- path$res[[best_pos]]
              max_score <- path_score[[best_pos]]
            }
          }
          if(length(best_path) != 1) {
            cluster_list <- V(mst_linker_graph_tmp)[best_path]$clusters
            names_list <- as.character(seq_along(cluster_list))
            names(cluster_list) <- names_list
            names(best_path) <- names_list
            for(i in names_list) {
              res <- lapply(cluster_list, intersect, cluster_list[[i]])
              if(length(intersect(unlist(cluster_list[as.character(which(as.numeric(names_list) < as.numeric(i)))]),
                                  unlist(cluster_list[as.character(which(as.numeric(names_list) > as.numeric(i)))]))) > 0) {
                if(length(setdiff(res[[i]], unique(unlist(res[names(res) != i])))) == 0) {
                  cluster_list <- cluster_list[names(cluster_list) != i]
                  names_list <- names_list[names_list != i]
                }
              }
            }
            best_path <- best_path[names_list]
          }
          # Induce a subgraph from input graph containing only
          #   nodes along this best path and edges between them
          pos_cluster <- V(sub_mig)[unique(unlist(V(mst_linker_graph_tmp)[best_path]$clusters))]$name
          tmp <- unlist(strsplit(pos_cluster, "cluster"))
          ttmp <- as.numeric(matrix(tmp, nrow = 2)[2, ])
          tmp_meta_nodes <- unlist(lapply(conn_comp_graph, igraph::vertex_attr, "name")[ttmp])
          tmp_border_nodes <- V(mst_linker_graph_tmp)[best_path]$name
          tmp_nodes <- c(tmp_border_nodes, tmp_meta_nodes)
          subgraph <- igraph::induced_subgraph(ig, vids = tmp_nodes)
          subgraph <- largest_connected_component(subgraph)
          type <- rep("desired", vcount(subgraph))
          names(type) <- V(subgraph)$name
          type[tmp_border_nodes[tmp_border_nodes %in% names(type)]] <- "linker"
          V(subgraph)$type <- type
          subgraph_list[[j]] <- subgraph
        }
        subgraph_list <- subgraph_list[order(unlist(lapply(subgraph_list, function(x) sum(igraph::vertex_attr(x, van)))), decreasing = TRUE)] # graphs are in descending order by sum of scores
        # IDEA (1.9.25): Somehow merge subgraphs in 'subgraph_list' that share nodes and that give a better net score than any individual subgraph.
        #   - Subgraphs will only have meta nodes in common.
        #   - For each pair, start with the highest scoring subgraph and see if merging with the other gives a higher score
        #   - If yes, merge them. If no, move to the next pair.
        g_tmp <- subgraph_list[[1]]
        score_tmp <- sum(igraph::vertex_attr(subgraph_list[[1]], van))
        ll <- ifelse(length(subgraph_list) > 1, 2, 1)
        for(l in ll:length(subgraph_list)) {
          all_nodes <- unique(c(V(g_tmp)$name, V(subgraph_list[[l]])$name))
          g_tmp2 <- igraph::induced_subgraph(ig, which(V(ig)$name %in% all_nodes))
          score_tmp2 <- sum(igraph::vertex_attr(g_tmp2, van))
          if(score_tmp2 >= score_tmp) {
            g_tmp <- g_tmp2
            score_tmp <- score_tmp2
          }
        }
        ### Take component with largest sum of scores
        if(!igraph::is_connected(g_tmp, mode = "weak")) {
          g_tmp <- getBestComponent(g_tmp, van)
        }
        # g_tmp <- largest_connected_component(g_tmp)
        g_tmp <- igraph::delete_vertex_attr(g_tmp, van)
        class(g_tmp) <- original_classes
        return(g_tmp)
      } else {
        if(!igraph::is_simple(linker_graph)) linker_graph <- igraph::simplify(linker_graph)
        mst_linker_graph <- igraph::mst(linker_graph, E(linker_graph)$weight) # MST of linker graph
        getPathScore <- function(path, g1_score, g1_clust, g2_score) {
          s1 <- g1_score[path] # scores from a path in the MST linker graph
          tmp <- unique(unlist(g1_clust[path])) # meta nodes that linker nodes along 'path' are connected to
          s2 <- g2_score[tmp] # scores of meta nodes that linker nodes along 'path' are connected to
          sum(s1, s2)
        }
        # Find highest scoring path (by vertex weights) in the linker MST plus attached meta-nodes (Step vii in dNetFind function documentation)
        max_score <- 0
        mst_linker_graph_score_tmp <- igraph::vertex_attr(mst_linker_graph, van)
        mst_linker_graph_clust_tmp <- V(mst_linker_graph)$clusters
        sub_mig_score_tmp <- igraph::vertex_attr(sub_mig, van)
        for(i in seq_len(vcount(mst_linker_graph))) {
          path <- igraph::all_shortest_paths(mst_linker_graph, from = V(mst_linker_graph)[i], mode = "all")
          path_score <- unlist(lapply(path$res, getPathScore, g1_score = mst_linker_graph_score_tmp, g1_clust = mst_linker_graph_clust_tmp, g2_score = sub_mig_score_tmp))
          best_pos <- which.max(path_score)
          if(path_score[[best_pos]] > max_score){
            best_path <- path$res[[best_pos]]
            max_score <- path_score[[best_pos]]
          }
        }
        if(length(best_path) != 1) {
          cluster_list <- V(mst_linker_graph)[best_path]$clusters
          names_list <- as.character(1:length(cluster_list))
          names(cluster_list) <- names_list
          names(best_path) <- names_list
          for(i in names_list) {
            res <- lapply(cluster_list, intersect, cluster_list[[i]])
            if(length(intersect(unlist(cluster_list[as.character(which(as.numeric(names_list) < as.numeric(i)))]),
                                unlist(cluster_list[as.character(which(as.numeric(names_list) > as.numeric(i)))]))) > 0) {
              if(length(setdiff(res[[i]], unique(unlist(res[names(res) != i])))) == 0) {
                cluster_list <- cluster_list[names(cluster_list) != i]
                names_list <- names_list[names_list != i]
              }
            }
          }
          best_path <- best_path[names_list]
        }
        # Induce a subgraph from input graph containing only
        #   nodes along this best path and edges between them
        pos_cluster <- V(sub_mig)[unique(unlist(V(mst_linker_graph)[best_path]$clusters))]$name
        tmp <- unlist(strsplit(pos_cluster, "cluster"))
        ttmp <- as.numeric(matrix(tmp, nrow = 2)[2, ])
        tmp_meta_nodes <- unlist(lapply(conn_comp_graph, igraph::vertex_attr, "name")[ttmp])
        tmp_border_nodes <- V(mst_linker_graph)[best_path]$name
        tmp_nodes <- c(tmp_border_nodes, tmp_meta_nodes)
        subgraph <- igraph::induced_subgraph(ig, vids = tmp_nodes)
        type <- rep("desired", vcount(subgraph))
        names(type) <- V(subgraph)$name
        type[tmp_border_nodes[tmp_border_nodes %in% names(type)]] <- "linker"
        V(subgraph)$type <- type
        if(!igraph::is_connected(subgraph, mode = "weak")) {
          subgraph <- getBestComponent(subgraph, van)
        }
        # subgraph <- largest_connected_component(subgraph)
        subgraph <- igraph::delete_vertex_attr(subgraph, van)
        class(subgraph) <- original_classes
        return(subgraph)
      }
    }
  }
}


#' @title Get the 'best' component of a disconnected graph
#' 
#' @description The 'best' is the component with the largest sum of a numeric vertex attribute given by 'v_attr_name'.
#' 
#' @param g igraph
#' @param v_attr_name vertex attribute name whose sum for each component will determine the
#' 
#' @return connected igraph
#'
#' @noRd
#'
getBestComponent <- function(g, v_attr_name) {
  if(!(is.character(v_attr_name) && length(v_attr_name) == 1) || !v_attr_name %in% igraph::vertex_attr_names(g) || !is.numeric(igraph::vertex_attr(g, v_attr_name))) {
    stop("'v_attr_name' must be a character string corresponding to a numeric vertex attribute of 'g'.")
  }
  clust <- igraph::components(g, mode = "weak")
  res <- numeric(length(unique(clust$membership)))
  for(i in seq_along(res)) {
    res[i] <- sum(igraph::vertex_attr(g, v_attr_name, which(clust$membership == i)))
  }
  subg <- igraph::induced_subgraph(g, which(clust$membership == which.max(res)))
  return(subg)
}




