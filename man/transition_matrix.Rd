% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RandomWalk.R
\name{transition_matrix}
\alias{transition_matrix}
\title{Construct a transition matrix}
\usage{
transition_matrix(
  network,
  network_hierarchy,
  normalize = c("degree", "penalized_degree"),
  k = 0.5,
  crosstalk_params = NULL,
  degree_bias = NULL,
  in_parallel = FALSE,
  n_cores = NULL
)
}
\arguments{
\item{network}{igraph object of class 'HMNMgraph' as a result from \code{create_integrated_network()}.}

\item{network_hierarchy}{igraph object of class 'hierarchy' as a result from \code{create_network_hierarchy()}.}

\item{normalize}{Adjacency matrix normalization method to construct transition matrix.}

\item{k}{Penalization factor for normalize="penalized_degree". Must be non-negative, with larger values resulting in a greater penalty for node degree, in an effort to mitigate degree bias.}

\item{crosstalk_params}{A named numeric vector containing the crosstalk parameters for each category in network hierarchy. If NULL (default), a uniform value of 0.5 is set. Hierarchicy categories not given in \emph{crosstalk_params} will be given this default value of 0.5.}

\item{degree_bias}{A character vector or list, or NULL (default). The character vector denotes the layers to which the degree bias mitigation method will be applied. The list must contain this character vector of layers (named 'layers') and a numeric scalar (named 'gamma') between 0 and 1 denoting the strength of degree bias mitigation. The default gamma value is 0.2. Set to NULL for no degree bias adjustment.}

\item{in_parallel}{Logical. Whether to run certain operations in parallel, using the \emph{parallel}, \emph{doParallel} and \emph{foreach} packages.}

\item{n_cores}{Numeric scalar or NULL (default). Number of cores to use during parallel processing. If NULL and in_parallel=TRUE, defaults to two-thirds the number of cores detected on machine.}
}
\value{
sparse matrix of class dgCMatrix
}
\description{
Constructs a transition matrix from a given network and network hierarchy. There are various options for normalization methods and degree bias adjustment.
}
\details{
\code{transition_matrix()} first column-normalizes intra-layer and inter-layer adjacency matrices individually. Then values in \strong{crosstalk_params} are applied to inter-layer adjacency matrices to control the amount of information sharing or 'crosstalk' between layers.

Values in \strong{crosstalk_params} represent the probability of the random walker to jump from a layer of the specified category to layers of other categories in the same sibling set of the hierarchy. The network hierarchy is used such that it is easier for information to spread between layers that lie closer together in the hierarchy.

For normalize="penalized_degree", the adjacency matrix is first transformed by \code{D^(-k) \%*\% A}, where \code{A} is an adjacency matrix, \code{D} is a diagonal matrix of columns sums of \code{A}, and \code{k} is a penalization factor. This transformed matrix is then column-normalized to get a transition matrix. This is equivalent to a biased random walk which penalizes transitions to nodes as a function of node degree, with penalization factor \emph{k} controlling the strength of this inverse relationship.

The \strong{degree_bias} argument can mitigate degree bias by applying a degree bias adjustment method to specific layers in the network (see \code{\link[=bistochastic_scaling]{bistochastic_scaling()}}).
}
\examples{
# Attach igraph package
library(igraph)

# Create integrated network
g <- create_integrated_network(network_layers = simple_network)

# Create transition matrix
tmat <- transition_matrix(network = g$network,
                          network_hierarchy = g$hierarchy,
                          normalize = "penalized_degree",
                          k = 0.5)

}
\seealso{
\code{\link[=create_integrated_network]{create_integrated_network()}}, \code{\link[=module_identification]{module_identification()}}
}
