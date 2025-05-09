% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_integrated_network.R
\name{create_network_hierarchy}
\alias{create_network_hierarchy}
\title{Create a network hierarchy}
\usage{
create_network_hierarchy(el)
}
\arguments{
\item{el}{a 2-column matrix representing a directed network hierarchy (col 1 = parent, col 2 = child).}
}
\value{
an igraph object with class "hierarchy"
}
\description{
Creates an igraph object with additional class "hierarchy" representing the network hierarchy.
}
\details{
The network hierarchy contains levels representing factors (e.g., tissue, molecule, data type), with each factor containing categories (e.g., tissue=\{lung, blood\}, molecule=\{protein, metabolite\}, data type=\{proteomic, metabolomic\}). All categories must have a unique name, even if they represent the same category at the same level, but different location (e.g., 'proteomic' and 'proteomic2' when there are multiple layers with proteomic data). The leaves in the hierarchy (i.e., categories with no children categories, bottom level) correspond to layers in the multilayer network. The nodes in each layer are defined by the categories of its ancestors in the hierarchy.

The purpose of the hierarchy is to control how information is spread through the multilayer network during Random Walk with Restart (RWR). Information is more readily shared between layers located closer together in the hierarchy. It is important to carefully consider the hierarchy structure in light of the research domain or experimental design.

For a hierarchy with M levels, the top level is level 1 and the bottom level is level M.
}
\examples{
# Attach igraph package
library(igraph)

# Inspect 'multilayer_hierarchy' provided in HMNM package
multilayer_hierarchy

# Create network hierarchy as an igraph object with additional class 'hierarchy'
net_hier <- create_network_hierarchy(multilayer_hierarchy)
class(net_hier)
if(interactive()) {
 plot(net_hier)
}

}
