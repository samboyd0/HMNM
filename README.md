
<!-- README.md is generated from README.Rmd. Please edit that file -->

# HMNM

<!-- badges: start -->
<!-- badges: end -->

The *HMNM* package (Hierarchical Multilayer Network Models) provides
functionality for constructing and analyzing hierarchical multilayer
networks. A multilayer network consists of distinct layers with nodes or
edges that differ from other layers by some characteristic.
Functionalities of *HMNM* include multilayer network construction,
network diffusion (Random Walk with Restart), and module identification.

In *HMNM* , each layer is defined by categories coming from one or more
factors that are organized hierarchically. The purpose of the hierarchy
is to govern how information flows between layers during diffusion
processes, namely Random Walk with Restart (RWR); information is shared
more easily between layers located closer in the hierarchy. This
hierarchy can have an arbitrary number of levels that correspond to
factors, where factors can describe node/edge types or node/edge
attribute types. One guiding principle for hierarchy design is to set
the factors for each level such that more similar layers are closer
together in the hierarchy.

The two main analysis methods available in *HMNM* are RWR and the AMEND
algorithm (Active Module identification with Experimental data and
Network Diffusion), both generalized for hierarchical multilayer
networks. Random Walk with Restart (RWR) is a powerful network diffusion
method that returns node scores reflecting proximity in the network to
other high scoring nodes. Node scores can be used for node ranking, link
prediction, or in machine learning algorithms. AMEND is an active module
identification method that finds a single, connected module of densely
connected nodes with large experimental values (e.g., ’omic
experiments). It uses RWR to create node weights, which are used to find
a maximum-weight connected subnetwork, with this subnetwork used as
input into RWR for the next iteration. This is continues until an
optimal subnetwork (i.e., module) is found. AMEND can accommodate
multilayer networks, making it a widely applicable tool. These complex
networks can include several node types, edge types, and data types.

## Installation

You can install *HMNM* from [GitHub](https://github.com/samboyd0/HMNM)
with:

``` r
devtools::install_github("samboyd0/HMNM", build_vignettes = TRUE)
```

## Vignette

A vignette is available that illustrates how to use the *HMNM* package.
It can be accessed with the following code.

``` r
vignette("use_HMNM", package = "HMNM")
```

## Author(s)

Sam Boyd (Channing Division of Network Medicine, Brigham & Women’s
Hospital, Harvard Medical School)

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/master/examples>. -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
