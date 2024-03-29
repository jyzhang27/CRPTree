% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testing.R
\name{one_tree_all_methods}
\alias{one_tree_all_methods}
\title{Computes multiple p-values for phylogenetic association testing}
\usage{
one_tree_all_methods(
  tree,
  nperm_obs = 500,
  nperm_null = 500,
  nperm_labels = 499,
  alt_methods = FALSE,
  path = ""
)
}
\arguments{
\item{tree}{A \code{phylo} object: a ranked, planar, partially labeled tree shape}

\item{nperm_obs}{Integer: number of planar permutations of the given tree}

\item{nperm_null}{Integer: number of simultaneous label planar permutations of the given tree}

\item{nperm_labels}{Integer: number of label permutations of the given tree}

\item{alt_methods}{Boolean: whether to run for additional methods}

\item{path}{If given, path to save observed distribution of S and null distributions of S and mu_hat}
}
\value{
Named vector containing the p-values and observed mu_hat
}
\description{
This function computes the p-value for phylogenetic association testing
Default: average S p-value, tree-average p-value
Additional: parsimony score, association index
}
