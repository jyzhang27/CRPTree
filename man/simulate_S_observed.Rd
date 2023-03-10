% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testing.R
\name{simulate_S_observed}
\alias{simulate_S_observed}
\title{Simulate observed distribution of S for a tree}
\usage{
simulate_S_observed(tree, nperm)
}
\arguments{
\item{tree}{A \code{phylo} object: a ranked, planar, partially labeled tree shape}

\item{nperm}{Integer: number of planar permutations of the given tree}
}
\value{
A vector of S values
}
\description{
This function simulates the observed distribution of S by randomly sampling planarities
}
