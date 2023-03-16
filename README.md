# CRPTree
Phylogenetic trait association for binary traits. [arxiv paper link](https://arxiv.org/abs/2302.09402)

## Installation

Install using `devtools::install_github("jyzhang27/CRPTree")`. We don't recommend `build_vignettes = TRUE` due to the speed of some computation. 

## Basic Usage

If you would just like to run phylogenetic trait association for your given tree (`phylo` object from `ape`), we require there to be no ties in the internal nodes (ie no two internal nodes are the same distance to the root). The binary label can either be directly specified in `tip.label` of the `phlyo` object, or can be provided by 
1. `term`: a string such that all tip labels containing `term` are one category.
2. `tip_corresponding`: a matrix with two columns, column 1 containing the tip labels of the tree, and column 2 containing values -1, -2 indicating the category. 

Then, applying `tree <- process_tree(tree, term)` or `tree <- process_tree(tree, tip_corresponding)` will convert any `phylo` object into the correct format. Other functions such as `compute_mu_hat(), one_tree_all_methods()` can be applied accordingly. See Vignettes for more detail. 

## Vignettes

1. [Fixed Trees](https://github.com/jyzhang27/CRPTree/blob/main/vignettes/crp_tree.Rmd): A vignette detailing how to generate trees from the CRP-Tree model and run our test fixed trees. Several enumerative methods are also given. 
2. [Bayesian Examples](https://github.com/jyzhang27/CRPTree/blob/main/vignettes/Bayesian.Rmd): A vignette showing how to run our tests on a posterior sample of trees. The provided example is from [1]. 

Note: In the vignettes, there are parameters for number of simulations that are set to smaller numbers for ease of runtime. If you would like to run our methods, leaving these as default will suffice. 

## References 

1. N.R. Faria et al.
[Genomic and epidemiological monitoring of Yellow Fever virus transmission potential](https://www.science.org/doi/10.1126/science.aat7115).
*Science*, 361(6405):894–899, 2018.
