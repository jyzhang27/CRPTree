# CRPTree
Phylogenetic trait association for binary traits. [arxiv paper link](https://arxiv.org/abs/2302.09402)

## Installation

Install using `devtools::install_github("jyzhang27/CRPTree")`. We don't recommend `build_vignettes = TRUE` due to the speed of some computation. 

## Vignettes

1. [Fixed Trees](https://github.com/jyzhang27/CRPTree/blob/main/vignettes/crp_tree.Rmd): A vignette detailing how to generate trees from the CRP-Tree model and run our test fixed trees. Several enumerative methods are also given. 
2. [Bayesian Examples](https://github.com/jyzhang27/CRPTree/blob/main/vignettes/Bayesian.Rmd): A vignette showing how to run our tests on a posterior sample of trees. The provided example is from [1]. 

Note: In the vignettes, there are parameters for number of simulations that are set to smaller numbers for ease of runtime. If you would like to run our methods, leaving these as default will suffice. 

## References 

1. N.R. Faria et al.
[Genomic and epidemiological monitoring of Yellow Fever virus transmission potential](https://www.science.org/doi/10.1126/science.aat7115).
*Science*, 361(6405):894–899, 2018.
