#' Computes p-values for a posterior distribution of trees
#'
#' This function computes S_mean, pval_tree, and pval_avg for a list of trees
#' It also processes the trees into desired format given a term or matrix indicating
#' the categories of tip labels
#'
#' @param tree_list A list of ranked, planar, partially labeled trees
#' @param tree_term String specifying all tip labels containing this term are one category, and the rest form the other category
#' @param tree_tip_correspondence Matrix with two columns: first column is the tip labels of the tree, second is the category
#' @param nperm_planar Integer: number of planar permutations of the given tree
#' @param nperm_labels Integer: number of label permutations of the given tree
#'
#' @return A matrix with the values of S_mean, pval_tree, and pval_avg for each tree
#' @export
posterior_trees_pval <- function(tree_list, tree_term, tree_tip_correspondence, nperm_planar=500, nperm_labels=499) {
  tree_list_processed <- lapply(tree_list, process_tree, tree_term, tree_tip_correspondence)
  ntrees <- length(tree_list_processed)
  tree_list_post <- t(parallel::mcmapply(one_tree_all_methods, tree_list_processed, rep(nperm_planar, ntrees),
                                         rep(nperm_planar, ntrees), rep(nperm_labels, ntrees), mc.cores=8))
  colnames(tree_list_post) <- c('S_mean', 'pval_tree', 'pval_avg')
  return(tree_list_post)
}

#' Applies BaTS on a posterior distribution of trees
#'
#' This function computes S_mean, pval_tree, and pval_avg for a list of trees
#' It also processes the trees into desired format given a term (only term)
#'
#' @param tree_list A list of trees, usually a posterior distribution of trees
#' @param tree_term String specifying all tip labels containing this term are one category, and the rest form the other category
#' @param nperm_obs Integer: number of planar permutations of the given tree to get S_mean
#' @param n_bats Integer: number of label permutations of the given tree
#'
#' @return A vector with the median values of S_mean across BaTS
#' @export
posterior_trees_bats <- function(tree_list, tree_term, n_bats = 500, nperm_obs = 500) {
  trees_with_label <- lapply(tree_list, process_tree)
  tip_labels <- trees_with_label[[1]]$tip.label
  tip_labels_binary <- ifelse(grepl(tree_term, tip_labels), -1, -2)
  tip_labels_matched <- cbind(tip_labels, tip_labels_binary)
  ntips <- length(tip_labels)

  tree_permutations <- t(replicate(n_bats, sample(1:length(tip_labels))))
  tree_bats <- parallel::mclapply(trees_with_label, one_tree_bats, tip_labels_base = tip_labels,
                        tip_labels_binary_base= tip_labels_binary, label_permutations= tree_permutations, mc.cores = 8)
  tree_bats <- simplify2array(tree_bats)
  tree_bats_medians <- apply(tree_bats, 1, median)
  return(tree_bats_medians)
}


###############################################################
# Helpers

one_tree_bats <- function(tree, tip_labels_base, tip_labels_binary_base, label_permutations, nperm_obs=500) {
  N <- tree$Nnode+1
  edge_matrix <- tree$edge
  all_cherries <- which(tabulate(edge_matrix[, 1][edge_matrix[, 2] <= N]) == 2 )
  permuting_nodes <- setdiff((N+1):(2*N-1), all_cherries)
  n_possible_permute_nodes <- length(permuting_nodes)

  tree_tip_labels <- tree$tip.label
  tree_binary_tip_labels <- tip_labels_binary_base[match(tree_tip_labels, tip_labels_base)]
  tree_label_permutations_all <- apply(label_permutations, 2, function(x) {tree_binary_tip_labels[x]})

  S_matrix <- matrix(nrow=nperm_obs, ncol=nrow(tree_label_permutations_all))
  for (i in 1:nperm_obs) {
    select_nodes <- rbinom(N-1, 1, 0.5) # choose the nodes to switch
    select_nodes[permuting_nodes - N] <- 0 # set the nodes with cherry to 0
    select_nodes <- which(select_nodes == 1) + N

    tree_new_planar <- reorder_branches(tree, select_nodes)
    tree_new_attachments <- compute_attachment_matrix(tree_new_planar, rep(0,N-1))
    tree_tip_order <- match(tree_new_planar$tip.label, tree_tip_labels)

    labels_new_all <- tree_label_permutations_all[, tree_tip_order]

    S_matrix[i,] <- apply(labels_new_all, 1, function(x){compute_S_given_attachment_matrix(tree_new_attachments, x)})

  }
  S0 <- colMeans(S_matrix)
  return(S0)
}
