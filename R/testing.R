#' Process a tree into a workable format
#'
#' This function converts a given tree into Newick format, with two tip labels, isochronous tips, and ranked internal nodes
#'
#' @param tree A \code{phylo} object
#' @param term String specifying all tip labels containing this term are one category, and the rest form the other category
#' @param tip_corresponding Matrix with two columns: first column is the tip labels of the tree, second is the category
#'
#' @return A \code{phylo} object of the desired format
#' @export
process_tree <- function(tree, term='', tip_corresponding=NULL) {
  tip_labels <- tree$tip.label
  unique_labels <- unique(tip_labels)

  # If the tip labels have names as the category indicator
  if (nchar(term) >0) {
    tree$tip.label <- ifelse(grepl(term, tip_labels), -1, -2)
  }

  # If there are two categories, ensure the labels are -1 and -2
  if (nchar(term) == 0 & length(unique_labels) == 2 & length(setdiff(unique_labels, c(-2,-1))) != 0) {
    tree$tip.label <-  replace(tip_labels, tip_labels == unique_labels[1], -1)
    tree$tip.label <-  replace(tip_labels, tip_labels == unique_labels[2], -2)
  }

  if (!is.null(tip_corresponding)) {
    if (nrow(tip_corresponding) != length(tip_labels) | ncol(tip_corresponding) != 2) {
      stop('Invalid tip label-category correspondence')
    }

    index <- match(tip_labels, tip_corresponding[,1])
    tree$tip.label <- tip_corresponding[index,2]
  }

  # Change to Newick form which is the only workable formula
  tree <- ape::read.tree(text=TreeTools::NewickTree(tree))

  # Extend the tips to be isochronous
  distances <- castor::get_all_distances_to_root(tree)
  tree <- castor::extend_tree_to_height(tree, new_height = ceiling(max(distances)))$tree

  # check there are no ties in internal nodes
  N <- tree$Nnode + 1
  distances <- castor::get_all_distances_to_root(tree)[(N+1): (2*N-1)]
  if (length(unique(distances)) < N-1) {
    stop('Ties exist in internal nodes')
  }

  return(tree)
}


### Actual testing methods

#' Simulate null distribution of S for a tree
#'
#' This function simulates the null distribution of S by randomly sampling planarities and labelings simultaneously
#'
#' @param tree A \code{phylo} object: a ranked, planar, partially labeled tree shape
#' @param nperm Integer: number of permutations of the given tree
#'
#' @return A vector of S values
#' @export
simulate_S_tree_null <- function(tree, nperm) {
  cherries_children <- get_cherries_and_children(tree)
  all_cherries <- cherries_children$cherries
  all_children <- cherries_children$children

  null_S_tree <- replicate(nperm, simulate_S_tree_one(tree, all_cherries, all_children))
  return(null_S_tree)
}

#' Simulate observed distribution of S for a tree
#'
#' This function simulates the observed distribution of S by randomly sampling planarities
#'
#' @param tree A \code{phylo} object: a ranked, planar, partially labeled tree shape
#' @param nperm Integer: number of planar permutations of the given tree
#'
#' @return A vector of S values
#' @export
simulate_S_observed <- function(tree, nperm) {
  N <- tree$Nnode+1
  no_permuting_nodes <- get_internal_nodes_cherry_same(tree)
  permuting_nodes <- setdiff((N+1):(2*N-1), no_permuting_nodes)
  n_possible_permute_nodes <- length(permuting_nodes)
  max_possible <- 2^n_possible_permute_nodes

  if (nperm >= max_possible) {
    s_stat <- vector(mode='integer', length=max_possible)
    for (i in 1:max_possible) {
      nodes_list <- as.integer(unlist(strsplit(R.utils::intToBin(i-1), split = "")))
      nodes_list <- c(rep(0, n_possible_permute_nodes - length(nodes_list)), nodes_list)
      select_nodes <- rep(0, N-1)
      select_nodes[permuting_nodes - N] <- nodes_list
      select_nodes <- which(select_nodes == 1) + N

      s_stat[i] <- count_same_attachments(tree, select_nodes)
    }
  } else {
    s_stat <- vector(mode='integer', length=nperm)
    for (i in 1:nperm) {
      select_nodes <- rbinom(N-1, 1, 0.5) # choose the nodes to switch
      select_nodes[no_permuting_nodes - N] <- 0 # set the nodes with cherry same to 0
      select_nodes <- which(select_nodes == 1) + N
      s_stat[i] <- count_same_attachments(tree, select_nodes)
    }
  }
  return(s_stat)
}

#' Estimates mu_hat for a tree
#'
#' This function estimates the mean number of attachments for a given tree shape
#' with a set number of permutations
#'
#' @param tree A \code{phylo} object: a ranked, planar, partially labeled tree shape
#' @param nperm Integer: number of planar permutations of the given tree
#'
#' @return An integer with the value of mu_hat
#' @export
compute_mu_hat <- function(tree, nperm) {
  return(mean(simulate_S_observed(tree, nperm)))
}

#' Simulates S for a grid of planarities and labels
#'
#' This function samples S for a grid of planarities and labels where
#' one row is one planar representation and one column is a partial labeling
#'
#' @param tree A \code{phylo} object: a ranked, planar, partially labeled tree shape
#' @param nperm_planar Integer: number of planar permutations of the given tree
#' @param nperm_labels Integer: number of label permutations of the given tree
#'
#' @return A matrix with the S values
#' @export
simulate_S_grid <- function(tree, nperm_planar, nperm_labels) {
  N <- tree$Nnode+1
  edge_matrix <- tree$edge
  all_cherries <- which(tabulate(edge_matrix[, 1][edge_matrix[, 2] <= N]) == 2 )
  permuting_nodes <- setdiff((N+1):(2*N-1), all_cherries)
  n_possible_permute_nodes <- length(permuting_nodes)
  max_possible_planar <- 2^n_possible_permute_nodes

  tree_nolabel <- tree
  tree_nolabel$tip.label <- 1:N

  tree_blues <- sum(tree$tip.label== -1)
  if (nperm_labels > choose(N, tree_blues)) {
    labels_all <- arrangements::permutations(freq = c(tree_blues, N-tree_blues), k= N, x= c(-1, -2))
  } else {
    labels_all <- t(replicate(nperm_labels, sample(tree$tip.label)))
  }

  if (nperm_planar > max_possible_planar) {
    S_matrix <- matrix(nrow=max_possible_planar, ncol=nrow(labels_all))

    for (i in 1:max_possible_planar) {
      nodes_list <- as.integer(unlist(strsplit(R.utils::intToBin(i-1), split = "")))
      nodes_list <- c(rep(0, n_possible_permute_nodes - length(nodes_list)), nodes_list)
      select_nodes <- rep(0, N-1)
      select_nodes[permuting_nodes - N] <- nodes_list
      select_nodes <- which(select_nodes == 1) + N

      tree_newlabel <- reorder_branches(tree_nolabel, select_nodes)
      tree_tip_order <- tree_newlabel$tip.label
      labels_new_all <- labels_all[, tree_tip_order]

      tree_new_attachments <- compute_attachment_matrix(tree_newlabel, rep(0,N-1))

      S_matrix[i,] <- apply(labels_new_all, 1, function(x){compute_S_given_attachment_matrix(tree_new_attachments, x)})
    }
  } else {
    S_matrix <- matrix(nrow=nperm_planar, ncol=nrow(labels_all))
    for (i in 1:nperm_planar) {
      select_nodes <- rbinom(N-1, 1, 0.5) # choose the nodes to switch
      select_nodes[all_cherries - N] <- 0 # set the nodes with cherry to 0
      select_nodes <- which(select_nodes == 1) + N

      tree_newlabel <- reorder_branches(tree_nolabel, select_nodes)
      tree_tip_order <- tree_newlabel$tip.label
      labels_new_all <- labels_all[, tree_tip_order]

      tree_new_attachments <- compute_attachment_matrix(tree_newlabel, rep(0,N-1))
      S_matrix[i,] <- apply(labels_new_all, 1, function(x){compute_S_given_attachment_matrix(tree_new_attachments, x)})

    }
  }
  return(S_matrix)
}

#' Simulate null distribution of mu_hat
#'
#' This function simulates the null distribution of mu_hat for a given tree
#'
#' @param tree A \code{phylo} object: a ranked, planar, partially labeled tree shape
#' @param nperm_planar Integer: number of planar permutations of the given tree
#' @param nperm_labels Integer: number of label permutations of the given tree
#'
#' @return A vector with the mu_hat values
#' @export
simulate_mu_hat_tree <- function(tree, nperm_planar, nperm_labels) {
  S_matrix <- simulate_S_grid(tree, nperm_planar, nperm_labels)
  mu_hat <- colMeans(S_matrix)
  return(mu_hat)
}

#' Computes multiple p-values for phylogenetic association testing
#'
#' This function computes the p-value for phylogenetic association testing
#' Default: average S p-value, tree-average p-value
#' Additional: parsimony score, association index
#'
#' @param tree A \code{phylo} object: a ranked, planar, partially labeled tree shape
#' @param nperm_obs Integer: number of planar permutations of the given tree
#' @param nperm_null Integer: number of simultaneous label planar permutations of the given tree
#' @param nperm_labels Integer: number of label permutations of the given tree
#' @param alt_methods Boolean: whether to run for additional methods
#' @param path If given, path to save observed distribution of S and null distributions of S and mu_hat
#'
#' @return Named vector containing the p-values and observed mu_hat
#' @export
one_tree_all_methods <- function(tree, nperm_obs=500, nperm_null=500, nperm_labels=499, alt_methods=FALSE, path=''){
  S_obs <- simulate_S_observed(tree, nperm_obs)
  mu_hat <- mean(S_obs)

  S_tree_null <- simulate_S_tree_null(tree, nperm_null)
  mu_hat_null <- simulate_mu_hat_tree(tree, nperm_null, nperm_labels)

  pval_avg <- (1+ sum(mu_hat_null >= mu_hat))/(1+length(mu_hat_null))
  pval_tree <- mean(rbind(outer(S_tree_null, S_obs, '>='), rep(TRUE, length(S_obs))))

  if (nchar(path) >0) {
    tree_name <- deparse(substitute(tree))
    write.csv(S_obs, paste(path, tree_name, '_S_obs.csv', sep=''))
    write.csv(S_tree_null, paste(path, tree_name, '_S_tree_null.csv', sep=''))
    write.csv(mu_hat_null, paste(path, tree_name, '_mu_hat_null.csv', sep=''))
  }

  if (alt_methods) {
    pval_ai <- tree_ai(tree, nperm_labels)$pval
    pval_ps <- tree_fitch(tree, nperm_labels)$pval

    result <- c(mu_hat, pval_tree, pval_avg, pval_ai, pval_ps)
    names(result) <- c('mu_hat', 'pval_tree', 'pval_avg', 'pval_ai', 'pval_ps')
    return(result)
  } else {
    result <- c(mu_hat, pval_tree, pval_avg)
    names(result) <- c('mu_hat', 'pval_tree', 'pval_avg')
    return(result)
  }
}

##########################################################################
## Helpers

# Simulate 1 S given ranked tree shape
# will sample labels and planar representation
simulate_S_tree_one <- function(tree, cherries, children) {
  N <- tree$Nnode + 1
  tree$tip.label <- sample(tree$tip.label)
  children_colors <- matrix(tree$tip.label[children], ncol=2, byrow=TRUE)
  same_cherries <- cherries[which(children_colors[,1] == children_colors[,2])]

  select_nodes <- rbinom(N-1, 1, 0.5)
  select_nodes[same_cherries - N] <- 0 # set the nodes with cherry same to 0
  select_nodes <- which(select_nodes == 1) + N

  return(count_same_attachments(tree, select_nodes))
}

## Functions for the other methods

one_tree_fitch <- function(tree, observed) {
  N <- tree$Nnode + 1
  tree_temp <- tree
  tree_temp$tip.label <- 1:N
  tips <- tree$tip.label
  if (observed) {
    tip_mat <- matrix(tips, ncol=1, byrow=TRUE)
  } else {
    tip_mat <- matrix(sample(tips), ncol=1, byrow=TRUE)
  }
  row.names(tip_mat) <- as.character(1:N)
  tree_phyDat <- phangorn::phyDat(tip_mat, type="USER", levels= unique(as.vector(tip_mat)))
  return(fitch(tree_temp, tree_phyDat))
}

# Parsimony score (fitch algorithm) via phangorn
# low is stronger
tree_fitch <- function(tree, nperm) {
  null <- replicate(nperm, one_tree_fitch(tree, FALSE))
  observed <- one_tree_fitch(tree, TRUE)
  pval <- (1+ sum(null<=observed))/(nperm+1)
  results <- list(observed, null, pval)
  names(results) <- c('observed', 'null', 'pval')
  return(results)
}

# Association index: low is stronger
one_tree_ai <- function(tree, observed) {
  if (!observed) {
    tree$tip.label <- sample(tree$tip.label)
  }
  all_subtrees <- castor::get_subtrees_at_nodes(tree, 1:(tree$Nnode))$subtrees
  num_tips <- sapply(all_subtrees, function(x){x$Nnode+1})
  prop_tips <- sapply(all_subtrees, function(x){max(table(x$tip.label))/ (x$Nnode+1)})

  return(sum((1-prop_tips)/(2^num_tips -1)))
}

tree_ai <- function(tree, nperm) {
  null <- replicate(nperm, one_tree_ai(tree, FALSE))
  observed <- one_tree_ai(tree, TRUE)
  pval <- (1+ sum(null<=observed))/(nperm+1)
  results <- list(observed, null, pval)
  names(results) <- c('observed', 'null', 'pval')
  return(results)
}
