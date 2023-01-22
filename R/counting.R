
#' Count the number of ranked partially labeled tree shapes with N tips and B blue
#'
#' @param N Number of tips
#' @param B Number of blue tips
#'
#' @return The number of ranked partially labeled tree shapes with N tips and B blue
#' @export
count_ranked_pl_trees <- function(N,B) {
  if (N<0 | B<0 | B>N) {
    return(0)
  } else if (N==1 | N==2) {
    return(1)
  } else if (N==3) {
    if (B==0 | B==3) {
      return(1)
    } else if (B==1 | B==2) {
      return(2)
    }
  } else if (B==0 | B==N){
    total <- 0
    for (k in 0:(N-2)) {
      total <- total + choose(N-2, k) * count_ranked_pl_trees(k+1, 0)* count_ranked_pl_trees(N-1-k, 0)/2
    }
    return(total)
  } else {
    total <- 0
    for (n in 1:(N-1)) {
      for (b in 0:min(B, n)) {
        total <- total + choose(N-2, n-1) * count_ranked_pl_trees(n,b) * count_ranked_pl_trees(N-n, B-b)/2
      }
    }
    return(total)
  }
}

#' Count the number of partial labelings on a given tree shape
#'
#' This function counts the number of possible partial labelings with N tips and B blue
#' on the given ranked tree shape
#'
#' @param tree A ranked tree shape
#' @param B Number of blue tips
#'
#' @return The number of possible partial labelings
#' @export
count_partial_labels <- function(tree, B) {
  ntips <- tree$Nnode+1

  # set the tip labels to be all unique
  # labels are arbitrary, only need uniqueness
  tree$tip.label <- 1:ntips

  # Base cases for recursion: N= 1, 2, 3,
  if (B < 0 | B > ntips) {
    # case where infeasible number of B
    return(0)
  } else if (ntips <= 2) {
    # tree with 2 leaves
    return(1)
  } else if (B == 0 | B == ntips) {
    # labels are all the same
    return(1)
  } else if (ntips <= 3) {
    # if 3 leaves and have B=1, 2, then its 2
    # B=0, 3 cases already taken care of above
    return(2)
  }

  # for more tips: N >= 4

  # prune the tree at the node
  sub1 <- get_subtree_at_node(tree, 2)$subtree
  n1 <- sub1$Nnode + 1
  n2 <- ntips - n1

  # case where 4 leaves
  if (ntips == 4) {
    if (n1 == 2) {
      if (B == 2) {
        return(3)
      } else {
        return(2)
      }
    } else {
      if (B == 2) {
        return(4)
      } else {
        return(3)
      }
    }
  }

  # if get single node as right-most subtree
  if (n2 == 1) {
    total <- count_partial_labels(sub1, B) + count_partial_labels(sub1, B-1)
    return(total)
  }

  # if two branch right-most subtree
  if (n2 == 2) {
    total <- count_partial_labels(sub1, B) + count_partial_labels(sub1, B-1)
    if (B >= 2) {
      total <- total + count_partial_labels(sub1, B-2)
    }
    return(total)
  }

  # get the other half of the tree
  # above case takes care of 1-tup subtree case
  sub2 <- get_subtree_with_tips(tree, omit_tips = sub1$tip.label)$subtree

  if (n1 == 1) {
    total <- count_partial_labels(sub2, B) + count_partial_labels(sub2, B-1)
    return(total)
  }

  if (n1 == 2) {
    total <- count_partial_labels(sub2, B) + count_partial_labels(sub2, B-1)
    if (B >= 2) {
      total <- total + count_partial_labels(sub2, B-2)
    }
    return(total)
  }

  # Every other case: recursively count based on cases where labellings are feasible
  total <- 0
  for (i in 0:max(n1, n2)) {
    if (n1-i >= 0 & B-i>=0 & n2+i-B >= 0 ) {
      total <- total + count_partial_labels(sub1, i) * count_partial_labels(sub2, B-i)
    }
  }
  return(total)

}

#' Compute expected value of S for a tree generated under the CRP-TREE model
#'
#' This function computes the expected value of S for a tree generated under CRP-TREE(N, B, alpha)
#'
#' @param N Number of tips
#' @param B Number of blue tips
#' @param alpha Numeric parameter for the tree model
#'
#' @return The expected value
#' @export
expected_value_function <- function(N, B, alpha) {
  if (alpha == 1) {
    return(expected_value_null(N, B))
  } else{
    return(expected_value_alt(N, B, alpha))
  }
}

#' Compute the MLE of alpha for a ranked, planar, partially labeled tree
#'
#' This function computes the MLE of alpha for a ranked, planar, partially labeled tree
#' generated under CRP-TREE(N,B, alpha)
#'
#' @param tree A ranked, planar, partially labeled tree
#'
#' @return A vector with the MLE and maximum log-likelihood value
#' @export
get_mle_1tree <- function(tree) {
  N <- tree$Nnode +1
  data <- get_tree_sufficient_stats(tree)
  S <- data[1]
  W <- data[2:length(data)]
  return(optimize(log_like_one, c(1,10000), tol= 1e-6, N=N, S=S, W=W, maximum=TRUE))
}

#' Computes the MLE of alpha for a list of ranked, planar, partially labeled trees
#'
#' This function computes the MLE of alpha for a list of ranked, planar, partially labeled trees
#' with the same N generated under CRP-TREE(N,B, alpha)
#'
#' @param tree_list A list of ranked, planar, partially labeled trees
#'
#' @return A vector with the MLE and maximum log-likelihood value
#' @export
get_mle_trees <- function(tree_list) {
  data_matrix <- sapply(tree_list, get_tree_sufficient_stats)
  Ss <- data_matrix[1, ]
  Ws <- data_matrix[2:nrow(data_matrix), ]
  Ns <- rep(nrow(data_matrix)+1, length(tree_list))

  suff_stats <- rbind(Ns, data_matrix)
  return(optimize(log_like_multiple, c(1,10000), tol= 1e-6, suff_stats, maximum=TRUE))
}

#' Applies LRT for a ranked, planar, partially labeled tree
#'
#' This function applies the LRT for testing alpha = 1 vs alpha > 1 for a
#' ranked, planar, partially labeled tree generated under CRP-TREE(N,B, alpha)
#'
#' @param tree A ranked, planar, partially labeled tree
#'
#' @return A vector with the MLE and p-value
#' @export
lrt_1tree <- function(tree) {
  mle <- get_mle_1tree(tree)
  alpha_mle <- mle$maximum
  alt_like <- mle$objective

  null_like <- -sum(log(2:(tree$Nnode)))

  D <- 2*(alt_like - null_like)
  pval_chi_squared <- 1-pchisq(max(D,0), 1)
  results <- c(alpha_mle, pval_chi_squared)
  names(results) <- c('alpha_MLE', 'pval')
  return(results)
}

#' Applies LRT for a list of ranked, planar, partially labeled trees
#'
#' This function applies the LRT for testing alpha = 1 vs alpha > 1 for a list of
#' ranked, planar, partially labeled trees generated under CRP-TREE(N,B, alpha)
#' with the same N and alpha
#'
#' @param tree A ranked, planar, partially labeled tree
#'
#' @return A vector with the MLE and p-value
#' @export
lrt_many_trees <- function(tree_list) {
  mle <- get_mle_trees(tree_list)
  alpha_mle <- mle$maximum
  alt_like <- mle$objective

  null_like <- -length(tree_list)*sum(log(2:(tree_list[[1]]$Nnode)))

  D <- 2*(alt_like - null_like)
  pval_chi_squared <- 1-pchisq(max(D,0), 1)
  results <- c(alpha_mle, pval_chi_squared)
  names(results) <- c('alpha_MLE', 'pval')
  return(results)
}


################################################################
# Helpers

# Expected value when alpha=1
expected_value_null <- function(N, B) {
  mean <- (N-2)*(B*(B-1) + (N-B) * (N-B-1)) / (N^2-N)
  return(mean)
}


expectation_helper_1 <- function(N, B, alpha, k, i) {
  value1 <- alpha * i/(k-1-i+alpha * i) * B/N * choose(k-1, i) * choose(N-k, B-i-1) / choose(N-1, B-1)
  value2 <- alpha * i/(k-1-i+alpha * i) * (N-B)/N * choose(k-1, i) * choose(N-k, N-B-i-1) / choose(N-1, N-B-1)
  return(value1+value2)
}

expectation_helper_1(500, 200, 2, 250, 100)

expectation_helper_2 <- function(N, B, alpha, k) {
  i <- 0:(k-1)
  terms <- sapply(i, expectation_helper_1, N=N, B=B, alpha=alpha, k=k)
  return(sum(terms))
}


# Expected value when alpha > 1
expected_value_alt <- function(N, B, alpha) {
  k <- 3:N
  terms <- sapply(k, expectation_helper_2, N=N, B=B, alpha=alpha)
  return(sum(terms))
}

log_like_one <- function(alpha, N, S, W) {
  K <- 3:N
  return(S*log(alpha) - sum(log(K-1-W+ alpha*W)))
}

log_like_multiple <- function(alpha, suff_stats) {
  log_like <- apply(suff_stats, 2, function(x) { log_like_one(alpha, x[1], x[2], x[3:length(x)]) })
  return(sum(log_like))
}
