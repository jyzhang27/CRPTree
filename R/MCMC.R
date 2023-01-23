
#' Compute the power of the testing methods conditional on a given tree shape
#'
#' This function computes the power of four testing methods (pval_tree, pval_avg, PS, AI) by
#' - sampling planar representations of the given tree using MCMC under specific alpha values
#' - computing the null thresholds
#' Also can return the number of acceptances in each step of the MCMC
#'
#' @param tree A \code{phylo} object: ranked and partially labeled with two categories
#' @param alphas Vector of alphas for which to compute the power under
#' @param sig_level Numeric value of the significance level
#' @param num_sample Integer: number of iterations to run the labeling chain of the MCMC
#' @param num_pl Integer: number of iterations to tun the planar chain of the MCMC
#' @param nperm_null Integer: number of planar permutations of the given tree
#' @param nperm_labels Integer: number of label permutations of the given tree
#' @param acceptances Boolean whether to return the number of acceptances as well (default \code{FALSE})
#'
#' @return A matrix with the power of each method under each alpha (and possible the number of acceptances)
#' @export
crp_tree_power <- function(tree, alphas, sig_level=0.05, num_sample= 500, num_pl=500, nperm_null = 300, nperm_labels=300, acceptances=FALSE) {
  tree_MCMC <- lapply(alphas, MCMC_sample_tree_cond_ranked, fixed_tree= tree, num_sample, num_pl)
  tree_null_thresh <- compute_null_thresholds(tree, nperm_null, nperm_labels, sig_level)
  tree_power_all <- t(sapply(tree_MCMC, compute_power, tree_null_thresh))
  tree_power_all <- cbind(alphas, tree_power_all)

  tree_MCMC_unlist <- unlist(tree_MCMC, recursive=FALSE)
  tree_acceptances_planar <-sapply(tree_MCMC_unlist[seq(5, length(tree_MCMC_unlist), by=6)], mean)
  tree_acceptances_labels <- unlist(tree_MCMC_unlist[seq(6, length(tree_MCMC_unlist), by=6)])
  if (acceptances) {
    result <- list(tree_power_all,tree_acceptances_planar, tree_acceptances_labels)
    names(result) <- c('Power', 'Acceptances_planar', 'Acceptances_labels')
    return(result)
  } else {
    return(tree_power_all)
  }
}

################################################################################
# Helpers


# Function to sample planar representations under CRP-TREE(alpha>1)
# conditional on a ranked partially labeled tree
# Inputs: tree, alpha value, number of samples
# g_func: proposal function, which can be 3 values: unif_one, unif_all, weighted
# Output: Returns probabilities, S values, and number of acceptances in MCMC chain
# Does not return the trees
MCMC_sample_tree_cond_ranked_pl <- function(fixed_tree, alpha, num_sample, g_func='unif_one') {
  N <- fixed_tree$Nnode + 1

  tree_initial <- fixed_tree

  # conditioned on only tree shape, will initialize with randomized labeling
  # so that the chain actually moves

  if (g_func == 'weighted') {
    num_same_cherries <- length(get_internal_nodes_cherry_same(tree_initial))
    num_same_cherries_weight_denom <- 2^(num_same_cherries-1)*(num_same_cherries+2) - num_same_cherries - 1
    num_same_cherries_weight <- (num_same_cherries+1 - (1:num_same_cherries))*choose(num_same_cherries, 1:num_same_cherries) / num_same_cherries_weight_denom
  }

  prob_S_initial <- pcrp_tree(tree_initial, alpha, FALSE)
  pcrp_initial <- prob_S_initial[[1]]
  s_initial <- prob_S_initial[[2]]

  pcrp_all <- vector('double', num_sample+1)
  s_all <- vector('double', num_sample+1)

  pcrp_all[1] <- pcrp_initial
  s_all[1] <- s_initial

  acceptances <- 0

  for (i in 1:num_sample) {
    no_permuting_nodes <- get_internal_nodes_cherry_same(tree_initial)
    permuting_nodes <- setdiff((N+1):(2*N-1), no_permuting_nodes)

    if (g_func== 'unif_all') {
      # switch all with uniform choice
      select_nodes <- rbinom(N-1, 1, 0.5) # choose the nodes to switch
      select_nodes[no_permuting_nodes - N] <- 0 # set the nodes with cherry same to 0
      select_nodes <- which(select_nodes == 1) + N
    } else if (g_func == 'unif_one') {
      # proposal where you only switch one
      select_nodes <- sample(permuting_nodes, 1)
    } else if (g_func == 'weighted') {
      # first select number of nodes to switch
      k <- sample(num_same_cherries, 1, prob = num_same_cherries_weight)
      select_nodes <- sort(sample(permuting_nodes, k))
    } else {
      stop('Invalid proposal function')
    }

    tree_candidate <- reorder_branches(tree_initial, select_nodes)

    prob_S_candidate <- pcrp_tree(tree_candidate, alpha, FALSE)
    pcrp_candidate <- prob_S_candidate[[1]]

    # store log-probabilities for precision reasons
    r <- pcrp_candidate - pcrp_initial
    u <- runif(1)

    # update if meet criteria
    if (log(u) <= r) {
      tree_initial <- tree_candidate
      pcrp_initial <- pcrp_candidate
      s_initial <- prob_S_candidate[[2]]
      acceptances <- acceptances + 1
    }

    pcrp_all[i+1] <- pcrp_initial
    s_all[i+1] <- s_initial
  }
  print(c('acceptances', acceptances))
  return(list(pcrp_all, s_all, acceptances))
}

# Function to sample planar, partially labeled representations uncer CRP-TREE(alpha>1)
# conditional on a ranked tree shape and a given value of $B$
# by simultaneously proposing a planar representation and relabeling
# Inputs: tree, alpha value, number of samples
# Output: Returns probabilities, S values, and number of acceptances in MCMC chain
# Does not return the trees
# Proposal function taken to be unif_one by default
MCMC_sample_tree_simult <- function(fixed_tree, alpha, num_sample) {
  N <- fixed_tree$Nnode + 1

  tree_initial <- fixed_tree
  tree_initial$tip.label <- sample(tree_initial$tip.label)
  pcrp_initial <- pcrp_tree(tree_initial, alpha, FALSE)
  s_initial <- pcrp_initial[[2]]
  pcrp_initial <- pcrp_initial[[1]]

  pcrp_all <- vector('double', num_sample+1)
  trees_all <- vector('list', num_sample+1)
  s_all <- vector('double', num_sample+1)

  trees_all[[1]] <-tree_initial
  pcrp_all[1] <- pcrp_initial[[1]]
  s_all[1] <- s_initial

  acceptances <- 0

  for (i in 1:num_sample) {
    all_cherries <- which(tabulate(tree_initial$edge[, 1][tree_initial$edge[, 2] <= N]) == 2 )
    permuting_nodes <- setdiff((N+1):(2*N-1), all_cherries)
    select_nodes <- sample(permuting_nodes, 1)

    tree_candidate <- reorder_branches(tree_initial, select_nodes)
    tree_candidate$tip.label <- sample(tree_candidate$tip.label)

    prob_S_candidate <- pcrp_tree(tree_candidate, alpha, FALSE)
    pcrp_candidate <- prob_S_candidate[[1]]

    r <- pcrp_candidate - pcrp_initial
    u <- runif(1)
    if (log(u)<= r) {
      tree_initial <- tree_candidate
      pcrp_initial <- pcrp_candidate
      s_initial <- prob_S_candidate[[2]]
      acceptances <- acceptances+1
    }

    trees_all[[i+1]] <- tree_initial
    pcrp_all[i+1] <- pcrp_initial
    s_all[i+1] <- s_initial
  }

  print(c('acceptances', acceptances))
  return(list(trees_all, pcrp_all, s_all, acceptances))
}


# Function to sample planar, partially labeled representations under CRP-TREE(alpha>1)
# conditional on a ranked tree shape and a given value of $B$
# by proposing relabelings first, and then sampling planar representations given these relabeled trees
# Inputs: tree, alpha value, number of samples
# g_func: proposal function, which can be 3 values: unif_one, unif_all, weighted
# Output: Returns probabilities, S values, and number of acceptances in MCMC chain
# Does not return the trees
# Need to run per alpha
MCMC_sample_tree_cond_ranked <- function(fixed_tree, alpha, num_sample, num_pl, g_func='unif_one') {
  N <- fixed_tree$Nnode + 1

  tree_initial <- fixed_tree
  tree_initial$tip.label <- sample(tree_initial$tip.label)
  prob_S_initial <- pcrp_tree(tree_initial, alpha, FALSE)
  pcrp_initial <- prob_S_initial[[1]]
  s_initial <- prob_S_initial[[2]]

  acceptances <- 0
  tree_labeled_sample <- vector('list', 1+num_pl)
  s_labeled_sample <- vector('integer', 1+num_pl)
  pcrp_labeled_sample <- vector('integer', 1+num_pl)

  tree_labeled_sample[[1]] <- tree_initial
  s_labeled_sample[1] <- s_initial
  pcrp_labeled_sample[1] <- pcrp_initial

  for (i in 1:num_pl) {
    # propose new labeling
    tree_candidate <- fixed_tree
    tree_candidate$tip.label <- sample(tree_candidate$tip.label)

    prob_S_candidate <- pcrp_tree(tree_candidate, alpha, FALSE)
    pcrp_candidate <- prob_S_candidate[[1]]

    r <- pcrp_candidate - pcrp_initial
    u <- runif(1)

    # update if meet criteria
    if (log(u) <= r) {
      tree_initial <- tree_candidate
      s_initial <- prob_S_candidate[[2]]
      acceptances <- acceptances + 1
    }

    tree_labeled_sample[[i+1]] <- tree_initial
    s_labeled_sample[i+1] <- s_initial
  }

  MCMC_planar_samples <- mclapply(tree_labeled_sample, MCMC_sample_tree_cond_ranked_pl, alpha,num_sample, g_func, mc.cores=8)
  MCMC_planar_samples <- unlist(MCMC_planar_samples, recursive=FALSE, use.names=FALSE)

  MCMC_pcrp_all <- unlist(MCMC_planar_samples[seq(1, length(MCMC_planar_samples), by=3)], use.names=FALSE)
  MCMC_s_all <- MCMC_planar_samples[seq(2, length(MCMC_planar_samples), by=3)]
  s_mean_all <- unlist(lapply(MCMC_s_all, mean), use.names=FALSE)
  MCMC_s_all <- unlist(MCMC_s_all, use.names=FALSE)
  MCMC_planar_acceptances <- unlist(MCMC_planar_samples[seq(3, length(MCMC_planar_samples), by=3)])

  print(c('acceptances', acceptances))

  result <- list(tree_labeled_sample, MCMC_pcrp_all, MCMC_s_all, s_mean_all, MCMC_planar_acceptances, acceptances)
  names(result) <- c('trees_labeled_MCMC', 'prob', 'S', 'S_mean', 'acceptances_planar', 'acceptances_labels')
  return(result)

}

# Function that computes the values of PS and AI for the sample of trees
# under the alternative
# Input: First element of the MCMC output
other_stats_alt <- function(tree_list_MCMC) {
  tree_MCMC_ai <- sapply(tree_list_MCMC, one_tree_ai, TRUE)
  tree_MCMC_ps <- sapply(tree_list_MCMC, one_tree_fitch, TRUE)
  # tree_MCMC_moran.i <- sapply(tree_list_MCMC, tree_moran.i, FALSE)
  result <- list(tree_MCMC_ai, tree_MCMC_ps)
  names(result) <- c('ai_alt', 'ps_alt')
  #result <- list(tree_MCMC_ai, tree_MCMC_ps, tree_MCMC_moran.i)
  #names(result) <- c('ai_alt', 'ps_alt', 'moran.i_alt')
  return(result)
}

# Function to compute the null thresholds of S, S_mean, PS, AI
# for a given tree shape, number of permutations, and significant levels
# Note S, S_mean: large = more extreme
# PS, AI: small = more extreme
# Returns 4 thresholds
compute_null_thresholds <- function(tree, nperm_null, nperm_labels, sig_level) {
  S_grid_null <- simulate_S_grid(tree, nperm_null, nperm_labels)
  S_mean_null <- colMeans(S_grid_null)
  S_null <- c(S_grid_null)
  S_thresh <- quantile(S_null, 1-sig_level)
  S_mean_thresh <- quantile(S_mean_null, 1-sig_level)

  ai_null <- tree_ai(tree, nperm_labels)$null
  ps_null <- tree_fitch(tree, nperm_labels)$null

  #because smaller is more significant
  ai_thresh <- quantile(ai_null, sig_level)
  ps_thresh <- quantile(ps_null, sig_level)

  thresh <- c(S_thresh, S_mean_thresh, ai_thresh, ps_thresh)
  names(thresh) <- c('S', 'S_mean', 'AI', 'PS')
  return(thresh)
}

# Function to compute the power for the 4 methods
# Need to run once per alpha
# Input: MCMC output, null_threshold output
compute_power <- function(tree_MCMC, null_thresh) {
  S_all_power <- mean(tree_MCMC$S >= null_thresh[1])
  S_mean_power <- mean(tree_MCMC$S_mean >= null_thresh[2])

  other_stats_alt <- other_stats_alt(tree_MCMC$trees_labeled_MCMC)
  ai_power <- mean(other_stats_alt$ai_alt <= null_thresh[3])
  ps_power <- mean(other_stats_alt$ps_alt <= null_thresh[4])

  result <- c(S_all_power, S_mean_power, ai_power, ps_power)
  names(result) <- c('S', 'S_mean', 'AI', 'PS')
  return(result)
}


