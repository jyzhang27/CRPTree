#' Generate a tree according to the CRP-TREE model
#'
#' This function generates a ranked, planar, partially labeled tree under
#' the CRP-TREE(alpha) model with N tips and B labeled blue
#'
#' @param N Number of tips
#' @param B Number of tips labeled blue
#' @param alpha parameter of the CRP-TREE model
#' @return A \code{phylo} object generated from CRP-TREE(alpha)
#'
#' @import ape
#' @import apTreeshape
#' @import phangorn
#' @import castor
#' @import TreeTools
#' @import phytools
#' @import parallel
#' @import arrangements
#' @import stepR
#' @import stats
#' @import utils
#'
#' @export
#' @examples rcrp_tree(10, 5, 2)
rcrp_tree <- function(N, B, alpha) {
  # Exceptions
  if (N < 0 | N %% 1 != 0) {
    stop('Invalid number of tips')
  }

  if (B > N | B < 0 | B %% 1 != 0) {
    stop('Invalid number of label type')
  }

  if (alpha < 1) {
    stop('Invalid alpha parameter')
  }

  R <- N-B
  tips <- c(rep(-1, B), rep(-2,R))
  # Create the random ordering:
  # B = -1, R = -2
  rtips <- sample(tips)

  # 2-tip tree
  if (N == 2) {
    edge_list = matrix(ncol=2, nrow=2)
    edge_list[1,] = c(N+1, 2)
    edge_list[2,] = c(N+1, 1)
    tree <- list(edge = edge_list, tip.label = c(rtips[2], rtips[1]), Nnode = N-1, edge.length=c(1,1))
    class(tree) <- 'phylo'
    return(tree)
  }

  # Calculate w_k
  cumulative_blue <- cumsum(rtips==-1)
  cumulative_red <- cumsum(rtips==-2)

  # New blue attached to either blue or red
  prob_blue <- alpha*cumulative_blue/(cumulative_blue*alpha+ cumulative_red)
  #prob_blue_red <- 1/(cumulative_red*alpha+ cumulative_blue)

  # New red attached to either blue or red
  #prob_red_blue <- 1/(cumulative_blue*alpha+ cumulative_red)
  prob_red<- alpha*cumulative_red/(cumulative_red*alpha+ cumulative_blue)

  edge_list = matrix(ncol=2, nrow=2*N-2)
  edge_list[1,] = c(N+1, 2)
  edge_list[2,] = c(N+1, 1)
  curr_root = N+1

  for (k in 3:N) {
    # k is the tip label to be added
    # so currently there are k-1 tips
    curr_tips = rtips[1:(k-1)]
    tip = rtips[k]

    # put the probabilities in order of root, same, different
    # if (tip == -1) {
    #    probs[k,] = prob_blue[k-1,]* c(cumulative_blue[k-1], cumulative_red[k-1])
    # } else{
    #    probs[k,] = prob_red[k-1,]* c(cumulative_red[k-1], cumulative_blue[k-1])
    # }

    # choose the type: 1=same, 2=different

    if (tip == -1) {
      prob <- prob_blue[k-1]
    } else {
      prob <- prob_red[k-1]
    }

    select_type <- rbinom(1,1, prob)

    if (select_type == 1) {
      # select node of same color
      selected_node <- sample(which(curr_tips == tip), size=1)

      # first edge we have to change the old edge that connected old parent and selected node
      index = which(edge_list[, 2] == selected_node)
      old_parent = edge_list[index, 1]
      edge_list[index, ] = c(old_parent,k+N-1)
      # second edge is new internal and new tip
      edge_list[2*k-3,] = c(k+N-1, k)
      # third edge is new internal and the one we are adding to
      edge_list[2*k-2,] = c(k+N-1, selected_node)
    } else {
      # select node of diff color
      selected_node <- sample(which(curr_tips != tip), size=1)

      # first edge we have to change the old edge that connected old parent and selected node
      index = which(edge_list[, 2] == selected_node)
      old_parent = edge_list[index, 1]
      edge_list[index, ] = c(old_parent,k+N-1)
      # second edge is new internal and new tip
      edge_list[2*k-3,] = c(k+N-1, k)
      # third edge is new internal and the one we are adding to
      edge_list[2*k-2,] = c(k+N-1, selected_node)
    }
  }

  # Create branch length 1 between each internal node
  temp_edge<-edge_list[,2]
  temp_edge[temp_edge<(N+1)]<-2*N
  branch_lengths<-abs(temp_edge-edge_list[,1])

  # Create desired format
  tree <- list(edge = edge_list, tip.label = 1:N, Nnode = N-1, edge.length=branch_lengths)
  class(tree) <- "phylo"
  tree <- read.tree(text=NewickTree(tree))
  tree$tip.label = rtips[strtoi(tree$tip.label)]

  return(tree)
}

#' Plot a tree with two categories of labelings
#'
#' This function plots a partially labeled tree with tip colors red and blue
#'
#' @param tree A \code{phylo} object to plot
#' @param node_labels Boolean on whether or not to highlight internal node labels
#' @param cex_tip Numeric for the scaling factor of the tip labels
#' @param cex_node Numeric for the scaling factor of the internal node labels
#' @param frame String specifying the shape of the internal node labels
#' @param bg String specifying the background color of internal node labels
#' @param align_tips Boolean specifying whether tips should be same length from root
#' @return NULL
#' @export
plot_crp_tree <- function(tree, node_labels=FALSE, cex_tip=1.5, cex_node= 0.5, frame = "circle", bg='darkseagreen1', align_tips=TRUE) {
  tip_labels <- tree$tip.label
  tree$tip.label <- ifelse(tip_labels == -1, '+', '*')
  tip_colors <- ifelse(tip_labels == -1, 'blue', 'red')
  plot(tree, direction='downwards', tip.color=tip_colors, cex=cex_tip,  align.tip.label=align_tips)
  if (node_labels) {
    nodelabels(frame =frame, bg=bg, cex=cex_node)
  }
}

#' Calculates the probability of a ranked, planar, partially labeled tree shape
#'
#' This function computes the probability of the ranked, planar, partially labeled
#' tree shape under the CRP-TREE model for a specific alpha value
#'
#' @param tree A \code{phylo} object: a ranked, planar, partially labeled tree shape
#' @param alpha Numeric value for alpha
#' @param only_prob Boolean specifying whether to return only the probability (default \code{TRUE}). If \code{FALSE}, returns a list with the log probability and sufficient statistics S and Ws
#' @return Numeric of the probability or list of log probability and sufficient statistics
#' @export
pcrp_tree <- function(tree, alpha, only_prob = TRUE) {
  tree_data <- get_tree_sufficient_stats(tree)
  S <- tree_data[1]
  W <- tree_data[2:length(tree_data)]
  B <- sum(tree$tip.labels == -1)
  N <- tree$Nnode + 1
  K <- 3:N
  product_terms <- K-1-W + alpha * W

  prob <- alpha^S /choose(N,B)
  #prob <- prob/prod(product_terms)
  log_prob <- log(prob) - sum(log(product_terms))
  if (!only_prob) {
    return(list(log_prob, S, W))
  } else {
    return(exp(log_prob))
  }
}

#' This function computes the conditional probability of partially labeled tree shape
# given the ranked tree shape under CRP-TREE(alpha=1)
#'
#' @param tree A \code{phylo} object: a ranked, partially labeled tree shape
#' @return Numeric of the probability
#' @export
pnull_tree <- function(tree) {
  num_same_cherries <- length(get_internal_nodes_cherry_same(tree))
  edge_matrix <- tree$edge
  num_cherries <- length(which(tabulate(edge_matrix[, 1][edge_matrix[, 2] <= (tree$Nnode+1)]) == 2 ))
  B <- sum(tree$tip.label == -1)
  prob_null <- (2^(num_cherries -num_same_cherries))/choose(tree$Nnode +1, B)
  return(prob_null)
}


#' This function computes the number of same attachments for the given tree
#'
#' @param tree A \code{phylo} object: a ranked, planar, partially labeled tree shape
#' @param nodes_to_switch vector of internal nodes to switch the left and right subtrees when calculating S (default empty)
#' @param tip_order Boolean for whether to return the order the tips were added or not
#' @return Integer: number of same attachments
#' @export
count_same_attachments <- function(tree, nodes_to_switch=c(), tip_order=FALSE) {
  N <- tree$Nnode +1
  internal_nodes <- (N+1):(2*N-1)

  if (length(nodes_to_switch) > 0) {
    if (min(nodes_to_switch) < N+1 | max(nodes_to_switch) > 2*N-1 | length(nodes_to_switch) > N-1) {
      stop('Invalid list of nodes to switch')
    }
    nodes_to_switch <- as.integer(internal_nodes %in% nodes_to_switch)
  } else {
    nodes_to_switch <- rep(0, N-1)
  }

  tree_colors <- tree$tip.label
  all_distances_pairs <- get_all_pairwise_distances(tree)
  distance_to_root <- all_distances_pairs[(N+1), internal_nodes]

  coalescent_order <- internal_nodes[order(distance_to_root, decreasing=TRUE)]
  edges <- tree$edge
  ordered_edges_by_internal_node <- edges[order(distance_to_root, decreasing=TRUE)]
  tip_order_added <- vector('integer', N)
  tip_attached_order <- vector('integer', N-2)

  count_same <- 0
  for (i in 1:(N-2)) {

    node <- coalescent_order[i]
    edges_with_node_as_parent <- which(edges[, 1] == node)
    children <- edges[edges_with_node_as_parent, 2]

    left_child <- edges[min(edges_with_node_as_parent), 2]
    right_child <- edges[max(edges_with_node_as_parent), 2]
    #print(paste('left', left_child,'right',right_child))

    num_children_are_tips <- (left_child <= N) + (right_child<= N)

    if (num_children_are_tips ==  2) {
      # its a cherry
      tip_to_add= left_child # the left tip was added
      tip_added_to = right_child
    } else if (num_children_are_tips == 1) {
      # one is tip and one is internal
      if (left_child <= N) {
        # left is tip, this is being added
        tip_to_add = left_child

        dist_to_node <- all_distances_pairs[, node]
        dist_difference <- dist_to_node - dist_to_node[left_child]
        same_dist <- which(abs(dist_difference) <= 1.5e-08)
        same_dist <- same_dist[(same_dist > left_child) & (same_dist <= N)] #only nodes to right and only tips
        tip_added_to <- setdiff(same_dist, tip_order_added)
      } else {
        # right is tip, it is being added to
        tip_added_to <- right_child

        dist_to_node <- all_distances_pairs[, node]
        dist_difference <- dist_to_node - dist_to_node[right_child]
        same_dist <- which(abs(dist_difference) <= 1.5e-08)
        same_dist <- same_dist[(same_dist < right_child) & (same_dist <=N)] #only nodes to left and only tips
        tip_to_add <- setdiff(same_dist, tip_order_added)
      }

    } else {
      # both are internal
      # tip to add in left child descendants
      # tip added to in right child descendants

      dist_to_left_child <- all_distances_pairs[1:N, left_child]
      dist_to_right_child <- all_distances_pairs[1:N, right_child]

      dist_difference_left <- dist_to_left_child - min(dist_to_left_child)
      dist_difference_right <- dist_to_right_child - min(dist_to_right_child)

      possible_tip_to_add <- which(abs(dist_difference_left) <= 1.5e-08)
      possible_tip_added_to <- which(abs(dist_difference_right) <= 1.5e-08)

      tip_to_add <- setdiff(possible_tip_to_add, tip_order_added)
      tip_added_to <- setdiff(possible_tip_added_to, tip_order_added)
    }

    # if switch, swap the L/R
    if (nodes_to_switch[node-N] == 1) {
      temp <- tip_added_to
      tip_added_to <- tip_to_add
      tip_to_add <- temp
    }

    # now we have to store the tip to add
    tip_order_added[N-i+1] <- tip_to_add

    tip_attached_order[i] <- tip_added_to

    if (tree_colors[tip_to_add] == tree_colors[tip_added_to]) {
      count_same = count_same + 1
    }
  }

  if (tip_order) {
    # write the first two added
    remaining <- setdiff(1:N, tip_order_added)
    if (nodes_to_switch[1] == 1) {
      tip_order_added[1] = min(remaining)
      tip_order_added[2] = max(remaining)
    } else {
      tip_order_added[1] = max(remaining)
      tip_order_added[2] = min(remaining)
    }

    tip_order_added<- order(tip_order_added)
    tip_attached_final <- vector('integer', N)
    for (j in 1:N) {
      if (tip_order_added[j] > 2) {
        tip_attached_final[j] <- tip_order_added[tip_attached_order[N-tip_order_added[j]+1]]
      } else {
        tip_attached_final[j] <- 0
      }
    }

    result <- list(tip_order_added, tip_attached_final, count_same)
    names(result) <- c('tip_order', 'tip_attached_order', 'num_same')
    # tip_order: from L-R on tree, the order in which tips were added
    # tip_attached_order:
    return(result)
  } else {
    return(count_same)
  }
}

#############################################################
# Helpers

# Get internal nodes that are cherries and all the children of cherries
get_cherries_and_children <- function(tree) {
  edge_matrix <- tree$edge
  n <- tree$Nnode + 1
  all_cherries <- which(tabulate(edge_matrix[, 1][edge_matrix[, 2] <= n]) == 2 )
  children <- sapply(all_cherries, function(x) { edge_matrix[which(edge_matrix[,1] == x), 2]} )
  result <- list(all_cherries, children)
  names(result) <- c('cherries', 'children')
  return(result)
}

# Get internal nodes that are cherries with same label
get_internal_nodes_cherry_same <- function(tree) {
  cherries_children <- get_cherries_and_children(tree)
  children <- cherries_children$children
  cherries <- cherries_children$cherries
  children_colors <- matrix(tree$tip.label[children], ncol=2, byrow=TRUE)
  same_cherries <- cherries[which(children_colors[,1] == children_colors[,2])]
  return(same_cherries)
}

# Create different planar representation of the tree
# by swapping left and right of given nodes
reorder_branches <- function(tree, nodes) {
  N <- tree$Nnode+1
  edges <- tree$edge
  branch_length <-tree$edge.length
  tip_labels <- tree$tip.label
  edges_w_length <- cbind(edges, branch_length)

  for (node in nodes) {
    #print(edges_w_length)
    edge_node <- which(edges_w_length [,1] == node)
    first_edge <- min(edge_node)
    second_edge <- max(edge_node)
    #print(c(first_edge, second_edge))
    temp_1 <- edges_w_length[first_edge:(second_edge-1), ]

    second_edge_child <- edges_w_length[second_edge,2]
    if (second_edge_child <= N) {
      temp_2 <- edges_w_length[second_edge, ]
    } else{
      second_edge_all_children <- getDescendants(tree, second_edge_child)
      second_edge_all_children_tips <- second_edge_all_children[second_edge_all_children <= N]
      second_edge_last_edge <- max(which(edges_w_length [,2] %in% second_edge_all_children))
      #print(second_edge_last_edge)
      temp_2 <- edges_w_length[second_edge:second_edge_last_edge, ]
    }
    dim_1 <- if (is.null(nrow(temp_1))) 1 else nrow(temp_1)
    dim_2 <- if (is.null(nrow(temp_2))) 1 else nrow(temp_2)

    edges_new_w_length <- edges_w_length
    edges_new_w_length[first_edge:(first_edge+dim_2-1), ] <- temp_2
    edges_new_w_length[(first_edge+dim_2):(second_edge+dim_2-1),] <- temp_1

    edges_w_length <- edges_new_w_length
    #print(edges_w_length)
  }

  edges_new <- edges_w_length[,1:2]

  new_tiporder <- edges_new[which(edges_new[,2] <= N), 2]
  edges_new[which(edges_new[,2] <= N), 2] <- 1:N

  edges_renumbered<- RenumberTree(parent= edges_new[,1], child=edges_new[,2])

  tree_new <- tree
  tree_new$edge <- edges_renumbered
  tree_new$edge.length <- edges_w_length[,3]
  tree_new$tip.label <- tip_labels[new_tiporder]
  #plot(tree_new, direction='downwards')

  return(tree_new)
}

# Computes the sufficient statistics of the tree
# Returns a vector with first entry S, and other entries W_k
get_tree_sufficient_stats <- function(tree) {
  N <- tree$Nnode + 1
  data <- count_same_attachments(tree, tip_order=TRUE)
  S <- data$num_same
  C <- tree$tip.label[order(data$tip_order)]
  W <- c()

  for (k in 3:N) {
    W[k-2] <- sum(C[1:(k-1)] == C[k])
  }
  return(c(S, W))
}

# Output the attachment matrix based on a tree and planar representation
compute_attachment_matrix <- function(tree, nodes_to_switch) {
  N <- tree$Nnode +1
  if (length(nodes_to_switch) != N-1) {
    stop('Invalid list of nodes to switch: too many or two few')
  }

  if (length(setdiff(unique(nodes_to_switch), c(0,1))) != 0) {
    stop('Invalid list of nodes to switch: some nodes not specified')
  }

  internal_nodes <- (N+1):(2*N-1)

  all_distances_pairs <- get_all_pairwise_distances(tree)
  distance_to_root <- all_distances_pairs[(N+1), internal_nodes]

  coalescent_order <- internal_nodes[order(distance_to_root, decreasing=TRUE)]
  edges <- tree$edge
  ordered_edges_by_internal_node <- edges[order(distance_to_root, decreasing=TRUE)]
  tip_order_added <- vector('integer', N)
  tip_order_added_to <- rep(0, N)


  for (i in 1:(N-2)) {
    node <- coalescent_order[i]
    edges_with_node_as_parent <- which(edges[, 1] == node)
    children <- edges[edges_with_node_as_parent, 2]

    left_child <- edges[min(edges_with_node_as_parent), 2]
    right_child <- edges[max(edges_with_node_as_parent), 2]

    num_children_are_tips <- (left_child <= N) + (right_child<= N)

    if (num_children_are_tips ==  2) {
      # its a cherry
      tip_to_add= left_child # the left tip was added
      tip_added_to = right_child
    } else if (num_children_are_tips == 1) {
      # one is tip and one is internal
      if (left_child <= N) {
        # left is tip, this is being added
        tip_to_add = left_child

        dist_to_node <- all_distances_pairs[, node]
        dist_difference <- dist_to_node - dist_to_node[left_child]
        same_dist <- which(abs(dist_difference) <= 1.5e-08)
        same_dist <- same_dist[(same_dist > left_child) & (same_dist <= N)] #only nodes to right and only tips

        tip_added_to <- setdiff(same_dist, tip_order_added)
      } else {
        # right is tip, it is being added to
        tip_added_to <- right_child

        dist_to_node <- all_distances_pairs[, node]
        dist_difference <- dist_to_node - dist_to_node[right_child]
        same_dist <- which(abs(dist_difference) <= 1.5e-08)
        same_dist <- same_dist[(same_dist < right_child) & (same_dist <=N)] #only nodes to left and only tips
        tip_to_add <- setdiff(same_dist, tip_order_added)
      }

    } else {
      # both are internal
      # tip to add in left child descendants
      # tip added to in right child descendants

      dist_to_left_child <- all_distances_pairs[1:N, left_child]
      dist_to_right_child <- all_distances_pairs[1:N, right_child]

      dist_difference_left <- dist_to_left_child - min(dist_to_left_child)
      dist_difference_right <- dist_to_right_child - min(dist_to_right_child)

      possible_tip_to_add <- which(abs(dist_difference_left) <= 1.5e-08)
      possible_tip_added_to <- which(abs(dist_difference_right) <= 1.5e-08)

      tip_to_add <- setdiff(possible_tip_to_add, tip_order_added)
      tip_added_to <- setdiff(possible_tip_added_to, tip_order_added)
    }

    # if switch, swap the L/R
    if (nodes_to_switch[node-N] == 1) {
      temp <- tip_added_to
      tip_added_to <- tip_to_add
      tip_to_add <- temp
    }

    # now we have to store the tip to add
    tip_order_added[N-i+1] <- tip_to_add
    tip_order_added_to[N-i+1] <- tip_added_to
  }

  # write the first two added
  remaining <- setdiff(1:N, tip_order_added)
  if (nodes_to_switch[1] == 1) {
    tip_order_added[1] = min(remaining)
    tip_order_added[2] = max(remaining)
  } else {
    tip_order_added[1] = max(remaining)
    tip_order_added[2] = min(remaining)
  }

  total <- cbind(tip_order_added, tip_order_added_to)
  return(total)
}

# Compute S using a matrix and node colors
compute_S_given_attachment_matrix <- function(matrix, tree_colors) {
  N <- length(tree_colors)
  tree_colors_order_added <- tree_colors[matrix[3:N, 1]]
  tree_colors_order_added_to <- tree_colors[matrix[3:N, 2]]
  S <- sum(tree_colors_order_added == tree_colors_order_added_to)
  return(S)
}

