# Functions to convert back and forth from 3 representations
# tree -> tables -> attachments -> tree

#' Find the "tables" associated with a CRP-TREE tree
#'
#' This function computes the table representation associated with a CRP-TREE tree
#'
#' @param tree A \code{phylo} object: a ranked, planar, partially labeled tree shape
#' @return List of tables where the numbers corresponds to the tips on the tree with rightmost = 1, leftmost = N
#' @export
get_crp_tables <- function(tree) {
  N <- tree$Nnode + 1
  data <- count_same_attachments(tree, tip_order=TRUE)
  colors <- tree$tip.label[order(data$tip_order)]
  attachments_matrix <- cbind(data$tip_order, data$tip_attached_order)
  attachments_matrix <- attachments_matrix[attachments_matrix[,2] != 0, ]
  colors_matrix <- cbind(colors[attachments_matrix[,1]], colors[attachments_matrix[,2]])
  same <- colors_matrix[,1] == colors_matrix[,2]
  attachments_df <- as.data.frame(attachments_matrix)
  names(attachments_df) <- c('tip_added', 'tip_attached_to')
  attachments_df$same <- same
  attachments_df <- attachments_df[order(attachments_df$tip_added), ]
  rownames(attachments_df) <- NULL

  true_positions <- matrix(nrow=N, ncol=2)
  true_positions[,1] <- 1:N

  if (colors[1] != colors[2]) {
    tables <- list({1}, {2})
    true_positions[1,2] <- 1
    true_positions[2,2] <- 2
  } else {
    tables <- list({c(2,1)})
    true_positions[1,2] <- 1
    true_positions[2,2] <- 1
  }

  for (i in 1:(N-2)) {
    if (attachments_df$same[i]) {
      table_ind <- true_positions[attachments_df$tip_attached_to[i],2]

      select_table <- tables[[table_ind]]
      seat_ind <- which(select_table == attachments_df$tip_attached_to[i])
      select_table <- append(select_table, attachments_df$tip_added[i], after=seat_ind-1)
      tables[[table_ind]] <- select_table
      true_positions[(i+2), 2] <- table_ind
    } else {
      tables <- append(tables, list({c(attachments_df$tip_added[i], attachments_df$tip_attached_to[i])}))
      true_positions[(i+2), 2] <- length(tables)
    }
  }

  return(tables)
}


#' Find the attachments and color type for a list of tables
#'
#' This function computes the attachments and type of attachment (same color or not) for a list of tables satisfying rules to correspond to a CRP-TREE tree
#'
#' @param tables_list A list of vectors indicating the list of tables
#' @return A matrix with three columns: attaching (int), attached_to (int), and same (boolean). Again: rightmost tip = 1, leftmost tip = N
#' @export
crp_tables_to_colors <- function(tables_list) {
  attachments_all <- crp_tables_to_attachments(tables_list)
  total_tables <- length(tables_list)

  start_same <- tail(duplicated(rbind(attachments_all, c(2,1))),1)>0

  if (start_same) {
    print('start same')
    # set to be true for now
    attachments_all$same <- rep(TRUE, nrow(attachments_all))
    all_other_attachments <- attachments_all[which(attachments_all$attached_to != 1 | attachments_all$attaching != 2), ]

    select_start_attachments <- attachments_all[which(attachments_all$attached_to == 1 & attachments_all$attaching == 2), ]
  }

  if (!start_same) {
    print('start different')
    start_attachment_1 <- attachments_all[which(attachments_all$attached_to == 1), ]
    start_attachment_1 <- start_attachment_1[which.min(start_attachment_1$attaching), ]
    start_attachment_2 <- attachments_all[which(attachments_all$attached_to == 2), ]
    start_attachment_2 <- start_attachment_2[which.min(start_attachment_2$attaching), ]


    all_other_attachments <- attachments_all[-as.numeric(c(rownames(start_attachment_1), rownames(start_attachment_2))) , ]
    # set to be true for now
    all_other_attachments$same <- rep(TRUE, nrow(all_other_attachments))

    # smallest attachments
    if (!(1 %in% tables_list)) {
      start_attachment_1$same <- TRUE
    } else {
      start_attachment_1$same <- FALSE
    }

    if (!(2 %in% tables_list)) {
      start_attachment_2$same <- TRUE
    } else {
      start_attachment_2$same <- FALSE
    }

    select_start_attachments <- rbind(start_attachment_1, start_attachment_2)
    select_start_attachments[3,1:2] <- c(2,1)
    select_start_attachments[3,3] <- FALSE

  }

  # Case either where (2,1) does exist, or cases where A_j != 1, 2
  for (i in 1:nrow(all_other_attachments)) {
    att <- as.numeric(all_other_attachments[i, 1:2])
    table_ind <- which(lengths(sapply(tables_list, function(x) which(x %in% att))) == 2)
    table <- tables_list[[table_ind]]

    if (length(table) == 2) {
      # if 2-element list

      all_other_attachments$same[i] <- FALSE
    } else {
      # greater than 2-element list
      # if no smaller element
      if(!any(table[-length(table)] < att[1])) {
        all_other_attachments$same[i] <- FALSE
      }
    }
  }

  all <- rbind(all_other_attachments, select_start_attachments)
  return(all)

}

#' Construct tree from a sequence of attachments and type of attachments
#'
#' This function constructs a tree object based on the given attachment matrix
#'
#' @param attachment_matrix  A matrix with three columns: attaching (int), attached_to (int), and same (boolean)
#' @return A ranked, planar, partially labeled tree
#' @export
create_tree_attachments <- function(attachment_matrix) {
  N <- max(attachment_matrix$attaching)
  attachment_matrix <- attachment_matrix[order(attachment_matrix$attaching), ]

  colors <- vector('integer', N)
  colors[1] <- -1
  for (i in 1:(N-1)) {
    colors[i+1] <- ifelse(attachment_matrix[i, 3], colors[attachment_matrix[i,2]], setdiff(c(-1,-2),colors[attachment_matrix[i,2]]))
  }

  edge_list = matrix(ncol=2, nrow=2*N-2)
  edge_list[1,] = c(N+1, 2)
  edge_list[2,] = c(N+1, 1)
  curr_root = N+1

  for (k in 3:N) {
    # k is the tip label to be added
    # so currently there are k-1 tips
    tip = colors[k]
    selected_node = attachment_matrix[(k-1), 2]

    # first edge we have to change the old edge that connected old parent and selected node
    index = which(edge_list[, 2] == selected_node)
    old_parent = edge_list[index, 1]
    edge_list[index, ] = c(old_parent,k+N-1)
    # second edge is new internal and new tip
    edge_list[2*k-3,] = c(k+N-1, k)
    # third edge is new internal and the one we are adding to
    edge_list[2*k-2,] = c(k+N-1, selected_node)
  }

  # Create branch length 1 between each internal node
  branch_lengths <- c()

  for (i in 1:(2*N-2)) {
    edge <- edge_list[i, ]
    if (all(edge > N)) {
      # both internal nodes
      branch_lengths[i] <- max(edge) - min(edge)
    } else {
      # one internal and one tip
      internal_node <- max(edge)
      branch_lengths[i] <- 2*N - internal_node
    }
  }

  # Create desired format
  tree <- list(edge = edge_list, tip.label = 1:N, Nnode = N-1, edge.length=branch_lengths)
  class(tree) <- "\code{phylo}"
  tree <- read.tree(text=NewickTree(tree))
  node_labels <- colors[strtoi(tree$tip.label)]
  tree$tip.label = node_labels

  return(tree)
}

#####################################################################
# Helpers

# Helper function to get the attachments for a single table
get_attachments_table <- function(table) {
  if (length(table) == 2 & table[1] > table[2]) {
    return(matrix(table,nrow=1, ncol=2))
  }
  which_less <- outer(table, table, '>')
  which_less[lower.tri(which_less)] <- FALSE
  which_less <- which_less[1:(length(table)-1), ]
  which_first <- max.col(which_less, "first")
  if (sum(which_first ==1) > 0) {
    stop('Invalid lists: infeasible ordering')
  }

  attached_to <- table[which_first]
  result <- cbind(table[1:(length(table)-1)], attached_to)
  return(result)
}

# Converts a list of tables to a matrix of attachments
# Again: rightmost tip = 1, leftmost tip = N
crp_tables_to_attachments <- function(tables_list) {
  # check element-related
  all_elements <- unlist(tables_list)
  N <- max(all_elements)
  total_tables <- length(tables_list)
  total_elements <- length(all_elements)
  elements_per_table <- lengths(tables_list)

  if (total_tables > N) {
    stop('Invalid set of lists: too many ')
  } else if (total_elements != N+total_tables -1 & total_elements != N+total_tables -2) {
    stop('Invalid set of lists: incorrect length')
  } else if (length(setdiff(unique(all_elements), 1:N)) != 0) {
    stop('Invalid set of list: missing elements')
  }

  # check per-list related
  one_element_tables <- tables_list[which(elements_per_table == 1)]
  multi_element_tables <- tables_list[which(elements_per_table != 1)]

  if (length(setdiff(unlist(one_element_tables), c(1,2))) != 0) {
    stop('Invalid lists: infeasible start')
  } else if (sum(duplicated(tables_list)) > 0) {
    stop('Invalid lists: repeat elements')
  }

  attachment_list_per_table <- lapply(multi_element_tables, get_attachments_table)
  attachment_list_all <-  as.data.frame(do.call(rbind, attachment_list_per_table))
  names(attachment_list_all)[1] <- 'attaching'

  # Cannot have case (2,1) exists and {1} is table
  if (tail(duplicated(rbind(attachment_list_all, c(2,1))),1)>0 & (1 %in% tables_list)) {
    stop('Invalid lists: invalid starting setup')
  }

  if (length(attachment_list_all$attaching) != N-2 & length(attachment_list_all$attaching) != N-1) {
    stop('Invalid lists: incorrect number of attachments')
  } else if (1 %in% attachment_list_all$attaching) {
    stop('Invalid lists: invalid indices for starting')
  } else if (length(unique(attachment_list_all$attaching)) != nrow(attachment_list_all)) {
    stop('Invalid lists: duplicate attachments')
  } else if (length(setdiff(unique(c(attachment_list_all$attaching, 2)), 2:N)) != 0) {
    stop('Invalid lists: exists invalid attachments')
  }

  return(attachment_list_all)
}

