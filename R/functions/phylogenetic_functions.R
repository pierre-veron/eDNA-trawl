### ----------------------------------------------------------------------------
#   Defines functions to use phylogenetic trees
#   These functions allow to extract edges from the tree (e.g. to color them)
#   Author: Pierre Veron, pierre.veron.2017@polytechnique.org
### ----------------------------------------------------------------------------

#' get_edges_tip_to_root
#' 
#' Extracts all the edge indices from a tip to the root of the tree. 
#' @param tree a phylogenetic tree (from class phylo)
#' @param tip int, the numer of a tip
#'
#' @return a vector containing all the indices of the edges linking the 
#' specified tip to the root of the tree.
#' The edges indices are the same as in the table in tree$edge
#' @examples
get_edges_tip_to_root <- function(tree, tip) {
  edg <- tree$edge
  sel <- c() 
  cursor <- which(edg[,2] == tip)
  while (length(cursor) > 0) {
    sel <- c(sel, cursor)
    cursor <- which(edg[,2] == edg[cursor, 1])
  }
  sel
}

#' get_edges_parent
#'
#' @param tree a phylogenetic tree (from class phylo)
#' @param tips a vector of tips indices
#'
#' @return a vector containing the sorted and unique edges that are parent from
#' the specified tips.
#' @export
#'
#' @examples
get_edges_parent <- function(tree, tips) {
  sort(unique(unlist(sapply(tips, function(tip) {
    get_edges_tip_to_root(tree, tip)
  }))))
}

#' get_edges_parent_from_label
#'
#' @param tree a phylogenetic tree (from class phylo)
#' @param tip_labels a vector containing some tip labels of tree
#'
#' @return a vector containing the indices of the edges that are parent from the 
#' tips specified by the labels
#' @export
#'
#' @examples
get_edges_parent_from_label <- function(tree, tip_labels) {
  tips <- which(tree$tip.label %in% tip_labels)
  get_edges_parent(tree, tips)
}