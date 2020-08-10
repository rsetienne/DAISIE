#' Single-branch phylogenetic tree
#'
#' @param tip.label Name of the unique tip
#' @param edge.length Length of the unique edge
#'
#' @return A phylogenetic tree in `phylo` format.
#'
#' @examples
#'
#' # A tree with one species with a 15Myr-old stem
#' DAISIE_single_branch("some species", edge.length = 15)
#'
#' @author Raphael Scherrer (github.com/rscherrer)
#'
#' @keywords internal
# Function to display a single branch
DAISIE_single_branch <- function(tip.label = "t1", edge.length = 1) {
  tree <- list(
    edge = matrix(c(2, 1), 1, 2),
    tip.label = tip.label,
    edge.length = edge.length,
    Nnode = 1
  )
  class(tree) <- "phylo"
  tree
}
