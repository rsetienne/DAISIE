#' Convert a simulated archipelago into a datalist for the likelihood
#'
#' @description Converts one replicate of \code{\link{DAISIE_A_sim}()} into
#'   the archipelago datalist format of \code{\link{DAISIE_A_loglik}()}.
#'   Optionally, endemic species are removed from the phylogenies at random
#'   and recorded as missing species (with their ranges).
#'
#' @param replicate One replicate of the output of
#'   \code{\link{DAISIE_A_sim}()}.
#' @param prop_missing Probability that an endemic species is missing from
#'   the phylogeny (default 0). At least one species of every clade is kept.
#' @param missing_range_known Logical; if \code{FALSE}, the ranges of the
#'   missing species are recorded as unknown (\code{NA}).
#'
#' @return Archipelago datalist, see \code{\link{DAISIE_A_loglik}()}.
#' @family DAISIE_A
#' @export
DAISIE_A_sim_to_datalist <- function(replicate,
                                     prop_missing = 0,
                                     missing_range_known = TRUE) {
  first <- replicate[[1]]
  datalist <- list(list(island_age = first$island_age,
                        not_present = first$not_present,
                        n_islands = first$n_islands))
  for (k in seq_along(replicate)[-1]) {
    x <- replicate[[k]]
    missing <- character(0)
    clades <- list()
    for (cl in x$clades) {
      nd <- cl$nodes
      tip_rows <- which(!is.na(nd$species))
      tip_match <- match(nd$species[tip_rows], cl$tips$species)
      nodes <- data.frame(
        node = nd$node,
        parent = nd$parent,
        age = pmax(0, nd$age),
        range = NA_character_,
        type = NA_character_,
        species = nd$species,
        stringsAsFactors = FALSE
      )
      nodes$range[tip_rows] <- cl$tips$range[tip_match]
      nodes$type[tip_rows] <- ifelse(cl$tips$type[tip_match] == "I", "I", "E")
      nodes$age[tip_rows] <- 0
      if (prop_missing > 0) {
        endemic <- tip_rows[nodes$type[tip_rows] == "E"]
        drop <- endemic[stats::runif(length(endemic)) < prop_missing]
        if (length(drop) == length(tip_rows)) drop <- drop[-1]
        if (length(drop) > 0) {
          rng <- nodes$range[drop]
          missing <- c(missing, if (missing_range_known) rng else
            rep(NA_character_, length(rng)))
          nodes <- DAISIE_A_prune_nodes(nodes, setdiff(tip_rows, drop))
        }
      }
      clades[[length(clades) + 1]] <- list(
        colonisation_time = cl$colonisation_time,
        nodes = nodes
      )
    }
    datalist[[k]] <- list(clades = clades, nonendemic_range = NA,
                          missing = missing)
  }
  return(datalist)
}

#' Prune tips from a node table
#'
#' @param nodes Node table (columns \code{node}, \code{parent}, \code{age}
#'   and tip information), with parents before children.
#' @param keep Row numbers of the tips to keep.
#' @return The node table of the reconstructed tree of the kept tips,
#'   renumbered such that parents precede children.
#' @keywords internal
DAISIE_A_prune_nodes <- function(nodes, keep) {
  nodes <- nodes[order(nodes$node), , drop = FALSE]
  nn <- nrow(nodes)
  pos <- stats::setNames(seq_len(nn), nodes$node)
  parent_row <- ifelse(nodes$parent == 0, 0L, pos[as.character(nodes$parent)])
  marked <- rep(FALSE, nn)
  for (r in keep) {
    w <- r
    while (w > 0 && !marked[w]) {
      marked[w] <- TRUE
      w <- parent_row[w]
    }
  }
  n_marked_children <- tabulate(parent_row[marked & parent_row > 0],
                                nbins = nn)
  is_tip <- !(seq_len(nn) %in% parent_row)
  kept <- (is_tip & seq_len(nn) %in% keep) | (n_marked_children == 2)
  kept_rows <- which(kept)
  new_parent <- vapply(kept_rows, function(r) {
    w <- parent_row[r]
    while (w > 0 && !kept[w]) w <- parent_row[w]
    as.integer(w)
  }, integer(1))
  out <- nodes[kept_rows, , drop = FALSE]
  out$parent <- ifelse(new_parent == 0, 0L, match(new_parent, kept_rows))
  out$node <- seq_along(kept_rows)
  rownames(out) <- NULL
  return(out)
}

#' Create a clade for the archipelago likelihood from a phylogeny
#'
#' @description Converts a dated, ultrametric phylogeny (class \code{phylo}
#'   of package ape) of the species of one colonisation event, with tip
#'   ranges and endemicity, into the clade format of
#'   \code{\link{DAISIE_A_loglik}()}. Requires the package ape.
#'
#' @param phy A rooted, bifurcating, ultrametric phylogeny with extant tips.
#' @param tip_ranges Named character vector (names are tip labels) with the
#'   range of every tip, e.g. \code{c(sp1 = "110", sp2 = "001")}.
#' @param tip_types Named character vector with \code{"E"} (endemic) or
#'   \code{"I"} (non-endemic) for every tip; default all endemic.
#' @param colonisation_time Colonisation time (age) of the clade.
#'
#' @return A clade: list with \code{colonisation_time} and \code{nodes}.
#' @family DAISIE_A
#' @export
DAISIE_A_phylo_to_clade <- function(phy,
                                    tip_ranges,
                                    tip_types = NULL,
                                    colonisation_time) {
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is needed for DAISIE_A_phylo_to_clade().")
  }
  ntip <- length(phy$tip.label)
  if (is.null(tip_types)) {
    tip_types <- stats::setNames(rep("E", ntip), phy$tip.label)
  }
  depth <- ape::node.depth.edgelength(phy)
  age <- max(depth) - depth
  root <- ntip + 1
  # breadth-first order: parents before children
  order_nodes <- root
  k <- 1
  while (k <= length(order_nodes)) {
    ch <- phy$edge[phy$edge[, 1] == order_nodes[k], 2]
    order_nodes <- c(order_nodes, ch)
    k <- k + 1
  }
  parent_of <- rep(0L, max(phy$edge))
  parent_of[phy$edge[, 2]] <- phy$edge[, 1]
  new_id <- match(seq_len(max(phy$edge)), order_nodes)
  nodes <- data.frame(
    node = seq_along(order_nodes),
    parent = ifelse(parent_of[order_nodes] == 0, 0L,
                    new_id[parent_of[order_nodes]]),
    age = pmax(0, age[order_nodes]),
    range = NA_character_,
    type = NA_character_,
    label = NA_character_,
    stringsAsFactors = FALSE
  )
  tips <- order_nodes <= ntip
  labels <- phy$tip.label[order_nodes[tips]]
  nodes$label[tips] <- labels
  nodes$range[tips] <- tip_ranges[labels]
  nodes$type[tips] <- tip_types[labels]
  nodes$age[tips] <- 0
  if (colonisation_time < max(nodes$age)) {
    stop("The colonisation time must be older than the crown age.")
  }
  list(colonisation_time = colonisation_time, nodes = nodes)
}
