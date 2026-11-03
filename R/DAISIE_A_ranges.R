# Utilities for ranges (sets of islands) in the archipelago model DAISIE_A.
#
# A range is stored internally as a 0/1 occupancy vector of length n_islands
# (one row of the occupancy matrix of the simulation), or as an integer vector
# of island indices ("set"). Its integer code is sum_i b_i 2^(i - 1), so that
# the 2^n - 1 ranges are numbered 1, ..., 2^n - 1 and island 1 is the least
# significant bit. Its string representation lists b_1 b_2 ... b_n, e.g.
# "110" for a lineage on islands 1 and 2 of a three-island archipelago.

#' Integer code of a range
#'
#' @param occ 0/1 occupancy vector, or 0/1 matrix with one range per row.
#' @return Integer code(s) \eqn{\sum_i b_i 2^{i-1}}.
#' @keywords internal
DAISIE_A_range_code <- function(occ) {
  if (is.matrix(occ)) {
    return(as.numeric(occ %*% (2^(seq_len(ncol(occ)) - 1))))
  }
  sum(occ * 2^(seq_along(occ) - 1))
}

#' String representation of a range
#'
#' @param occ 0/1 occupancy vector, or 0/1 matrix with one range per row.
#' @return Character string(s) such as \code{"101"} (island 1 first).
#' @keywords internal
DAISIE_A_range_string <- function(occ) {
  if (is.matrix(occ)) {
    if (nrow(occ) == 0) return(character(0))
    return(apply(occ, 1, function(x) paste(x, collapse = "")))
  }
  paste(occ, collapse = "")
}

#' Occupancy vector from an integer range code
#'
#' @param code Integer code of the range.
#' @param n_islands Number of islands.
#' @return 0/1 vector of length \code{n_islands}.
#' @keywords internal
DAISIE_A_code_to_occ <- function(code, n_islands) {
  as.numeric(bitwAnd(as.integer(code), as.integer(2^(seq_len(n_islands) - 1))) > 0)
}

#' Is a set of islands connected in the adjacency graph?
#'
#' @param set Integer vector of island indices.
#' @param adjacency Logical adjacency matrix.
#' @return \code{TRUE} if the subgraph induced by \code{set} is connected.
#' @keywords internal
DAISIE_A_is_connected <- function(set, adjacency) {
  if (length(set) <= 1) return(TRUE)
  reached <- set[1]
  repeat {
    nb <- set[colSums(adjacency[reached, set, drop = FALSE]) > 0]
    new_reached <- union(reached, nb)
    if (length(new_reached) == length(reached)) break
    reached <- new_reached
  }
  return(length(reached) == length(set))
}

#' Connected components of a set of islands
#'
#' @inheritParams DAISIE_A_is_connected
#' @return List of integer vectors, one per connected component.
#' @keywords internal
DAISIE_A_components <- function(set, adjacency) {
  comps <- list()
  left <- set
  while (length(left) > 0) {
    reached <- left[1]
    repeat {
      nb <- left[colSums(adjacency[reached, left, drop = FALSE]) > 0]
      new_reached <- union(reached, nb)
      if (length(new_reached) == length(reached)) break
      reached <- new_reached
    }
    comps[[length(comps) + 1]] <- sort(reached)
    left <- setdiff(left, reached)
  }
  return(comps)
}

