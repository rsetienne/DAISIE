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

#' Is a range allowed under the state-space restrictions?
#'
#' @param set Integer vector of island indices.
#' @param pars A \code{DAISIE_A_pars} object.
#' @return \code{TRUE} if the range satisfies \code{max_range_size} and, if
#'   requested, \code{connected_ranges}.
#' @keywords internal
DAISIE_A_range_allowed <- function(set, pars) {
  length(set) >= 1 &&
    length(set) <= pars$max_range_size &&
    (!pars$connected_ranges || DAISIE_A_is_connected(set, pars$adjacency))
}

#' All allowed ranges of an archipelago
#'
#' @param pars A \code{DAISIE_A_pars} object.
#' @return A list with \code{sets} (list of integer vectors), \code{codes}
#'   (integer codes) and \code{strings}.
#' @family DAISIE_A
#' @export
#' @examples
#' pars <- DAISIE_A_create_pars(n_islands = 4, lac = 1, mu = 1, gam = 1,
#'   laa = 1, connectivity = "line", connected_ranges = TRUE)
#' length(DAISIE_A_ranges(pars)$sets) # 10 = 4 * 5 / 2 instead of 15
DAISIE_A_ranges <- function(pars) {
  n <- pars$n_islands
  sets <- list()
  codes <- c()
  for (code in seq_len(2^n - 1)) {
    occ <- DAISIE_A_code_to_occ(code, n)
    set <- which(occ == 1)
    if (DAISIE_A_range_allowed(set, pars)) {
      sets[[length(sets) + 1]] <- set
      codes <- c(codes, code)
    }
  }
  strings <- vapply(sets, function(s) {
    occ <- rep(0, n)
    occ[s] <- 1
    paste(occ, collapse = "")
  }, character(1))
  return(list(sets = sets, codes = codes, strings = strings))
}

#' Admissible vicariance splits of a range
#'
#' @description Enumerates the admissible bipartitions \eqn{\{A, R
#'   \setminus A\}} of a range \eqn{R} according to \code{pars$vic_mode} and
#'   \code{pars$connected_ranges}. Each bipartition is returned once, as the
#'   part \eqn{A} that contains the island with the lowest index.
#'
#' @param set Integer vector of island indices (the range \eqn{R}), sorted.
#' @param pars A \code{DAISIE_A_pars} object.
#' @return List of integer vectors (the parts \eqn{A}).
#' @keywords internal
DAISIE_A_admissible_splits <- function(set, pars) {
  r <- length(set)
  if (r < 2) return(list())
  adj <- pars$adjacency
  if (pars$vic_mode == "single") {
    if (r == 2) {
      splits <- list(set[1])
    } else {
      splits <- lapply(set, function(i) {
        if (i == set[1]) i else setdiff(set, i)
      })
    }
  } else {
    if (r > 20) {
      stop("Enumerating vicariance splits of ranges with more than 20 ",
           "islands is not supported; use vic_mode = 'all' or 'single'.")
    }
    others <- set[-1]
    bits <- as.integer(2^(seq_along(others) - 1))
    splits <- vector("list", 2^(r - 1) - 1)
    k <- 0
    for (mask in 0:(2^(r - 1) - 2)) {
      k <- k + 1
      splits[[k]] <- c(set[1], others[bitwAnd(as.integer(mask), bits) > 0])
    }
    if (pars$vic_mode == "connected") {
      if (DAISIE_A_is_connected(set, adj)) {
        keep <- vapply(splits, function(A) {
          DAISIE_A_is_connected(A, adj) &&
            DAISIE_A_is_connected(setdiff(set, A), adj)
        }, logical(1))
      } else {
        comps <- DAISIE_A_components(set, adj)
        keep <- vapply(splits, function(A) {
          all(vapply(comps, function(cmp) {
            all(cmp %in% A) || !any(cmp %in% A)
          }, logical(1)))
        }, logical(1))
      }
      splits <- splits[keep]
    }
  }
  if (pars$connected_ranges) {
    keep <- vapply(splits, function(A) {
      DAISIE_A_is_connected(A, adj) &&
        DAISIE_A_is_connected(setdiff(set, A), adj)
    }, logical(1))
    splits <- splits[keep]
  }
  return(splits)
}

#' Cached admissible vicariance splits
#'
#' @inheritParams DAISIE_A_admissible_splits
#' @param cache An environment used to store splits per range.
#' @return List of integer vectors, see
#'   \code{\link{DAISIE_A_admissible_splits}()}.
#' @keywords internal
DAISIE_A_cached_splits <- function(set, pars, cache) {
  key <- paste0("s", paste(set, collapse = "_"))
  splits <- cache[[key]]
  if (is.null(splits)) {
    splits <- DAISIE_A_admissible_splits(set, pars)
    assign(key, splits, envir = cache)
  }
  return(splits)
}

#' Number of admissible vicariance splits of a range
#'
#' @description Returns \eqn{s(R)}, the number of admissible bipartitions of
#'   the range \eqn{R} (see \code{\link{DAISIE_A_admissible_splits}()}), using
#'   closed forms where possible: \eqn{2^{|R|-1} - 1} if all bipartitions are
#'   admissible and, for \code{vic_mode = "single"}, \eqn{|R|} (1 if
#'   \eqn{|R| = 2}).
#'
#' @inheritParams DAISIE_A_cached_splits
#' @return Number of admissible splits (0 for a single-island range).
#' @keywords internal
DAISIE_A_n_splits <- function(set, pars, cache = NULL) {
  r <- length(set)
  if (r < 2) return(0)
  if (!pars$connected_ranges) {
    if (pars$vic_mode == "all") return(2^(r - 1) - 1)
    if (pars$vic_mode == "single") return(if (r == 2) 1 else r)
  }
  if (pars$vic_mode != "single") {
    sub <- pars$adjacency[set, set]
    if (all(sub[upper.tri(sub)])) return(2^(r - 1) - 1)
  }
  if (is.null(cache)) return(length(DAISIE_A_admissible_splits(set, pars)))
  return(length(DAISIE_A_cached_splits(set, pars, cache)))
}

#' Weights of the admissible vicariance splits of a range
#'
#' @description The weight of the split of \eqn{R} into \eqn{A} and
#'   \eqn{R \setminus A} is \eqn{w(A) = \exp(-\kappa\, c(A))}, with
#'   \eqn{c(A) = \sum_{i \in A} \sum_{j \in R \setminus A}
#'   (C_{ij} + C_{ji})} the connectivity cut by the split and \eqn{\kappa}
#'   = \code{pars$vic_kappa}. For \eqn{\kappa = 0} all weights are 1.
#'
#' @inheritParams DAISIE_A_n_splits
#' @return Numeric vector of weights, in the order of
#'   \code{\link{DAISIE_A_admissible_splits}()}.
#' @keywords internal
DAISIE_A_split_weights <- function(set, pars, cache = NULL) {
  splits <- if (is.null(cache)) {
    DAISIE_A_admissible_splits(set, pars)
  } else {
    DAISIE_A_cached_splits(set, pars, cache)
  }
  kappa <- DAISIE_A_vic_kappa(pars)
  if (kappa == 0 || length(splits) == 0) return(rep(1, length(splits)))
  key <- paste0("w", paste(set, collapse = "_"))
  if (!is.null(cache) && !is.null(cache[[key]])) return(cache[[key]])
  Cs <- pars$C + t(pars$C)
  w <- vapply(splits, function(A) {
    exp(-kappa * sum(Cs[A, setdiff(set, A), drop = FALSE]))
  }, numeric(1))
  if (!is.null(cache)) assign(key, w, envir = cache)
  return(w)
}

#' Total vicariance weight of a range
#'
#' @description \eqn{W(R) = \sum_A w(A)} over the admissible splits; equal to
#'   the number of admissible splits \eqn{s(R)} if \code{vic_kappa = 0}
#'   (computed in closed form where possible).
#' @inheritParams DAISIE_A_n_splits
#' @return \eqn{W(R)} (0 for a single-island range).
#' @keywords internal
DAISIE_A_vic_weight_total <- function(set, pars, cache = NULL) {
  if (length(set) < 2) return(0)
  if (DAISIE_A_vic_kappa(pars) == 0) {
    return(DAISIE_A_n_splits(set, pars, cache))
  }
  return(sum(DAISIE_A_split_weights(set, pars, cache)))
}

#' Vicariance events of a range and their rates
#'
#' @description For a lineage with range \eqn{R}, every admissible split
#'   \eqn{A} occurs at rate \eqn{\lambda^v w(A)} (\code{vic_per = "split"})
#'   or \eqn{\lambda^v w(A) / W(R)} (\code{vic_per = "lineage"}).
#' @inheritParams DAISIE_A_n_splits
#' @return A list with the admissible \code{splits}, their \code{rates} and
#'   the \code{total} rate.
#' @keywords internal
DAISIE_A_vic_rates <- function(set, pars, cache = NULL) {
  if (length(set) < 2 || pars$lav == 0) {
    return(list(splits = list(), rates = numeric(0), total = 0))
  }
  splits <- if (is.null(cache)) {
    DAISIE_A_admissible_splits(set, pars)
  } else {
    DAISIE_A_cached_splits(set, pars, cache)
  }
  w <- DAISIE_A_split_weights(set, pars, cache)
  rates <- if (DAISIE_A_vic_per(pars) == "split") {
    pars$lav * w
  } else {
    pars$lav * w / sum(w)
  }
  return(list(splits = splits, rates = rates, total = sum(rates)))
}

#' Sample a vicariance split of a range
#'
#' @description The split is drawn with probability proportional to its
#'   weight \eqn{w(A)} (uniformly if \code{vic_kappa = 0}).
#'
#' @inheritParams DAISIE_A_cached_splits
#' @return Integer vector: one daughter range \eqn{A}; the other daughter
#'   range is \code{setdiff(set, A)}.
#' @keywords internal
DAISIE_A_sample_split <- function(set, pars, cache) {
  r <- length(set)
  if (DAISIE_A_vic_kappa(pars) > 0) {
    splits <- DAISIE_A_cached_splits(set, pars, cache)
    w <- DAISIE_A_split_weights(set, pars, cache)
    k <- if (length(splits) == 1) 1L else sample.int(length(splits), 1, prob = w)
    return(splits[[k]])
  }
  if (r > 12 && pars$vic_mode != "single") {
    sub <- pars$adjacency[set, set]
    complete <- all(sub[upper.tri(sub)])
  }
  if (r > 12 && pars$vic_mode != "single" &&
      (complete || (pars$vic_mode == "all" && !pars$connected_ranges))) {
    # All bipartitions are admissible. Rejection sampling: uniform over
    # ordered (A, R \ A), hence uniform over unordered bipartitions.
    repeat {
      pick <- stats::runif(r) < 0.5
      if (any(pick) && !all(pick)) break
    }
    return(set[pick])
  }
  splits <- DAISIE_A_cached_splits(set, pars, cache)
  if (length(splits) == 0) {
    stop("No admissible vicariance split for range ",
         paste(set, collapse = ","), ".")
  }
  k <- if (length(splits) == 1) 1L else sample.int(length(splits), 1)
  return(splits[[k]])
}

#' Which islands can be lost without disconnecting a range?
#'
#' @inheritParams DAISIE_A_cached_splits
#' @return Logical vector of the same length as \code{set}: \code{TRUE} if the
#'   range minus that island is still allowed (always \code{TRUE} unless
#'   \code{pars$connected_ranges}).
#' @keywords internal
DAISIE_A_removable <- function(set, pars, cache) {
  if (!pars$connected_ranges || length(set) <= 2) {
    return(rep(TRUE, length(set)))
  }
  key <- paste0("r", paste(set, collapse = "_"))
  res <- cache[[key]]
  if (is.null(res)) {
    res <- vapply(set, function(i) {
      DAISIE_A_is_connected(setdiff(set, i), pars$adjacency)
    }, logical(1))
    assign(key, res, envir = cache)
  }
  return(res)
}
