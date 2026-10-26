#' Calculate the event rates of the archipelago model
#'
#' @description Internal function. Computes the total rate of each class of
#'   events of the DAISIE_A model in the current state, together with the
#'   per-candidate weights needed to sample which lineage/island is affected.
#'
#'   With \eqn{n_i} the number of species (of the clade) on island \eqn{i},
#'   \eqn{f_i = \max(0, 1 - n_i/K_i)} the diversity-dependence factor and
#'   \eqn{e_i} an indicator of whether island \eqn{i} has emerged:
#'   \itemize{
#'     \item colonisation of island \eqn{i} by a mainland species:
#'       \eqn{\gamma_i f_i e_i}; for a mainland species that is present as
#'       a non-endemic lineage with range \eqn{R}, only colonisation of
#'       \eqn{i \notin R} (range expansion) or, if \eqn{R = \{i\}},
#'       of \eqn{i} itself (reset of the colonisation time) counts;
#'     \item extinction of the population of a lineage on island \eqn{i}:
#'       \eqn{\mu_i} (global extinction if it is the last island);
#'     \item anagenesis of a non-endemic lineage: \eqn{\lambda^a};
#'     \item in-situ cladogenesis of the population on island \eqn{i}:
#'       \eqn{\lambda^c_i f_i};
#'     \item dispersal of a lineage with range \eqn{R} to island
#'       \eqn{j \notin R}: \eqn{\sum_{i \in R} D_{ij} f_j e_j};
#'     \item vicariance of a lineage with \eqn{|R| \ge 2}: \eqn{\lambda^v}
#'       (\code{vic_per = "lineage"}) or \eqn{\lambda^v s(R)}, with
#'       \eqn{s(R)} the number of admissible splits (\code{vic_per = "split"});
#'     \item island-level anagenesis of the population on island
#'       \eqn{i \in R} of a lineage with \eqn{|R| \ge 2}:
#'       \eqn{\lambda^{ia}_i} (not diversity dependent: the number of
#'       species on every island stays the same).
#'   }
#'   Events that would create a range that is not allowed
#'   (\code{max_range_size}, \code{connected_ranges}) have rate zero.
#'
#' @param st State list of the simulation.
#' @param timeval Current (forward) time.
#' @param pars A \code{DAISIE_A_pars} object.
#' @param emerge Emergence times of the islands (forward time).
#' @param mainland_spec IDs of the mainland species.
#' @param split_cache Environment used for caching range computations.
#'
#' @return A named list with the class totals \code{immig}, \code{ext},
#'   \code{ana}, \code{clado}, \code{disp}, \code{vic}, \code{isl_ana}, their sum
#'   \code{total}, and the weights used for sampling.
#' @keywords internal
DAISIE_A_update_rates <- function(st,
                                  timeval,
                                  pars,
                                  emerge,
                                  mainland_spec,
                                  split_cache) {
  n <- pars$n_islands
  spec <- st$spec
  L <- nrow(spec)
  occ <- spec[, .DA_NFIX + seq_len(n), drop = FALSE]
  exists <- as.numeric(emerge <= timeval)
  n_i <- colSums(occ)
  f <- pmax(0, 1 - n_i / pars$K)
  range_size <- rowSums(occ)
  type <- spec[, .DA_TYPE]

  # Which islands can be added to the range of each lineage
  add_ok <- (1 - occ) * as.numeric(range_size < pars$max_range_size)
  if (pars$connected_ranges && L > 0) {
    add_ok <- add_ok * ((occ %*% pars$adjacency) > 0)
  }

  # Colonisation from the mainland (always to a single island)
  col_isl <- pars$gam * f * exists
  w_col <- rep(sum(col_isl), length(mainland_spec))
  I_rows <- which(type == 1)
  if (length(I_rows) > 0) {
    for (r in I_rows) {
      m <- match(spec[r, .DA_ID], mainland_spec)
      # recolonisation of the island of a single-island lineage (reset of
      # the colonisation time) or of a new island (range expansion without
      # reset), if allowed; recolonisation of an island in a widespread
      # range has no effect and is not counted
      w_col[m] <- sum(col_isl * (DAISIE_A_reset_islands(occ[r, ]) +
                                   add_ok[r, ]))
    }
  }
  immig <- sum(w_col)

  # Local (and global) extinction
  rem_ok <- occ
  if (pars$connected_ranges && L > 0) {
    for (r in which(range_size >= 3)) {
      set <- which(occ[r, ] == 1)
      rem_ok[r, set] <- as.numeric(DAISIE_A_removable(set, pars, split_cache))
    }
  }
  ext <- sum(colSums(rem_ok) * pars$mu)

  # Anagenesis of non-endemic lineages
  ana <- pars$laa * length(I_rows)

  # In-situ cladogenesis (split or budding)
  clado_isl <- pars$lac * f
  clado <- sum(n_i * clado_isl)

  # Dispersal (range expansion)
  Wd <- NULL
  disp <- 0
  if (L > 0 && any(pars$D > 0)) {
    Wd <- (occ %*% pars$D) * add_ok * rep(f * exists, each = L)
    disp <- sum(Wd)
  }

  # Vicariance (between-island speciation): rate lav per lineage, or lav per
  # admissible split of its range
  wide <- which(range_size >= 2)
  vic_w <- rep(1, length(wide))
  if (DAISIE_A_vic_per(pars) == "split" && pars$lav > 0 && length(wide) > 0) {
    vic_w <- vapply(wide, function(r) {
      DAISIE_A_n_splits(which(occ[r, ] == 1), pars, split_cache)
    }, numeric(1))
  }
  vic <- pars$lav * sum(vic_w)

  # Island-level anagenesis: the population on island i of a widespread
  # lineage becomes a new endemic species on i; the parent keeps R \ {i}
  laa_isl <- DAISIE_A_laa_isl(pars)
  Wi <- NULL
  isl_ana <- 0
  if (L > 0 && any(laa_isl > 0)) {
    Wi <- rem_ok * as.numeric(range_size >= 2) * rep(laa_isl, each = L)
    isl_ana <- sum(Wi)
  }

  total <- immig + ext + ana + clado + disp + vic + isl_ana
  rates <- list(
    immig = immig,
    ext = ext,
    ana = ana,
    clado = clado,
    disp = disp,
    vic = vic,
    isl_ana = isl_ana,
    total = total,
    w_col = w_col,
    col_isl = col_isl,
    add_ok = add_ok,
    rem_ok = rem_ok,
    clado_isl = clado_isl,
    Wd = Wd,
    Wi = Wi,
    wide = wide,
    vic_w = vic_w
  )
  return(rates)
}

#' Islands whose recolonisation resets the colonisation time
#'
#' @description The colonisation time of a non-endemic lineage is reset only
#'   if it occupies a single island and that island is recolonised.
#'
#' @param occ 0/1 occupancy vector of the non-endemic lineage.
#' @return 0/1 vector: \code{occ} if the lineage occupies one island, zeros
#'   otherwise.
#' @keywords internal
DAISIE_A_reset_islands <- function(occ) {
  if (sum(occ) == 1) occ else 0 * occ
}
