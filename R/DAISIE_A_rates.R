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
#'       \eqn{\gamma_i f_i e_i};
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
#'     \item founder-event speciation from range \eqn{R} to island
#'       \eqn{j \notin R}: \eqn{\sum_{i \in R} J_{ij} f_j e_j}.
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
#'   \code{ana}, \code{clado}, \code{disp}, \code{vic}, \code{jump}, their sum
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
      # recolonisation of an island in the range (reset of colonisation
      # time) or of a new island (range expansion), if allowed
      w_col[m] <- sum(col_isl * (occ[r, ] + add_ok[r, ]))
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

  # Founder-event (jump) speciation
  Wj <- NULL
  jump <- 0
  if (L > 0 && any(pars$J > 0)) {
    Wj <- (occ %*% pars$J) * (1 - occ) * rep(f * exists, each = L)
    jump <- sum(Wj)
  }

  total <- immig + ext + ana + clado + disp + vic + jump
  rates <- list(
    immig = immig,
    ext = ext,
    ana = ana,
    clado = clado,
    disp = disp,
    vic = vic,
    jump = jump,
    total = total,
    w_col = w_col,
    col_isl = col_isl,
    add_ok = add_ok,
    rem_ok = rem_ok,
    clado_isl = clado_isl,
    Wd = Wd,
    Wj = Wj,
    wide = wide,
    vic_w = vic_w
  )
  return(rates)
}
