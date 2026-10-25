#' Update the state of the archipelago given the sampled event
#'
#' @description Internal function. Makes the sampled event of the DAISIE_A
#'   model happen: it updates the species matrix (including the ranges), the
#'   species IDs and the genealogy.
#'
#'   Conventions (all following DAISIE where the archipelago has one island):
#'   \itemize{
#'     \item Colonisation: a mainland species colonises one island. If it is
#'       already present as a non-endemic lineage on that island only, that
#'       lineage is replaced by the new colonist (its colonisation time is
#'       reset, as in DAISIE). If it is present as a non-endemic lineage
#'       elsewhere, the island is added to its range and the colonisation
#'       time is not reset. (Recolonisation of an island that is part of a
#'       widespread range has no effect and has rate zero.)
#'     \item Extinction: the population on one island goes extinct; if this
#'       was the last island, the species goes extinct.
#'     \item Anagenesis: a non-endemic lineage (on all its islands) becomes a
#'       new endemic species.
#'     \item In-situ cladogenesis on island \eqn{i}: a single-island lineage
#'       splits into two new endemic species (\code{clado_single = "split"})
#'       or buds (\code{"bud"}); a widespread lineage buds: it keeps its
#'       identity and range, and a new endemic species on \eqn{\{i\}} appears.
#'     \item Dispersal: island \eqn{j} is added to the range.
#'     \item Vicariance: the range is split into two disjoint ranges, each
#'       occupied by a new endemic species. The lineage is chosen with weight
#'       1 (\code{vic_per = "lineage"}) or proportional to its number of
#'       admissible splits (\code{vic_per = "split"}); the split is then
#'       chosen uniformly among the admissible splits.
#'     \item Founder event: the lineage keeps its identity and range, and a
#'       new endemic species appears on island \eqn{j} outside the range.
#'   }
#'
#' @param st State list.
#' @param possible_event Integer event class, see
#'   \code{\link{DAISIE_A_sample_event}()}.
#' @param rates Output of \code{\link{DAISIE_A_update_rates}()}.
#' @param timeval Current (forward) time.
#' @param pars A \code{DAISIE_A_pars} object.
#' @param split_cache Environment for caching vicariance splits.
#'
#' @return The updated state list.
#' @keywords internal
DAISIE_A_update_state <- function(st,
                                  possible_event,
                                  rates,
                                  timeval,
                                  pars,
                                  split_cache) {
  n <- pars$n_islands
  occ_cols <- .DA_NFIX + seq_len(n)
  ##########################################
  # COLONISATION
  if (possible_event == 1) {
    colonist <- DAISIE_A_sample_one(st$mainland_spec, rates$w_col)
    isitthere <- which(st$spec[, .DA_ID] == colonist)
    testit::assert(length(isitthere) <= 1)
    if (length(isitthere) == 0) {
      w_isl <- rates$col_isl
    } else {
      w_isl <- rates$col_isl *
        (DAISIE_A_reset_islands(st$spec[isitthere, occ_cols]) +
           rates$add_ok[isitthere, ])
    }
    island <- if (n == 1) 1L else DAISIE_A_sample_one(seq_len(n), w_isl)
    if (length(isitthere) == 0) {
      st$n_clades <- st$n_clades + 1L
      st <- DAISIE_A_add_node(st, parent = 0, time = timeval, kind = 1,
                              side = 0, clade = st$n_clades,
                              species = colonist)
      new_occ <- rep(0, n)
      new_occ[island] <- 1
      st$spec <- rbind(st$spec, c(colonist, colonist, timeval, 1, 0,
                                  st$n_gen, 0, st$n_clades, new_occ))
    } else if (st$spec[isitthere, .DA_NFIX + island] == 0) {
      # range expansion of the non-endemic lineage, no reset
      st$spec[isitthere, .DA_NFIX + island] <- 1
    } else {
      # recolonisation of the only island of the non-endemic lineage
      r <- isitthere
      testit::assert(sum(st$spec[r, occ_cols]) == 1)
      old_clade <- st$spec[r, .DA_CLADE]
      # the old non-endemic lineage is replaced by the new colonist
      st <- DAISIE_A_add_node(st, parent = st$spec[r, .DA_NODE],
                              time = timeval, kind = 4,
                              side = st$spec[r, .DA_SIDE], clade = old_clade,
                              species = st$spec[r, .DA_ID])
      st$n_clades <- st$n_clades + 1L
      st <- DAISIE_A_add_node(st, parent = 0, time = timeval, kind = 1,
                              side = 0, clade = st$n_clades,
                              species = colonist)
      new_occ <- st$spec[r, occ_cols]
      st$spec[r, ] <- c(colonist, colonist, timeval, 1, 0, st$n_gen, 0,
                        st$n_clades, new_occ)
      st <- DAISIE_A_clade_singleton(st, old_clade)
    }
  }
  ##########################################
  # EXTINCTION (local, or global if the last population is lost)
  if (possible_event == 2) {
    idx <- which(rates$rem_ok == 1, arr.ind = TRUE)
    k <- DAISIE_A_sample_one(seq_len(nrow(idx)), pars$mu[idx[, 2]])
    r <- idx[k, 1]
    island <- idx[k, 2]
    if (sum(st$spec[r, occ_cols]) == 1) {
      clade <- st$spec[r, .DA_CLADE]
      st <- DAISIE_A_add_node(st, parent = st$spec[r, .DA_NODE],
                              time = timeval, kind = 3,
                              side = st$spec[r, .DA_SIDE], clade = clade,
                              species = st$spec[r, .DA_ID])
      st$spec <- st$spec[-r, , drop = FALSE]
      st <- DAISIE_A_clade_singleton(st, clade)
    } else {
      st$spec[r, .DA_NFIX + island] <- 0
    }
  }
  ##########################################
  # ANAGENESIS
  if (possible_event == 3) {
    immi_specs <- which(st$spec[, .DA_TYPE] == 1)
    # we only allow immigrants to undergo anagenesis
    if (length(immi_specs) == 1) {
      anagenesis <- immi_specs
    }
    if (length(immi_specs) > 1) {
      anagenesis <- DDD::sample2(immi_specs, 1)
    }
    st$maxspecID <- st$maxspecID + 1
    st$spec[anagenesis, .DA_TYPE] <- 2
    st$spec[anagenesis, .DA_ID] <- st$maxspecID
    st$spec[anagenesis, .DA_ANA] <- 1
  }
  ##########################################
  # IN-SITU CLADOGENESIS (within one island)
  if (possible_event == 4) {
    idx <- which(st$spec[, occ_cols, drop = FALSE] == 1, arr.ind = TRUE)
    k <- DAISIE_A_sample_one(seq_len(nrow(idx)), rates$clado_isl[idx[, 2]])
    r <- idx[k, 1]
    island <- idx[k, 2]
    parent <- st$spec[r, ]
    st <- DAISIE_A_add_node(st, parent = parent[.DA_NODE], time = timeval,
                            kind = 2, side = parent[.DA_SIDE],
                            clade = parent[.DA_CLADE])
    v <- st$n_gen
    if (sum(parent[occ_cols]) == 1 && pars$clado_single == "split") {
      # as in DAISIE: two new (endemic) species replace the parent
      st$spec[r, c(.DA_ID, .DA_TYPE, .DA_ANA, .DA_NODE, .DA_SIDE)] <-
        c(st$maxspecID + 1, 3, 0, v, 1)
      daughter_b <- parent
      daughter_b[c(.DA_ID, .DA_TYPE, .DA_ANA, .DA_NODE, .DA_SIDE)] <-
        c(st$maxspecID + 2, 3, 0, v, 2)
      st$spec <- rbind(st$spec, daughter_b)
      st$maxspecID <- st$maxspecID + 2
    } else {
      # budding: the parent keeps its identity (and range), one new endemic
      # species appears on this island only
      st$spec[r, c(.DA_NODE, .DA_SIDE)] <- c(v, 1)
      bud <- parent
      bud[c(.DA_ID, .DA_TYPE, .DA_ANA, .DA_NODE, .DA_SIDE)] <-
        c(st$maxspecID + 1, 3, 0, v, 2)
      bud[occ_cols] <- 0
      bud[.DA_NFIX + island] <- 1
      st$spec <- rbind(st$spec, bud)
      st$maxspecID <- st$maxspecID + 1
    }
  }
  ##########################################
  # DISPERSAL (range expansion)
  if (possible_event == 5) {
    idx <- which(rates$Wd > 0, arr.ind = TRUE)
    k <- DAISIE_A_sample_one(seq_len(nrow(idx)), rates$Wd[idx])
    st$spec[idx[k, 1], .DA_NFIX + idx[k, 2]] <- 1
  }
  ##########################################
  # VICARIANCE (between-island speciation)
  if (possible_event == 6) {
    wide <- rates$wide
    r <- if (length(wide) == 1) wide else DAISIE_A_sample_one(wide, rates$vic_w)
    parent <- st$spec[r, ]
    set <- which(parent[occ_cols] == 1)
    part_a <- DAISIE_A_sample_split(set, pars, split_cache)
    part_b <- setdiff(set, part_a)
    st <- DAISIE_A_add_node(st, parent = parent[.DA_NODE], time = timeval,
                            kind = 2, side = parent[.DA_SIDE],
                            clade = parent[.DA_CLADE])
    v <- st$n_gen
    occ_a <- rep(0, n)
    occ_a[part_a] <- 1
    occ_b <- rep(0, n)
    occ_b[part_b] <- 1
    st$spec[r, c(.DA_ID, .DA_TYPE, .DA_ANA, .DA_NODE, .DA_SIDE)] <-
      c(st$maxspecID + 1, 3, 0, v, 1)
    st$spec[r, occ_cols] <- occ_a
    daughter_b <- parent
    daughter_b[c(.DA_ID, .DA_TYPE, .DA_ANA, .DA_NODE, .DA_SIDE)] <-
      c(st$maxspecID + 2, 3, 0, v, 2)
    daughter_b[occ_cols] <- occ_b
    st$spec <- rbind(st$spec, daughter_b)
    st$maxspecID <- st$maxspecID + 2
  }
  ##########################################
  # FOUNDER-EVENT (JUMP) SPECIATION
  if (possible_event == 7) {
    idx <- which(rates$Wj > 0, arr.ind = TRUE)
    k <- DAISIE_A_sample_one(seq_len(nrow(idx)), rates$Wj[idx])
    r <- idx[k, 1]
    island <- idx[k, 2]
    parent <- st$spec[r, ]
    st <- DAISIE_A_add_node(st, parent = parent[.DA_NODE], time = timeval,
                            kind = 2, side = parent[.DA_SIDE],
                            clade = parent[.DA_CLADE])
    v <- st$n_gen
    st$spec[r, c(.DA_NODE, .DA_SIDE)] <- c(v, 1)
    founder <- parent
    founder[c(.DA_ID, .DA_TYPE, .DA_ANA, .DA_NODE, .DA_SIDE)] <-
      c(st$maxspecID + 1, 3, 0, v, 2)
    founder[occ_cols] <- 0
    founder[.DA_NFIX + island] <- 1
    st$spec <- rbind(st$spec, founder)
    st$maxspecID <- st$maxspecID + 1
  }
  rownames(st$spec) <- NULL
  return(st)
}

#' Relabel the last surviving cladogenetic species of a clade
#'
#' @description As in DAISIE, if only one species of a clade is left and it
#'   is of cladogenetic origin, it becomes an anagenetic species with
#'   anagenetic origin "Clado_extinct".
#'
#' @param st State list.
#' @param clade Clade (colonisation event) ID.
#' @return Updated state list.
#' @keywords internal
DAISIE_A_clade_singleton <- function(st, clade) {
  rows <- which(st$spec[, .DA_CLADE] == clade)
  if (length(rows) == 1 && st$spec[rows, .DA_TYPE] == 3) {
    st$spec[rows, .DA_TYPE] <- 2
    st$spec[rows, .DA_ANA] <- 2
  }
  return(st)
}
