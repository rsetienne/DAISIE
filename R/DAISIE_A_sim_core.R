# Core of the archipelago simulation (DAISIE_A).
#
# State of the simulation
# -----------------------
# st$spec  numeric matrix, one row per extant lineage (species) on the
#          archipelago, with columns
#            1 id      species ID (mainland ID for non-endemics, as in DAISIE)
#            2 anc     mainland ancestor
#            3 col_t   colonisation time (forward time since island birth)
#            4 type    1 = I (non-endemic), 2 = A (anagenetic/relict endemic),
#                      3 = C (cladogenetic endemic)
#            5 ana     anagenetic origin: 0 = NA, 1 = "Immig_parent",
#                      2 = "Clado_extinct"
#            6 node    genealogy node at the start of the current branch
#            7 side    side of that branch at its parent node (1 = A, 2 = B)
#            8 clade   ID of the colonisation event the lineage descends from
#            9..(8+n)  occupancy (0/1) of islands 1..n: the range
# st$gen   genealogy (full tree, including extinct lineages) with columns
#            node, parent, time, kind, side, clade, species
#          kind: 1 = colonisation (root), 2 = speciation, 3 = extinction,
#                4 = replaced by recolonisation, 5 = extant tip
#
# The row order of st$spec, the species IDs and the order of random number
# calls follow DAISIE_sim_core_cr(), so that for a single island the
# simulation reproduces DAISIE_sim_cr() exactly for the same seed.

.DA_ID <- 1L
.DA_ANC <- 2L
.DA_COLT <- 3L
.DA_TYPE <- 4L
.DA_ANA <- 5L
.DA_NODE <- 6L
.DA_SIDE <- 7L
.DA_CLADE <- 8L
.DA_NFIX <- 8L

#' Internal function of the DAISIE_A simulation: one clade (or island-wide)
#'
#' @description Simulates the archipelago model for \code{mainland_n}
#'   mainland species with a Gillespie algorithm. With \code{mainland_n = 1}
#'   this simulates one clade under clade-specific diversity dependence.
#'
#' @param time Simulation time (age of the oldest island).
#' @param mainland_n Number of mainland species.
#' @param pars A \code{DAISIE_A_pars} object, see
#'   \code{\link{DAISIE_A_create_pars}()}.
#' @param split_cache Environment used to cache admissible vicariance splits.
#'
#' @return A list as returned by \code{\link{DAISIE_create_island}()}
#'   (\code{stt_table}, \code{branching_times}, \code{stac},
#'   \code{missing_species}, and \code{all_colonisations} if relevant), where
#'   the archipelago is treated as a single island, with the additional
#'   archipelago elements \code{stt_islands}, \code{species} and
#'   \code{clades}, see \code{\link{DAISIE_A_sim}()}.
#' @keywords internal
DAISIE_A_sim_core <- function(time,
                              mainland_n = 1,
                              pars,
                              split_cache = new.env()) {
  total_time <- time
  testit::assert(DAISIE_A_are_pars(pars))
  if (all(pars$gam == 0)) {
    stop("Archipelago has no species and the rate of
    colonisation is zero. Archipelago cannot be colonised.")
  }
  n <- pars$n_islands
  emerge <- DAISIE_A_emergence_times(pars, total_time)
  submerge <- DAISIE_A_submergence_times(pars, total_time)

  st <- DAISIE_A_init_state(n_islands = n, mainland_n = mainland_n)
  stt <- DAISIE_A_init_stt(total_time = total_time, n_islands = n)

  timeval <- 0
  while (timeval < total_time) {
    rates <- DAISIE_A_update_rates(
      st = st,
      timeval = timeval,
      pars = pars,
      emerge = emerge,
      mainland_spec = st$mainland_spec,
      split_cache = split_cache,
      submerge = submerge
    )
    if (rates$total != 0) {
      dt <- stats::rexp(1, rates$total)
      new_time <- timeval + dt
    } else {
      new_time <- total_time
    }
    later <- c(emerge[emerge > timeval], submerge[submerge > timeval])
    next_island_event <- if (length(later) > 0) min(later) else Inf
    if (next_island_event < total_time && new_time >= next_island_event) {
      # An island emerges or submerges before the next event: the rates
      # change. Because the process is memoryless, we can move to that time
      # and draw a new waiting time. A submerging island loses all its
      # populations.
      timeval <- next_island_event
      sinking <- which(submerge == timeval)
      if (length(sinking) > 0) {
        st <- DAISIE_A_submerge(st, sinking, timeval)
        stt <- DAISIE_A_add_stt_row(stt, st, total_time - timeval, n)
      }
      next
    }
    timeval <- new_time
    if (timeval < total_time) {
      possible_event <- DAISIE_A_sample_event(rates)
      st <- DAISIE_A_update_state(
        st = st,
        possible_event = possible_event,
        rates = rates,
        timeval = timeval,
        pars = pars,
        split_cache = split_cache
      )
      stt <- DAISIE_A_add_stt_row(stt, st, total_time - timeval, n)
    }
  }
  stt_table <- stt$table[seq_len(stt$n), , drop = FALSE]
  stt_islands <- stt$islands[seq_len(stt$n), , drop = FALSE]
  stt_table <- rbind(stt_table, c(0, stt_table[nrow(stt_table), 2:4]))
  stt_islands <- rbind(stt_islands,
                       c(0, stt_islands[nrow(stt_islands), -1]))
  island <- DAISIE_A_create_island(
    st = st,
    stt_table = stt_table,
    stt_islands = stt_islands,
    total_time = total_time,
    mainland_n = mainland_n,
    pars = pars
  )
  return(island)
}

#' Emergence times of the islands in forward time
#'
#' @param pars A \code{DAISIE_A_pars} object.
#' @param total_time Simulation time.
#' @return Numeric vector with, for each island, the (forward) time at which
#'   it emerges; 0 if it is as old as the simulation.
#' @keywords internal
DAISIE_A_emergence_times <- function(pars, total_time) {
  if (is.null(pars$island_ages)) return(rep(0, pars$n_islands))
  if (any(pars$island_ages > total_time + 1e-12)) {
    stop("Island ages cannot be larger than the simulation time.")
  }
  pmax(0, total_time - pars$island_ages)
}

#' Submergence times of the islands in forward time
#'
#' @param pars A \code{DAISIE_A_pars} object.
#' @param total_time Simulation time.
#' @return Numeric vector with, for each island, the (forward) time at which
#'   it submerges; \code{Inf} if it still exists.
#' @keywords internal
DAISIE_A_submergence_times <- function(pars, total_time) {
  if (is.null(pars$submergence_ages)) return(rep(Inf, pars$n_islands))
  sub <- pars$submergence_ages
  if (any(sub >= total_time, na.rm = TRUE)) {
    stop("Submergence ages must be smaller than the simulation time.")
  }
  out <- total_time - sub
  out[is.na(out)] <- Inf
  return(out)
}

#' Submergence of islands in the simulation
#'
#' @description All populations on the submerging islands are lost: lineages
#'   lose these islands from their range, and lineages confined to them go
#'   extinct (with the DAISIE singleton rule applied to their clades).
#'
#' @param st State list.
#' @param islands Indices of the submerging islands.
#' @param timeval Current (forward) time.
#' @return Updated state list.
#' @keywords internal
DAISIE_A_submerge <- function(st, islands, timeval) {
  for (i in islands) {
    if (nrow(st$spec) == 0) next
    on_i <- which(st$spec[, .DA_NFIX + i] == 1)
    if (length(on_i) == 0) next
    size <- rowSums(st$spec[, .DA_NFIX + seq_len(st$n_islands), drop = FALSE])
    dying <- on_i[size[on_i] == 1]
    st$spec[on_i, .DA_NFIX + i] <- 0
    for (r in dying) {
      st <- DAISIE_A_add_node(st, parent = st$spec[r, .DA_NODE],
                              time = timeval, kind = 3,
                              side = st$spec[r, .DA_SIDE],
                              clade = st$spec[r, .DA_CLADE],
                              species = st$spec[r, .DA_ID])
    }
    if (length(dying) > 0) {
      clades <- unique(st$spec[dying, .DA_CLADE])
      st$spec <- st$spec[-dying, , drop = FALSE]
      for (cl in clades) st <- DAISIE_A_clade_singleton(st, cl)
    }
  }
  rownames(st$spec) <- NULL
  return(st)
}

#' Initial (empty) state of the DAISIE_A simulation
#'
#' @param n_islands Number of islands.
#' @param mainland_n Number of mainland species.
#' @return The state list \code{st}.
#' @keywords internal
DAISIE_A_init_state <- function(n_islands, mainland_n) {
  spec <- matrix(numeric(0), nrow = 0, ncol = .DA_NFIX + n_islands)
  gen <- matrix(NA_real_, nrow = 64, ncol = 7)
  colnames(gen) <- c("node", "parent", "time", "kind", "side", "clade",
                     "species")
  list(
    spec = spec,
    gen = gen,
    n_gen = 0L,
    n_clades = 0L,
    maxspecID = mainland_n,
    mainland_spec = seq_len(mainland_n),
    n_islands = n_islands
  )
}

#' Add a node to the genealogy
#'
#' @param st State list.
#' @param parent,time,kind,side,clade,species Node attributes.
#' @return Updated state; the new node has ID \code{st$n_gen}.
#' @keywords internal
DAISIE_A_add_node <- function(st, parent, time, kind, side, clade,
                              species = NA) {
  st$n_gen <- st$n_gen + 1L
  if (st$n_gen > nrow(st$gen)) {
    st$gen <- rbind(st$gen, matrix(NA_real_, nrow = nrow(st$gen),
                                   ncol = ncol(st$gen)))
  }
  st$gen[st$n_gen, ] <- c(st$n_gen, parent, time, kind, side, clade,
                          species)
  return(st)
}

#' Initialise the species-through-time tables
#'
#' @param total_time Simulation time.
#' @param n_islands Number of islands.
#' @return List with the archipelago-level table (\code{Time}, \code{nI},
#'   \code{nA}, \code{nC}, as in DAISIE) and the per-island table.
#' @keywords internal
DAISIE_A_init_stt <- function(total_time, n_islands) {
  tab <- matrix(NA_real_, nrow = 64, ncol = 4)
  colnames(tab) <- c("Time", "nI", "nA", "nC")
  tab[1, ] <- c(total_time, 0, 0, 0)
  isl <- matrix(NA_real_, nrow = 64, ncol = 1 + 4 * n_islands)
  colnames(isl) <- c("Time", paste0(rep(c("nI_", "nA_", "nC_", "nSIE_"),
                                        n_islands),
                                    rep(seq_len(n_islands), each = 4)))
  isl[1, ] <- c(total_time, rep(0, 4 * n_islands))
  list(table = tab, islands = isl, n = 1L)
}

#' Per-island species counts
#'
#' @param spec The species matrix of the state.
#' @param n_islands Number of islands.
#' @return Vector \code{nI_1, nA_1, nC_1, nSIE_1, nI_2, ...}: numbers of
#'   non-endemic, anagenetic, cladogenetic species and single-island endemics
#'   on each island.
#' @keywords internal
DAISIE_A_island_counts <- function(spec, n_islands) {
  if (nrow(spec) == 0) return(rep(0, 4 * n_islands))
  occ <- spec[, .DA_NFIX + seq_len(n_islands), drop = FALSE]
  type <- spec[, .DA_TYPE]
  single_end <- type != 1 & rowSums(occ) == 1
  counts <- rbind(
    colSums(occ[type == 1, , drop = FALSE]),
    colSums(occ[type == 2, , drop = FALSE]),
    colSums(occ[type == 3, , drop = FALSE]),
    colSums(occ[single_end, , drop = FALSE])
  )
  as.vector(counts)
}

#' Add a row to the species-through-time tables
#'
#' @param stt List of STT tables.
#' @param st State list.
#' @param time_bp Time before present.
#' @param n_islands Number of islands.
#' @return Updated \code{stt}.
#' @keywords internal
DAISIE_A_add_stt_row <- function(stt, st, time_bp, n_islands) {
  stt$n <- stt$n + 1L
  if (stt$n > nrow(stt$table)) {
    stt$table <- rbind(stt$table, matrix(NA_real_, nrow = nrow(stt$table),
                                         ncol = ncol(stt$table)))
    stt$islands <- rbind(stt$islands,
                         matrix(NA_real_, nrow = nrow(stt$islands),
                                ncol = ncol(stt$islands)))
  }
  type <- st$spec[, .DA_TYPE]
  stt$table[stt$n, ] <- c(time_bp,
                          length(which(type == 1)),
                          length(which(type == 2)),
                          length(which(type == 3)))
  stt$islands[stt$n, ] <- c(time_bp,
                            DAISIE_A_island_counts(st$spec, n_islands))
  return(stt)
}

#' Sample one element with (optional) weights
#'
#' @description Uses \code{DDD::sample2()} without weights when all weights
#'   are equal, so that the random number calls are identical to those of the
#'   standard DAISIE simulation.
#'
#' @param x Candidates.
#' @param w Non-negative weights.
#' @return One element of \code{x}.
#' @keywords internal
DAISIE_A_sample_one <- function(x, w = NULL) {
  if (is.null(w) || all(w == w[1])) {
    return(DDD::sample2(x, 1))
  }
  return(DDD::sample2(x, 1, prob = w))
}

#' Sample which type of event happens
#'
#' @param rates Rates as returned by \code{\link{DAISIE_A_update_rates}()}.
#' @return Integer: 1 colonisation, 2 extinction (local or global),
#'   3 anagenesis, 4 in-situ cladogenesis, 5 dispersal, 6 vicariance,
#'   7 island-level anagenesis.
#' @keywords internal
DAISIE_A_sample_event <- function(rates) {
  if (rates$disp + rates$vic + rates$isl_ana == 0) {
    # identical call as DAISIE_sample_event_cr()
    possible_event <- sample(
      x = 1:4,
      size = 1,
      replace = FALSE,
      prob = c(rates$immig, rates$ext, rates$ana, rates$clado)
    )
  } else {
    possible_event <- sample(
      x = 1:7,
      size = 1,
      replace = FALSE,
      prob = c(rates$immig, rates$ext, rates$ana, rates$clado,
               rates$disp, rates$vic, rates$isl_ana)
    )
  }
  return(possible_event)
}
