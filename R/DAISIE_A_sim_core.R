# Core of the archipelago simulation (DAISIE_A, version 10: divergence states).
#
# DAISIE_A_sim_core_R() is the reference implementation in R;
# src/DAISIE_A_sim.cpp implements the same algorithm in C++ with the same
# sequence of random numbers, so that both give identical output for the same
# seed (DAISIE_A_sim(..., use_cpp = TRUE / FALSE)). Any change here must be
# made there too.
#
# A lineage is a list with
#   id, anc, colt, type (1 = I non-endemic, 2 = A anagenetic, 3 = C
#   cladogenetic), ana (0, 1 = "Immig_parent", 2 = "Clado_extinct"), node
#   (genealogy node at the start of its current branch), side (1 = A, 2 = B),
#   clade (colonization event);
#   occ   0/1 vector of length n: the range;
#   S     (n + 1) x (n + 1) matrix of divergence states of the tracked pairs
#         (node n + 1 is the mainland), -1 = not tracked;
#   Tau   matching matrix with the time of the first step 0 -> 1 after the
#         last reset (-Inf while the state is 0);
#   s     within-island states (protracted option), stau their 0 -> 1 times.
# Tracked pairs: all pairs of adjacent nodes in the range (mainland and island
# i are adjacent if gam_i > 0, and only for non-endemics), plus gap pairs
# between parts of a range that are not connected by adjacent pairs.
#
# The genealogy (st$gen) has columns node, parent, time (as reported: the
# split, or with protracted = TRUE the start of divergence), kind
# (1 colonization, 2 speciation, 3 extinction, 4 replaced by recolonization,
# 5 extant tip), side, clade, species and time_split (time of the event).

.DA_ID <- 1L
.DA_ANC <- 2L
.DA_COLT <- 3L
.DA_TYPE <- 4L
.DA_ANA <- 5L
.DA_NODE <- 6L
.DA_SIDE <- 7L
.DA_CLADE <- 8L
.DA_NFIX <- 8L
.DA_GEN_COLS <- c("node", "parent", "time", "kind", "side", "clade",
                  "species", "time_split")

#' Settings of the DAISIE_A simulation core
#'
#' @param pars A \code{DAISIE_A_pars} object.
#' @param total_time Simulation time.
#' @param protracted Logical, see \code{\link{DAISIE_A_sim}()}.
#' @return List of settings used by both implementations of the core.
#' @keywords internal
DAISIE_A_core_cfg <- function(pars, total_time, protracted) {
  n <- pars$n_islands
  n1 <- n + 1L
  A1 <- matrix(FALSE, n1, n1)
  A1[seq_len(n), seq_len(n)] <- pars$adjacency
  A1[n1, seq_len(n)] <- pars$gam > 0
  A1[seq_len(n), n1] <- pars$gam > 0
  diag(A1) <- FALSE
  list(n = n, lac = pars$lac, mu = pars$mu, K = pars$K, gam = pars$gam,
       beta = pars$beta, alpha = pars$alpha, k = pars$k, D = pars$D, A1 = A1,
       split_single = pars$clado_single == "split",
       protracted = isTRUE(protracted),
       emerge = DAISIE_A_emergence_times(pars, total_time),
       submerge = DAISIE_A_submergence_times(pars, total_time),
       total_time = total_time)
}

#' Internal function of the DAISIE_A simulation: one clade (or island-wide)
#'
#' @description Simulates the archipelago model for \code{mainland_n}
#'   mainland species with a Gillespie algorithm and converts the result into
#'   the output of \code{\link{DAISIE_A_sim}()}. With \code{mainland_n = 1}
#'   this simulates one clade under clade-specific diversity dependence.
#'
#' @param time Simulation time (age of the oldest island).
#' @param mainland_n Number of mainland species.
#' @param pars A \code{DAISIE_A_pars} object, see
#'   \code{\link{DAISIE_A_create_pars}()}.
#' @param protracted Logical, see \code{\link{DAISIE_A_sim}()}.
#' @param use_cpp Logical; use the C++ implementation (default) or the R
#'   reference implementation. Both give the same result for the same seed.
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
                              protracted = FALSE,
                              use_cpp = TRUE) {
  testit::assert(DAISIE_A_are_pars(pars))
  if (all(pars$gam == 0)) {
    stop("Archipelago has no species and the rate of
    colonization is zero. Archipelago cannot be colonized.")
  }
  core <- if (use_cpp) {
    DAISIE_A_sim_core_cpp(time = time, mainland_n = mainland_n, pars = pars,
                          protracted = protracted)
  } else {
    DAISIE_A_sim_core_R(time = time, mainland_n = mainland_n, pars = pars,
                        protracted = protracted)
  }
  DAISIE_A_create_island(core = core, total_time = time,
                         mainland_n = mainland_n, pars = pars,
                         protracted = protracted)
}

#' C++ implementation of the DAISIE_A simulation core
#'
#' @inheritParams DAISIE_A_sim_core
#' @return See \code{\link{DAISIE_A_sim_core_R}()}.
#' @keywords internal
DAISIE_A_sim_core_cpp <- function(time, mainland_n, pars, protracted) {
  cfg <- DAISIE_A_core_cfg(pars, time, protracted)
  cfg$mainland_n <- as.integer(mainland_n)
  out <- .Call("DAISIE_A_sim_cpp", cfg, PACKAGE = "DAISIE")
  colnames(out$gen) <- .DA_GEN_COLS
  colnames(out$stt_table) <- c("Time", "nI", "nA", "nC")
  colnames(out$stt_islands) <- DAISIE_A_stt_island_names(cfg$n)
  out
}

#' R reference implementation of the DAISIE_A simulation core
#'
#' @inheritParams DAISIE_A_sim_core
#' @return A list with \code{spec} (one row per extant lineage: ID, mainland
#'   ancestor, colonization time, type, anagenetic origin, genealogy node,
#'   side, clade, and the range as 0/1 columns), \code{gen} (genealogy),
#'   \code{stt_table}, \code{stt_islands} and \code{n_events}.
#' @keywords internal
DAISIE_A_sim_core_R <- function(time, mainland_n, pars, protracted) {
  total_time <- time
  cfg <- DAISIE_A_core_cfg(pars, total_time, protracted)
  n <- cfg$n
  st <- list(lin = list(),
             gen = matrix(NA_real_, nrow = 64, ncol = 8,
                          dimnames = list(NULL, .DA_GEN_COLS)),
             n_gen = 0L, n_clades = 0L, maxspecID = mainland_n,
             mainland_n = mainland_n)
  stt <- DAISIE_A_init_stt(total_time = total_time, n_islands = n)
  timeval <- 0
  n_events <- 0
  debug <- isTRUE(getOption("DAISIE_A.debug", FALSE))
  while (timeval < total_time) {
    ev <- da_events(st, cfg, timeval)
    total <- if (ev$n > 0) ev$cs[ev$n] else 0
    if (total > 0) {
      new_time <- timeval + stats::rexp(1, total)
    } else {
      new_time <- total_time
    }
    later <- c(cfg$emerge[cfg$emerge > timeval],
               cfg$submerge[cfg$submerge > timeval])
    next_island_event <- if (length(later) > 0) min(later) else Inf
    if (next_island_event < total_time && new_time >= next_island_event) {
      # An island emerges or submerges before the next event: the rates
      # change. Because the process is memoryless, we move to that time and
      # draw a new waiting time.
      timeval <- next_island_event
      sinking <- which(cfg$submerge == timeval)
      for (i in sinking) st <- da_submerge(st, i, timeval, cfg)
      if (debug) da_check_invariants(st, cfg)
      if (length(sinking) > 0) {
        stt <- da_stt_row(stt, st, total_time - timeval, cfg)
      }
      next
    }
    timeval <- new_time
    if (timeval < total_time) {
      target <- stats::runif(1) * total
      e <- 1L
      while (ev$cs[e] <= target) e <- e + 1L
      st <- da_apply(st, ev, e, timeval, cfg)
      if (debug) da_check_invariants(st, cfg)
      n_events <- n_events + 1
      stt <- da_stt_row(stt, st, total_time - timeval, cfg)
    }
  }
  stt_table <- stt$table[seq_len(stt$n), , drop = FALSE]
  stt_islands <- stt$islands[seq_len(stt$n), , drop = FALSE]
  stt_table <- rbind(stt_table, c(0, stt_table[nrow(stt_table), 2:4]))
  stt_islands <- rbind(stt_islands, c(0, stt_islands[nrow(stt_islands), -1]))
  L <- length(st$lin)
  spec <- matrix(NA_real_, nrow = L, ncol = .DA_NFIX + n)
  for (r in seq_len(L)) {
    x <- st$lin[[r]]
    spec[r, ] <- c(x$id, x$anc, x$colt, x$type, x$ana, x$node, x$side,
                   x$clade, x$occ)
  }
  list(spec = spec, gen = st$gen[seq_len(st$n_gen), , drop = FALSE],
       stt_table = stt_table, stt_islands = stt_islands,
       n_events = n_events)
}

#' Events of the DAISIE_A simulation and their rates
#'
#' @description Lists all events with a positive rate, in a fixed order that
#'   the C++ implementation follows too: colonizations (per mainland species
#'   and island), then per lineage its local extinctions, within-island
#'   speciation (or within-island steps of the protracted option), dispersal
#'   events and divergence steps. Events that would change nothing (gene flow
#'   between populations whose pair is at 0) are left out.
#'
#' @param st State list.
#' @param cfg Settings, see \code{\link{DAISIE_A_core_cfg}()}.
#' @param timeval Current (forward) time.
#' @return List with the number of events \code{n}, their cumulative rates
#'   \code{cs} and the event descriptors \code{type} (1 colonization, 2
#'   extinction, 3 within-island, 4 dispersal, 5 divergence step), \code{r}
#'   (lineage, or mainland species for colonization), \code{a}, \code{b}.
#' @keywords internal
da_events <- function(st, cfg, timeval) {
  n <- cfg$n
  n1 <- n + 1L
  k <- cfg$k
  exists <- cfg$emerge <= timeval & cfg$submerge > timeval
  ni <- numeric(n)
  for (L in st$lin) ni <- ni + L$occ
  f <- pmax(0, 1 - ni / cfg$K)
  ne_row <- integer(st$mainland_n)
  for (r in seq_along(st$lin)) {
    if (st$lin[[r]]$type == 1) ne_row[st$lin[[r]]$id] <- r
  }
  cap <- 64L
  rate <- numeric(cap)
  type <- integer(cap)
  er <- integer(cap)
  ea <- integer(cap)
  eb <- integer(cap)
  m_ev <- 0L
  add <- function(x, ty, r, a, b = 0L) {
    m_ev <<- m_ev + 1L
    if (m_ev > length(rate)) {
      rate <<- c(rate, numeric(length(rate)))
      type <<- c(type, integer(length(type)))
      er <<- c(er, integer(length(er)))
      ea <<- c(ea, integer(length(ea)))
      eb <<- c(eb, integer(length(eb)))
    }
    rate[m_ev] <<- x
    type[m_ev] <<- ty
    er[m_ev] <<- r
    ea[m_ev] <<- a
    eb[m_ev] <<- b
  }
  for (m in seq_len(st$mainland_n)) {
    for (i in seq_len(n)) {
      if (!exists[i] || cfg$gam[i] == 0) next
      r <- ne_row[m]
      if (r == 0L) {
        x <- cfg$gam[i] * f[i]
      } else {
        L <- st$lin[[r]]
        if (L$occ[i] == 0L) {
          x <- cfg$gam[i] * f[i]
        } else if (sum(L$occ) == 1L) {
          x <- cfg$gam[i]
        } else if (L$S[i, n1] > 0 || (cfg$protracted && L$s[i] > 0L)) {
          x <- cfg$gam[i]
        } else {
          x <- 0
        }
      }
      if (x > 0) add(x, 1L, m, i)
    }
  }
  for (r in seq_along(st$lin)) {
    L <- st$lin[[r]]
    R <- which(L$occ == 1L)
    for (i in R) {
      if (cfg$mu[i] > 0) add(cfg$mu[i], 2L, r, i)
    }
    for (i in R) {
      if (cfg$protracted) {
        x <- if (L$s[i] == k - 1L) cfg$beta[i] * f[i] else cfg$beta[i]
      } else {
        x <- cfg$lac[i] * f[i]
      }
      if (x > 0) add(x, 3L, r, i)
    }
    for (i in R) {
      for (j in seq_len(n)) {
        if (j == i || cfg$D[i, j] == 0 || !exists[j]) next
        if (L$occ[j] == 0L) {
          x <- cfg$D[i, j] * f[j]
        } else if (L$S[i, j] > 0 || (cfg$protracted && L$s[j] > 0L)) {
          x <- cfg$D[i, j]
        } else {
          x <- 0
        }
        if (x > 0) add(x, 4L, r, i, j)
      }
    }
    if (cfg$alpha > 0) {
      nodes <- da_nodes(L, n)
      for (a in nodes) {
        for (b in nodes) {
          if (a < b && L$S[a, b] >= 0 && L$S[a, b] < k) {
            add(cfg$alpha, 5L, r, a, b)
          }
        }
      }
    }
  }
  cs <- numeric(m_ev)
  acc <- 0
  for (e in seq_len(m_ev)) {
    acc <- acc + rate[e]
    cs[e] <- acc
  }
  list(n = m_ev, cs = cs, type = type[seq_len(m_ev)], r = er[seq_len(m_ev)],
       a = ea[seq_len(m_ev)], b = eb[seq_len(m_ev)])
}

#' Make a sampled event of the DAISIE_A simulation happen
#'
#' @param st State list.
#' @param ev Events, see \code{\link{da_events}()}.
#' @param e Index of the event.
#' @param t Current (forward) time.
#' @param cfg Settings.
#' @return Updated state list.
#' @keywords internal
da_apply <- function(st, ev, e, t, cfg) {
  n <- cfg$n
  n1 <- n + 1L
  r <- ev$r[e]
  a <- ev$a[e]
  b <- ev$b[e]
  type <- ev$type[e]
  if (type == 1L) {
    # colonization of island a by mainland species r
    m <- r
    i <- a
    r <- 0L
    for (q in seq_along(st$lin)) {
      if (st$lin[[q]]$type == 1 && st$lin[[q]]$id == m) r <- q
    }
    if (r == 0L || sum(st$lin[[r]]$occ) == 1L && st$lin[[r]]$occ[i] == 1L) {
      if (r > 0L) {
        # recolonization of the only island of the non-endemic lineage: it is
        # replaced by the new colonist (its colonization time is reset)
        old <- st$lin[[r]]
        st <- da_add_node(st, parent = old$node, time = t, kind = 4,
                          side = old$side, clade = old$clade,
                          species = old$id, time_split = t)
      }
      st$n_clades <- st$n_clades + 1L
      st <- da_add_node(st, parent = 0, time = t, kind = 1, side = 0,
                        clade = st$n_clades, species = m, time_split = t)
      L <- da_new_lineage(id = m, anc = m, colt = t, type = 1, ana = 0,
                          node = st$n_gen, side = 0, clade = st$n_clades,
                          n = n)
      L$occ[i] <- 1L
      L$S[i, n1] <- 0
      L$S[n1, i] <- 0
      if (r == 0L) {
        st$lin[[length(st$lin) + 1L]] <- L
      } else {
        st$lin[[r]] <- L
        st <- da_clade_singleton(st, old$clade)
      }
    } else {
      L <- st$lin[[r]]
      if (L$occ[i] == 0L) {
        L <- da_add_island(L, i, n1, cfg)
      } else {
        # recolonization of an island of a widespread non-endemic lineage
        L <- da_set_pair(L, i, n1, 0, -Inf)
        if (cfg$protracted) {
          L$s[i] <- 0L
          L$stau[i] <- -Inf
        }
        L <- da_drop_gaps(L, cfg)
      }
      st$lin[[r]] <- L
    }
  } else if (type == 2L) {
    # local extinction on island a
    L <- st$lin[[r]]
    if (sum(L$occ) == 1L) {
      st <- da_add_node(st, parent = L$node, time = t, kind = 3,
                        side = L$side, clade = L$clade, species = L$id,
                        time_split = t)
      st$lin[[r]] <- NULL
      st <- da_clade_singleton(st, L$clade)
    } else {
      st$lin[[r]] <- da_remove_island(L, a, cfg)
      st <- da_check_split(st, r, t, cfg)
    }
  } else if (type == 3L) {
    if (cfg$protracted) {
      st <- da_within_step(st, r, a, t, cfg)
    } else {
      st <- da_within_speciation(st, r, a, t, cfg)
    }
  } else if (type == 4L) {
    # dispersal from island a to island b
    L <- st$lin[[r]]
    if (L$occ[b] == 0L) {
      L <- da_add_island(L, b, a, cfg)
    } else {
      L <- da_set_pair(L, a, b, 0, -Inf)
      if (cfg$protracted) {
        L$s[b] <- 0L
        L$stau[b] <- -Inf
      }
      L <- da_drop_gaps(L, cfg)
    }
    st$lin[[r]] <- L
  } else {
    # divergence step of pair (a, b)
    L <- st$lin[[r]]
    new_s <- L$S[a, b] + 1
    L <- da_set_pair(L, a, b, new_s, if (new_s == 1) t else L$Tau[a, b])
    st$lin[[r]] <- L
    if (new_s == cfg$k) st <- da_check_split(st, r, t, cfg)
  }
  st
}

#' Within-island speciation (budding, or split of a single-island lineage)
#'
#' @inheritParams da_apply
#' @param r Lineage.
#' @param i Island.
#' @return Updated state list.
#' @keywords internal
da_within_speciation <- function(st, r, i, t, cfg) {
  L <- st$lin[[r]]
  st <- da_add_node(st, parent = L$node, time = t, kind = 2, side = L$side,
                    clade = L$clade, species = NA, time_split = t)
  v <- st$n_gen
  if (sum(L$occ) == 1L && cfg$split_single) {
    # as in DAISIE: two new endemic species replace the parent
    st <- da_two_new_species(st, L, i, v)
    st$lin[[r]] <- st$new_a
    st$lin[[length(st$lin) + 1L]] <- st$new_b
  } else {
    # budding: the parent keeps its identity and range, one new endemic
    # species appears on this island
    L$node <- v
    L$side <- 1
    st$lin[[r]] <- L
    st$maxspecID <- st$maxspecID + 1
    B <- da_new_lineage(id = st$maxspecID, anc = L$anc, colt = L$colt,
                        type = 3, ana = 0, node = v, side = 2,
                        clade = L$clade, n = cfg$n)
    B$occ[i] <- 1L
    st$lin[[length(st$lin) + 1L]] <- B
  }
  st$new_a <- NULL
  st$new_b <- NULL
  st
}

#' Within-island step of the protracted option
#'
#' @description The within-island state of the population on island \code{i}
#'   advances by one; when it reaches \code{k}, the population becomes two new
#'   species on island \code{i}. A lineage on more islands keeps the other
#'   islands (node 1, dated at the step 0 -> 1); the two new species are
#'   separated by node 2 (at the completing step).
#'
#' @inheritParams da_within_speciation
#' @return Updated state list.
#' @keywords internal
da_within_step <- function(st, r, i, t, cfg) {
  L <- st$lin[[r]]
  L$s[i] <- L$s[i] + 1L
  if (L$s[i] == 1L) L$stau[i] <- t
  if (L$s[i] < cfg$k) {
    st$lin[[r]] <- L
    return(st)
  }
  if (sum(L$occ) == 1L) {
    st <- da_add_node(st, parent = L$node, time = t, kind = 2,
                      side = L$side, clade = L$clade, species = NA,
                      time_split = t)
    st <- da_two_new_species(st, L, i, st$n_gen)
    st$lin[[r]] <- st$new_a
    st$lin[[length(st$lin) + 1L]] <- st$new_b
  } else {
    t_start <- max(L$stau[i], st$gen[L$node, "time"])
    st <- da_add_node(st, parent = L$node, time = t_start, kind = 2,
                      side = L$side, clade = L$clade, species = NA,
                      time_split = t)
    v1 <- st$n_gen
    st <- da_add_node(st, parent = v1, time = t, kind = 2, side = 2,
                      clade = L$clade, species = NA, time_split = t)
    st <- da_two_new_species(st, L, i, st$n_gen)
    L$node <- v1
    L$side <- 1
    st$lin[[r]] <- da_remove_island(L, i, cfg)
    st$lin[[length(st$lin) + 1L]] <- st$new_a
    st$lin[[length(st$lin) + 1L]] <- st$new_b
    st <- da_check_split(st, r, t, cfg)
  }
  st$new_a <- NULL
  st$new_b <- NULL
  st
}

#' Two new endemic species on one island, descending from genealogy node v
#'
#' @param st State list.
#' @param L Parent lineage.
#' @param i Island.
#' @param v Genealogy node.
#' @return State list with the new lineages in \code{new_a} and \code{new_b}.
#' @keywords internal
da_two_new_species <- function(st, L, i, v) {
  n <- length(L$occ)
  for (side in 1:2) {
    st$maxspecID <- st$maxspecID + 1
    x <- da_new_lineage(id = st$maxspecID, anc = L$anc, colt = L$colt,
                        type = 3, ana = 0, node = v, side = side,
                        clade = L$clade, n = n)
    x$occ[i] <- 1L
    if (side == 1) st$new_a <- x else st$new_b <- x
  }
  st
}

#' Submergence of an island in the simulation
#'
#' @param st State list.
#' @param i Island.
#' @param t Current (forward) time.
#' @param cfg Settings.
#' @return Updated state list.
#' @keywords internal
da_submerge <- function(st, i, t, cfg) {
  r <- 1L
  while (r <= length(st$lin)) {
    L <- st$lin[[r]]
    if (L$occ[i] == 0L) {
      r <- r + 1L
      next
    }
    if (sum(L$occ) == 1L) {
      st <- da_add_node(st, parent = L$node, time = t, kind = 3,
                        side = L$side, clade = L$clade, species = L$id,
                        time_split = t)
      st$lin[[r]] <- NULL
      st <- da_clade_singleton(st, L$clade)
      next
    }
    st$lin[[r]] <- da_remove_island(L, i, cfg)
    st <- da_check_split(st, r, t, cfg)
    r <- r + 1L
  }
  st
}

#' New lineage of the DAISIE_A simulation
#'
#' @param id,anc,colt,type,ana,node,side,clade Lineage attributes.
#' @param n Number of islands.
#' @return Lineage list with an empty range and no tracked pairs.
#' @keywords internal
da_new_lineage <- function(id, anc, colt, type, ana, node, side, clade, n) {
  list(id = id, anc = anc, colt = colt, type = type, ana = ana, node = node,
       side = side, clade = clade, occ = integer(n),
       S = matrix(-1, n + 1L, n + 1L), Tau = matrix(-Inf, n + 1L, n + 1L),
       s = integer(n), stau = rep(-Inf, n))
}

#' Nodes (islands of the range, and the mainland for a non-endemic)
#'
#' @param L Lineage.
#' @param n Number of islands.
#' @return Integer vector of node indices (mainland = n + 1).
#' @keywords internal
da_nodes <- function(L, n) {
  v <- which(L$occ == 1L)
  if (L$type == 1) v <- c(v, n + 1L)
  v
}

#' Set the state of a pair
#'
#' @param L Lineage.
#' @param a,b Nodes.
#' @param s State (-1: not tracked).
#' @param tau Time of the step 0 -> 1.
#' @return Updated lineage.
#' @keywords internal
da_set_pair <- function(L, a, b, s, tau) {
  L$S[a, b] <- s
  L$S[b, a] <- s
  L$Tau[a, b] <- tau
  L$Tau[b, a] <- tau
  L
}

#' Connected components (breadth-first, in increasing node order)
#'
#' @param nodes Nodes, in increasing order.
#' @param E Logical matrix of edges.
#' @return Component label of each node (0 for nodes not in \code{nodes});
#'   components are numbered in order of their smallest node.
#' @keywords internal
da_comp <- function(nodes, E) {
  comp <- integer(nrow(E))
  nc <- 0L
  for (v in nodes) {
    if (comp[v] != 0L) next
    nc <- nc + 1L
    comp[v] <- nc
    queue <- v
    h <- 1L
    while (h <= length(queue)) {
      u <- queue[h]
      h <- h + 1L
      for (w in nodes) {
        if (comp[w] == 0L && E[u, w]) {
          comp[w] <- nc
          queue <- c(queue, w)
        }
      }
    }
  }
  comp
}

#' Order of pair labels (state, then time of the step 0 -> 1)
#'
#' @param s1,t1,s2,t2 States and times.
#' @return \code{TRUE} if label 1 is smaller than label 2.
#' @keywords internal
da_less <- function(s1, t1, s2, t2) {
  s1 < s2 || (s1 == s2 && t1 < t2)
}

#' Path values of all pairs of nodes of a lineage
#'
#' @description The path value of two nodes is the minimum, over all paths of
#'   tracked pairs between them, of the largest label (state, then time) on
#'   the path (Floyd-Warshall on the bottleneck semiring).
#'
#' @param L Lineage.
#' @param nodes Its nodes.
#' @return List with matrices \code{s} and \code{t}.
#' @keywords internal
da_bottleneck <- function(L, nodes) {
  n1 <- nrow(L$S)
  big <- .Machine$integer.max
  bs <- matrix(big, n1, n1)
  bt <- matrix(Inf, n1, n1)
  for (a in nodes) {
    for (b in nodes) {
      if (a == b) {
        bs[a, b] <- 0
        bt[a, b] <- -Inf
      } else if (L$S[a, b] >= 0) {
        bs[a, b] <- L$S[a, b]
        bt[a, b] <- L$Tau[a, b]
      }
    }
  }
  for (m in nodes) {
    for (a in nodes) {
      if (bs[a, m] == big) next
      for (b in nodes) {
        if (bs[m, b] == big) next
        if (da_less(bs[a, m], bt[a, m], bs[m, b], bt[m, b])) {
          cs <- bs[m, b]
          ct <- bt[m, b]
        } else {
          cs <- bs[a, m]
          ct <- bt[a, m]
        }
        if (da_less(cs, ct, bs[a, b], bt[a, b])) {
          bs[a, b] <- cs
          bt[a, b] <- ct
        }
      }
    }
  }
  list(s = bs, t = bt)
}

#' Add an island to the range of a lineage (copy rule)
#'
#' @description The new population on island \code{j} is a copy of its
#'   source \code{src} (an island, or the mainland \code{n + 1}): its state
#'   with every node \code{l} of the range adjacent to \code{j} is that of the
#'   pair (\code{src}, \code{l}), or its path value if that pair is not
#'   tracked, and 0 with the source itself. Gap pairs that are bridged by the
#'   new island are dropped.
#'
#' @param L Lineage.
#' @param j New island.
#' @param src Source node.
#' @param cfg Settings.
#' @return Updated lineage.
#' @keywords internal
da_add_island <- function(L, j, src, cfg) {
  nodes <- da_nodes(L, cfg$n)
  B <- da_bottleneck(L, nodes)
  for (l in nodes) {
    if (!cfg$A1[j, l]) next
    if (l == src) {
      L <- da_set_pair(L, j, l, 0, -Inf)
    } else if (L$S[src, l] >= 0) {
      L <- da_set_pair(L, j, l, L$S[src, l], L$Tau[src, l])
    } else {
      L <- da_set_pair(L, j, l, B$s[src, l], B$t[src, l])
    }
  }
  L$occ[j] <- 1L
  L$s[j] <- 0L
  L$stau[j] <- -Inf
  da_drop_gaps(L, cfg)
}

#' Drop gap pairs between nodes that are connected by adjacent pairs below k
#'
#' @description A gap pair is removed when its two nodes are connected again
#'   through adjacent pairs that are all below \code{k} (e.g. after the gap
#'   has been refilled, or after a reset), so that removing it can never
#'   disconnect the network of pairs below \code{k}. Called after range
#'   expansions and resets.
#'
#' @param L Lineage.
#' @param cfg Settings.
#' @return Updated lineage.
#' @keywords internal
da_drop_gaps <- function(L, cfg) {
  nodes <- da_nodes(L, cfg$n)
  comp <- da_comp(nodes, cfg$A1 & (L$S >= 0) & (L$S < cfg$k))
  for (a in nodes) {
    for (b in nodes) {
      if (a < b && L$S[a, b] >= 0 && !cfg$A1[a, b] && comp[a] == comp[b]) {
        L <- da_set_pair(L, a, b, -1, -Inf)
      }
    }
  }
  L
}

#' Remove an island from the range of a lineage (gap rule)
#'
#' @description Drops the pairs with island \code{l}; every two of its
#'   neighbors (nodes with a tracked pair with \code{l}) that are no longer
#'   connected by tracked pairs get a gap pair with the larger of their two
#'   labels with \code{l}.
#'
#' @param L Lineage (with at least two islands).
#' @param l Island.
#' @param cfg Settings.
#' @return Updated lineage.
#' @keywords internal
da_remove_island <- function(L, l, cfg) {
  nodes <- da_nodes(L, cfg$n)
  nb <- nodes[nodes != l & L$S[l, nodes] >= 0]
  ls <- L$S[l, ]
  lt <- L$Tau[l, ]
  L$S[l, ] <- -1
  L$S[, l] <- -1
  L$Tau[l, ] <- -Inf
  L$Tau[, l] <- -Inf
  L$occ[l] <- 0L
  L$s[l] <- 0L
  L$stau[l] <- -Inf
  comp <- da_comp(nodes[nodes != l], L$S >= 0)
  m <- length(nb)
  if (m >= 2) {
    for (x in 1:(m - 1)) {
      for (y in (x + 1):m) {
        a <- nb[x]
        b <- nb[y]
        if (comp[a] == comp[b]) next
        if (da_less(ls[a], lt[a], ls[b], lt[b])) {
          L <- da_set_pair(L, a, b, ls[b], lt[b])
        } else {
          L <- da_set_pair(L, a, b, ls[a], lt[a])
        }
      }
    }
  }
  L
}

#' Restrict a lineage to a set of nodes
#'
#' @param L Lineage.
#' @param P Nodes to keep.
#' @param cfg Settings.
#' @return Lineage with only the islands of \code{P} and the pairs within
#'   \code{P}.
#' @keywords internal
da_restrict <- function(L, P, cfg) {
  keep <- logical(cfg$n + 1L)
  keep[P] <- TRUE
  L$S[!keep, ] <- -1
  L$S[, !keep] <- -1
  L$Tau[!keep, ] <- -Inf
  L$Tau[, !keep] <- -Inf
  isl <- keep[seq_len(cfg$n)]
  L$occ[!isl] <- 0L
  L$s[!isl] <- 0L
  L$stau[!isl] <- -Inf
  L
}

#' Latest step 0 -> 1 among the tracked pairs between two sets of nodes
#'
#' @param L Lineage.
#' @param P,Q Sets of nodes.
#' @return Time (\code{-Inf} if there is no tracked pair).
#' @keywords internal
da_cut_tau <- function(L, P, Q) {
  m <- -Inf
  for (a in P) {
    for (b in Q) {
      if (L$S[a, b] >= 0 && L$Tau[a, b] > m) m <- L$Tau[a, b]
    }
  }
  m
}

#' Random permutation (Fisher-Yates, with the same random numbers as C++)
#'
#' @param m Length.
#' @return Permutation of \code{1:m}.
#' @keywords internal
da_perm <- function(m) {
  x <- seq_len(m)
  if (m >= 2) {
    for (j in m:2) {
      idx <- floor(stats::runif(1) * j) + 1
      tmp <- x[j]
      x[j] <- x[idx]
      x[idx] <- tmp
    }
  }
  x
}

#' Split a lineage if its network of pairs below k is disconnected
#'
#' @param st State list.
#' @param r Lineage.
#' @param t Current (forward) time.
#' @param cfg Settings.
#' @return Updated state list.
#' @keywords internal
da_check_split <- function(st, r, t, cfg) {
  L <- st$lin[[r]]
  nodes <- da_nodes(L, cfg$n)
  comp <- da_comp(nodes, (L$S >= 0) & (L$S < cfg$k))
  nc <- max(comp[nodes])
  if (nc == 1L) return(st)
  pieces <- lapply(seq_len(nc), function(x) nodes[comp[nodes] == x])
  da_split(st, r, pieces, t, cfg)
}

#' Split a lineage into its pieces
#'
#' @description An endemic lineage becomes one new species per piece. For a
#'   non-endemic lineage the piece with the mainland stays non-endemic (if it
#'   has islands) and the other pieces become new endemic species; if the
#'   mainland piece has no islands, the lineage becomes endemic (anagenesis),
#'   as one species or, with several island pieces, as one new species per
#'   piece. Three or more species from one event are resolved into binary
#'   nodes of zero length, in random order.
#'
#' @param st State list.
#' @param r Lineage.
#' @param pieces List of node sets.
#' @param t Current (forward) time.
#' @param cfg Settings.
#' @return Updated state list.
#' @keywords internal
da_split <- function(st, r, pieces, t, cfg) {
  n1 <- cfg$n + 1L
  L <- st$lin[[r]]
  if (L$type == 1) {
    pm <- which(vapply(pieces, function(p) n1 %in% p, logical(1)))
    others <- pieces[-pm]
    if (length(pieces[[pm]]) > 1L) {
      ord <- if (length(others) >= 2L) others[da_perm(length(others))] else others
      return(da_peel(st, r, ord, length(ord), pieces[[pm]], t, cfg))
    }
    # anagenesis: the mainland separates from all islands
    isl <- unlist(others)
    tau_m <- da_cut_tau(L, n1, isl)
    st$lin[[r]] <- da_restrict(L, isl, cfg)
    if (cfg$protracted) st <- da_set_clade_colt(st, L$clade, tau_m)
    if (length(others) == 1L) {
      st$maxspecID <- st$maxspecID + 1
      st$lin[[r]]$id <- st$maxspecID
      st$lin[[r]]$type <- 2
      st$lin[[r]]$ana <- 1
      return(st)
    }
    pieces <- others
  }
  m <- length(pieces)
  ord <- if (m >= 3L) pieces[da_perm(m)] else pieces
  st <- da_peel(st, r, ord, m - 1L, NULL, t, cfg)
  st$maxspecID <- st$maxspecID + 1
  st$lin[[r]]$id <- st$maxspecID
  st$lin[[r]]$type <- 3
  st$lin[[r]]$ana <- 0
  st
}

#' Peel pieces off a lineage, one binary node per piece
#'
#' @param st State list.
#' @param r Lineage.
#' @param ord Pieces in order.
#' @param n_peel Number of pieces to peel off (the first \code{n_peel}).
#' @param extra Nodes that stay with the lineage besides the remaining pieces.
#' @param t Current (forward) time.
#' @param cfg Settings.
#' @return Updated state list.
#' @keywords internal
da_peel <- function(st, r, ord, n_peel, extra, t, cfg) {
  for (j in seq_len(n_peel)) {
    P <- ord[[j]]
    Q <- c(if (j < length(ord)) unlist(ord[(j + 1):length(ord)]), extra)
    L <- st$lin[[r]]
    t_node <- if (cfg$protracted) {
      max(da_cut_tau(L, P, Q), st$gen[L$node, "time"])
    } else {
      t
    }
    st <- da_add_node(st, parent = L$node, time = t_node, kind = 2,
                      side = L$side, clade = L$clade, species = NA,
                      time_split = t)
    v <- st$n_gen
    B <- da_restrict(L, P, cfg)
    st$maxspecID <- st$maxspecID + 1
    B$id <- st$maxspecID
    B$type <- 3
    B$ana <- 0
    B$node <- v
    B$side <- 2
    L <- da_restrict(L, Q, cfg)
    L$node <- v
    L$side <- 1
    st$lin[[r]] <- L
    st$lin[[length(st$lin) + 1L]] <- B
  }
  st
}

#' Date the colonization of a clade at the start of its divergence from the
#'   mainland (protracted option)
#'
#' @param st State list.
#' @param clade Clade.
#' @param tau Time.
#' @return Updated state list.
#' @keywords internal
da_set_clade_colt <- function(st, clade, tau) {
  for (q in seq_along(st$lin)) {
    if (st$lin[[q]]$clade == clade) st$lin[[q]]$colt <- tau
  }
  root <- which(st$gen[seq_len(st$n_gen), "kind"] == 1 &
                  st$gen[seq_len(st$n_gen), "clade"] == clade)
  st$gen[root, "time"] <- tau
  st
}

#' Add a node to the genealogy
#'
#' @param st State list.
#' @param parent,time,kind,side,clade,species,time_split Node attributes.
#' @return Updated state; the new node has ID \code{st$n_gen}.
#' @keywords internal
da_add_node <- function(st, parent, time, kind, side, clade, species,
                        time_split) {
  st$n_gen <- st$n_gen + 1L
  if (st$n_gen > nrow(st$gen)) {
    st$gen <- rbind(st$gen, matrix(NA_real_, nrow = nrow(st$gen),
                                   ncol = ncol(st$gen)))
  }
  st$gen[st$n_gen, ] <- c(st$n_gen, parent, time, kind, side, clade,
                          species, time_split)
  st
}

#' Relabel the last surviving cladogenetic species of a clade
#'
#' @description As in DAISIE, if only one species of a clade is left and it
#'   is of cladogenetic origin, it becomes an anagenetic species with
#'   anagenetic origin "Clado_extinct".
#'
#' @param st State list.
#' @param clade Clade (colonization event) ID.
#' @return Updated state list.
#' @keywords internal
da_clade_singleton <- function(st, clade) {
  rows <- which(vapply(st$lin, function(x) x$clade == clade, logical(1)))
  if (length(rows) == 1 && st$lin[[rows]]$type == 3) {
    st$lin[[rows]]$type <- 2
    st$lin[[rows]]$ana <- 2
  }
  st
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

#' Column names of the per-island species-through-time table
#'
#' @param n_islands Number of islands.
#' @return Character vector.
#' @keywords internal
DAISIE_A_stt_island_names <- function(n_islands) {
  c("Time", paste0(rep(c("nI_", "nA_", "nC_", "nSIE_"), n_islands),
                   rep(seq_len(n_islands), each = 4)))
}

#' Initialize the species-through-time tables
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
  colnames(isl) <- DAISIE_A_stt_island_names(n_islands)
  isl[1, ] <- c(total_time, rep(0, 4 * n_islands))
  list(table = tab, islands = isl, n = 1L)
}

#' Add a row to the species-through-time tables if the counts changed
#'
#' @param stt List of STT tables.
#' @param st State list.
#' @param time_bp Time before present.
#' @param cfg Settings.
#' @return Updated \code{stt}.
#' @keywords internal
da_stt_row <- function(stt, st, time_bp, cfg) {
  n <- cfg$n
  arch <- numeric(3)
  isl <- numeric(4 * n)
  for (L in st$lin) {
    arch[L$type] <- arch[L$type] + 1
    single <- L$type != 1 && sum(L$occ) == 1L
    for (i in which(L$occ == 1L)) {
      isl[4 * (i - 1) + L$type] <- isl[4 * (i - 1) + L$type] + 1
      if (single) isl[4 * (i - 1) + 4] <- isl[4 * (i - 1) + 4] + 1
    }
  }
  if (all(stt$table[stt$n, 2:4] == arch) &&
      all(stt$islands[stt$n, -1] == isl)) {
    return(stt)
  }
  stt$n <- stt$n + 1L
  if (stt$n > nrow(stt$table)) {
    stt$table <- rbind(stt$table, matrix(NA_real_, nrow = nrow(stt$table),
                                         ncol = ncol(stt$table)))
    stt$islands <- rbind(stt$islands,
                         matrix(NA_real_, nrow = nrow(stt$islands),
                                ncol = ncol(stt$islands)))
  }
  stt$table[stt$n, ] <- c(time_bp, arch)
  stt$islands[stt$n, ] <- c(time_bp, isl)
  stt
}

#' Check the invariants of the state of the DAISIE_A simulation
#'
#' @description Used in testing (\code{options(DAISIE_A.debug = TRUE)}):
#'   stops if a lineage has pairs outside its range, an untracked adjacent
#'   pair, a gap pair within a part connected by adjacent pairs, a state
#'   outside \code{0..k}, a disconnected network of pairs below \code{k}, or if
#'   a mainland species has more than one non-endemic lineage.
#'
#' @param st State list.
#' @param cfg Settings.
#' @return \code{TRUE}, invisibly.
#' @keywords internal
da_check_invariants <- function(st, cfg) {
  n <- cfg$n
  n1 <- n + 1L
  ne <- integer(0)
  for (r in seq_along(st$lin)) {
    L <- st$lin[[r]]
    fail <- function(msg) stop("lineage ", r, ": ", msg)
    if (sum(L$occ) < 1) fail("empty range")
    if (L$type == 1) ne <- c(ne, L$id)
    nodes <- da_nodes(L, n)
    innode <- logical(n1)
    innode[nodes] <- TRUE
    if (any(L$S[!innode, ] >= 0) || any(L$S[, !innode] >= 0)) {
      fail("pair outside the range")
    }
    if (!isSymmetric(L$S) || !isSymmetric(L$Tau)) fail("asymmetric")
    if (any(diag(L$S) >= 0)) fail("diagonal pair")
    tr <- L$S >= 0
    if (any(L$S[tr] > cfg$k)) fail("state above k")
    if (any((L$S[tr] == 0) != (L$Tau[tr] == -Inf))) fail("tau inconsistent")
    if (any(L$Tau[!tr] != -Inf)) fail("tau of untracked pair")
    adj <- cfg$A1[nodes, nodes, drop = FALSE]
    if (any(adj & !tr[nodes, nodes, drop = FALSE])) fail("untracked adjacent pair")
    comp_adj <- da_comp(nodes, cfg$A1 & tr & L$S < cfg$k)
    for (a in nodes) for (b in nodes) {
      if (a < b && tr[a, b] && !cfg$A1[a, b] && comp_adj[a] == comp_adj[b]) {
        fail("gap pair within a part connected by adjacent pairs below k")
      }
    }
    comp <- da_comp(nodes, tr & L$S < cfg$k)
    if (max(comp[nodes]) != 1) fail("network of pairs below k disconnected")
    if (cfg$protracted) {
      if (any(L$s[L$occ == 1L] >= cfg$k)) fail("within-island state at k")
      if (any(L$s[L$occ == 0L] != 0)) fail("within-island state outside range")
    }
    if (st$gen[L$node, "time"] > st$gen[L$node, "time_split"] + 1e-12 &&
        st$gen[L$node, "kind"] != 1) {
      fail("node dated after its split")
    }
  }
  if (any(duplicated(ne))) stop("two non-endemic lineages of one mainland species")
  invisible(TRUE)
}
