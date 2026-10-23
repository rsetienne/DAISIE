#' Simulate archipelagos under the archipelago DAISIE model (DAISIE_A)
#'
#' @description
#' Simulates the assembly of an archipelago of \eqn{n} islands by colonisation
#' from the mainland, dispersal between islands, local and global extinction,
#' anagenesis, within-island speciation and between-island (vicariant)
#' speciation, with clade-specific diversity dependence on each island. Every
#' lineage has a range, i.e. the set of islands on which it occurs, coded as a
#' 0/1 vector \eqn{(b_1, \ldots, b_n)} (one of \eqn{2^n - 1} states). The
#' model is described in detail in the accompanying report; the rates are set
#' with \code{\link{DAISIE_A_create_pars}()}.
#'
#' Recolonisation: the colonisation time of a non-endemic lineage is reset
#' only when the mainland species recolonises the single island on which the
#' lineage occurs; colonisation of another island expands the range without
#' reset, and recolonisation of an island in a widespread range has no effect.
#'
#' Each of the \code{M} mainland species is simulated independently (as in
#' the clade-specific (CS) version of DAISIE). With a single island
#' (\code{n_islands = 1}) the model is identical to DAISIE and, for the same
#' random seed, this function returns the same archipelago-level output as
#' \code{\link{DAISIE_sim_cr}()}.
#'
#' @param time Length of the simulation (age of the oldest island).
#' @param M Number of mainland species (size of the mainland pool).
#' @param pars A \code{DAISIE_A_pars} object, see
#'   \code{\link{DAISIE_A_create_pars}()}.
#' @param replicates Number of archipelago replicates.
#' @param sample_freq Number of time points at which the species-through-time
#'   tables are sampled; \code{Inf} keeps all events.
#' @param cond Condition on the number of colonising lineages present on the
#'   archipelago at the end of the simulation (as in DAISIE; 0 = none).
#' @param plot_sims Logical; if \code{TRUE}, plots the archipelago-level
#'   species-through-time plot with \code{\link{DAISIE_plot_sims}()}.
#' @param verbose Logical; print progress messages.
#'
#' @return A list with one element per replicate. Each replicate is a list in
#'   the standard DAISIE format, in which the archipelago is treated as one
#'   island, extended with archipelago information:
#'   \itemize{
#'     \item \code{[[1]]}: \code{island_age}, \code{not_present},
#'       \code{stt_all} (as in DAISIE, for the whole archipelago),
#'       \code{n_islands}, \code{island_names}, and \code{stt_islands}: a list
#'       with one matrix per island with columns \code{Time}, \code{nI}
#'       (non-endemic species), \code{nA} (anagenetic endemics), \code{nC}
#'       (cladogenetic endemics), \code{nSIE} (single-island endemics) and
#'       \code{present} (number of mainland species with descendants on the
#'       island).
#'     \item \code{[[2]]}, \code{[[3]]}, ...: one element per mainland species
#'       with descendants on the archipelago, with the DAISIE elements
#'       \code{branching_times}, \code{stac}, \code{missing_species} (and
#'       \code{all_colonisations} after recolonisation) and the archipelago
#'       elements
#'       \itemize{
#'         \item \code{species}: data frame with one row per extant species:
#'           \code{species} ID, \code{mainland_ancestor},
#'           \code{colonisation_time} (before present), \code{type}
#'           (\code{"I"} non-endemic, \code{"A"} anagenetic or relict
#'           endemic, \code{"C"} cladogenetic endemic),
#'           \code{anagenetic_origin}, \code{range} (e.g. \code{"101"} for a
#'           species on islands 1 and 3), \code{range_code}
#'           (\eqn{\sum_i b_i 2^{i-1}}), \code{n_islands} and \code{clade}
#'           (colonisation event).
#'         \item \code{clades}: list with, for every colonisation event with
#'           extant descendants: \code{colonisation_time},
#'           \code{crown_age}, \code{branching_times} (colonisation time
#'           followed by the node ages of the reconstructed tree),
#'           \code{newick} (reconstructed tree with the stem as root edge;
#'           tip labels are species IDs), \code{nodes} (the same tree as a
#'           node table with columns \code{node}, \code{parent} (0 for the
#'           crown), \code{age} and \code{species}, used by
#'           \code{\link{DAISIE_A_sim_to_datalist}()}), \code{tips} (data frame with type
#'           and range of each tip) and \code{mixed} (\code{TRUE} if the clade
#'           contains both a non-endemic species and endemics budded off from
#'           it).
#'       }
#'   }
#'
#' @author Rampal S. Etienne
#' @family DAISIE_A
#' @export
#' @examples
#' pars <- DAISIE_A_create_pars(
#'   n_islands = 3, lac = 0.5, mu = 0.2, K = Inf, gam = 0.02, laa = 0.5,
#'   lav = 0.2, disp = 0.3, connectivity = "line"
#' )
#' set.seed(42)
#' sims <- DAISIE_A_sim(time = 2, M = 50, pars = pars, replicates = 2,
#'   verbose = FALSE)
#' sims[[1]][[1]]$stt_islands$island_1
DAISIE_A_sim <- function(time,
                         M,
                         pars,
                         replicates,
                         sample_freq = 25,
                         cond = 0,
                         plot_sims = FALSE,
                         verbose = TRUE) {
  total_time <- time
  if (!DAISIE_A_are_pars(pars)) {
    stop("'pars' must be created with DAISIE_A_create_pars().")
  }
  testit::assert(
    !is.na(sample_freq) && !is.null(sample_freq) && sample_freq >= 1
  )
  split_cache <- new.env()
  island_replicates <- list()
  for (rep in 1:replicates) {
    full_list <- list()
    if (cond == 0) {
      number_present <- -1
    } else {
      number_present <- 0
    }
    while (number_present < cond) {
      for (m_spec in 1:M) {
        full_list[[m_spec]] <- DAISIE_A_sim_core(
          time = total_time,
          mainland_n = 1,
          pars = pars,
          split_cache = split_cache
        )
      }
      stac_vec <- vapply(full_list, function(x) x$stac, numeric(1))
      present <- which(stac_vec != 0)
      number_present <- length(present)
    }
    if (verbose == TRUE) {
      message("Archipelago replicate ", rep)
    }
    island_replicates[[rep]] <- DAISIE_A_format(
      full_list = full_list,
      total_time = total_time,
      M = M,
      sample_freq = sample_freq,
      pars = pars,
      verbose = verbose
    )
  }
  if (plot_sims == TRUE) {
    DAISIE_plot_sims(
      island_replicates = island_replicates,
      sample_freq = sample_freq
    )
  }
  return(island_replicates)
}

#' Format the output of the archipelago simulation of one replicate
#'
#' @param full_list List with the output of \code{\link{DAISIE_A_sim_core}()}
#'   for each mainland species.
#' @param total_time Simulation time.
#' @param M Number of mainland species.
#' @param sample_freq Number of sampling points of the STT tables.
#' @param pars A \code{DAISIE_A_pars} object.
#' @param verbose Logical.
#'
#' @return One replicate, see \code{\link{DAISIE_A_sim}()}.
#' @keywords internal
DAISIE_A_format <- function(full_list,
                            total_time,
                            M,
                            sample_freq,
                            pars,
                            verbose = FALSE) {
  daisie_fields <- c("stt_table", "branching_times", "stac",
                     "missing_species", "all_colonisations")
  daisie_list <- lapply(full_list, function(x) x[intersect(names(x),
                                                           daisie_fields)])
  formatted <- DAISIE_format_CS(
    island_replicates = list(daisie_list),
    time = total_time,
    M = M,
    sample_freq = sample_freq,
    verbose = FALSE
  )[[1]]
  stac_vec <- vapply(full_list, function(x) x$stac, numeric(1))
  present <- which(stac_vec != 0)
  for (i in seq_along(present)) {
    x <- full_list[[present[i]]]
    formatted[[1 + i]]$species <- x$species
    formatted[[1 + i]]$clades <- x$clades
    if (!is.null(x$mixed_clade)) formatted[[1 + i]]$mixed_clade <- TRUE
  }
  formatted[[1]]$n_islands <- pars$n_islands
  formatted[[1]]$island_names <- pars$island_names
  formatted[[1]]$stt_islands <- DAISIE_A_format_stt_islands(
    full_list = full_list,
    total_time = total_time,
    sample_freq = sample_freq,
    pars = pars
  )
  return(formatted)
}

#' Per-island species-through-time tables summed over clades
#'
#' @inheritParams DAISIE_A_format
#' @return A named list with one matrix per island, with columns
#'   \code{Time}, \code{nI}, \code{nA}, \code{nC}, \code{nSIE} and
#'   \code{present}.
#' @keywords internal
DAISIE_A_format_stt_islands <- function(full_list,
                                        total_time,
                                        sample_freq,
                                        pars) {
  n <- pars$n_islands
  if (is.infinite(sample_freq)) {
    times <- sort(unique(unlist(lapply(full_list, function(x) {
      x$stt_islands[, "Time"]
    }))), decreasing = TRUE)
  } else {
    times <- rev(seq(from = 0, to = total_time, length.out = sample_freq + 1))
  }
  sums <- matrix(0, nrow = length(times), ncol = 4 * n)
  present <- matrix(0, nrow = length(times), ncol = n)
  for (x in full_list) {
    stt <- x$stt_islands
    row_index <- findInterval(-times, -stt[, "Time"])
    vals <- stt[row_index, -1, drop = FALSE]
    sums <- sums + vals
    for (i in seq_len(n)) {
      present[, i] <- present[, i] +
        as.numeric(rowSums(vals[, 4 * (i - 1) + 1:3, drop = FALSE]) > 0)
    }
  }
  out <- list()
  for (i in seq_len(n)) {
    m <- cbind(times, sums[, 4 * (i - 1) + 1:4, drop = FALSE], present[, i])
    colnames(m) <- c("Time", "nI", "nA", "nC", "nSIE", "present")
    out[[pars$island_names[i]]] <- m
  }
  return(out)
}
