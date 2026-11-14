#' Simulate trait-dependent islands, one clade per mainland species
#'
#' @description
#' Simulates island data sets under a trait-dependent model, with optional
#' hidden states. Species colonise from the mainland, go extinct, speciate on
#' the island (cladogenesis and anagenesis) and change state; all rates can
#' depend on the state of the species. Each mainland species is simulated as
#' its own clade (clade-specific, "CS"), and the clades are then put together
#' into one data set per replicate, in the same form as a real island data set.
#'
#' @details
#' \strong{States.} With \code{num_observed_states} observed states and
#' \code{num_hidden_states} hidden states there are
#' \code{n = num_observed_states * num_hidden_states} states. They are numbered
#' 1 to n, observed state first, then hidden state. With 2 observed and 2 hidden
#' states: 1 = 0A, 2 = 0B, 3 = 1A, 4 = 1B. Without hidden states
#' (\code{num_hidden_states = 1}): 1 = state 0, 2 = state 1, and so on.
#' All parameters ending in a number (\code{immig_rate1}, \code{K2}, ...) and
#' the elements of \code{mainland} use this numbering.
#'
#' \strong{How the simulation runs.} For each mainland species, the function
#' runs \code{TRAISIE_sim_core()} on an island where only that species can
#' colonise, with the state of its group in \code{mainland}. Each island starts
#' empty at \code{total_time} million years ago. At each step all rates are
#' computed, the waiting time to the next event is drawn from an exponential
#' distribution with the sum of all rates, and one event is chosen with
#' probability proportional to its rate, until the present. For each state X:
#' \itemize{
#'   \item colonisation: \code{immig_rateX} if the mainland species is in
#'     state X, 0 otherwise. If the species colonises again while it is still
#'     an unchanged immigrant on the island, its colonisation time is reset;
#'   \item extinction: \code{ext_rateX} times the number of island species in X;
#'   \item anagenesis: \code{ana_rateX} times the number of immigrants in X
#'     (the immigrant becomes an endemic species);
#'   \item cladogenesis: \code{clado_rateX} times the number of island species
#'     in X; both new species keep state X;
#'   \item change of state from X to Y: the transition rate from X to Y
#'     (see \code{trans_rate} below) times the number of island species in X.
#' }
#' The clades are then formatted with \code{DAISIE_format_CS()}.
#'
#' @param total_time Island age in million years. The simulation runs from
#'   \code{total_time} to the present (0).
#' @param mainland List with one whole number per state (length n, in the
#'   state order above): the number of mainland species in each state, for
#'   example \code{list(M1 = 100, M2 = 150)} for two states. Each mainland
#'   species gives one clade, so the total is the number of clades simulated.
#'   Use whole numbers (for example with \code{round()}). The names are not used.
#' @param trait_pars Named list of model parameters, built for n states:
#'   \describe{
#'     \item{\code{immig_rateX}}{Colonisation rate per mainland species in state X.}
#'     \item{\code{ext_rateX}}{Extinction rate per island species in state X.}
#'     \item{\code{ana_rateX}}{Anagenesis rate per immigrant species in state X.}
#'     \item{\code{clado_rateX}}{Cladogenesis rate per island species in state X.}
#'     \item{\code{trans_rate1} to \code{trans_rate(n^2)}}{Rates of change of
#'       state, one number each, given row by row: \code{trans_rate1} to
#'       \code{trans_rate(n)} are the changes from state 1 to states 1 to n,
#'       the next n are from state 2, and so on. The change from state i to
#'       state j is \code{trans_rate((i - 1) * n + j)}. The changes from a state
#'       to itself (\code{trans_rate1}, \code{trans_rate(n + 2)}, ...) must be 0.}
#'     \item{\code{KX}}{Carrying capacity for state X. Use \code{Inf} for no
#'       limit.}
#'     \item{\code{p}}{n x n matrix of 0 and 1. \code{p[i, j] = 1} means that
#'       when an immigrant changes from state i to state j it also becomes an
#'       endemic species (anagenesis at the change of state); \code{0} means it
#'       stays an immigrant. Species that are already endemic are not affected.
#'       The diagonal must be 0. Use \code{matrix(0, n, n)} for no anagenesis at
#'       changes of state.}
#'   }
#' @param replicates Number of island data sets to simulate.
#' @param sample_freq Number of time points at which the number of species
#'   on the whole island (\code{stt_all}) is recorded. Use \code{Inf} to keep
#'   every event. Default 100.
#' @param cond Minimum number of clades that must be present on the island
#'   at the present. A replicate with fewer clades is simulated again. Default
#'   0 (keep every replicate).
#' @param verbose If \code{TRUE}, print a message after each replicate.
#' @param files_to_write Leave at 0 (default). Values above 0 are not
#'   working: the function would try to load files that it never saves.
#' @param num_observed_states Number of observed states (at least 1).
#' @param num_hidden_states Number of hidden states; 1 if hidden states are
#'   not used.
#'
#' @returns
#' A list with one element per replicate. Each element is an island data set
#' in the same form as a real one:
#' \itemize{
#'   \item element 1: \code{island_age}, \code{not_present} (number of mainland
#'     species with no descendants on the island at the present) and
#'     \code{stt_all} (number of immigrant, anagenetic and cladogenetic species
#'     on the whole island through time, columns \code{Time}, \code{nI},
#'     \code{nA}, \code{nC}, and \code{present}, the number of clades with at
#'     least one species);
#'   \item elements 2, 3, ...: one per clade present at the present, with
#'     \code{branching_times}, \code{stac}, \code{missing_species},
#'     \code{traits}, \code{root_state}, \code{sampling_fraction} and
#'     \code{phylogeny}. \code{traits} gives the observed state (0, 1, ...) of
#'     each species; the hidden part is dropped. \code{root_state} marks with a
#'     1 the observed state of the mainland ancestor.
#' }
#'
#' @section Warning:
#' When the same mainland species colonises more than once, and an earlier
#' colonisation has already become endemic, the clade keeps only one of the
#' colonisations. This matches the likelihood, which treats \code{stac} 2, 3
#' and 5 the same way.
#'
#' @seealso \code{\link{TRAISIE_sim_core}}
#'
#' @examples
#' \dontrun{
#' # 2 observed states, no hidden states (n = 2)
#' set.seed(21)
#' trait_pars <- list(immig_rate1 = 0.09, ext_rate1 = 0.95, ana_rate1 = 1.4, clado_rate1 = 0.64,
#'                    immig_rate2 = 0.09, ext_rate2 = 0.35, ana_rate2 = 0.4, clado_rate2 = 0.32,
#'                    trans_rate1 = 0,   trans_rate2 = 1.6,   # from state 1 to states 1, 2
#'                    trans_rate3 = 2.1, trans_rate4 = 0,     # from state 2 to states 1, 2
#'                    K1 = Inf, K2 = Inf,
#'                    p = matrix(0, 2, 2))                    # no anagenesis at changes of state
#' sim <- TRAISIE_sim_CS(total_time = 4,
#'                       mainland = list(M1 = 100, M2 = 150),
#'                       trait_pars = trait_pars,
#'                       replicates = 1,
#'                       num_observed_states = 2,
#'                       num_hidden_states = 1)
#' sim[[1]][[1]]$not_present   # mainland species with no descendants on the island
#' sim[[1]][[2]]               # first clade present
#' }
#' @export
TRAISIE_sim_CS <- function(total_time,
                           mainland,
                           trait_pars,
                           replicates,
                           sample_freq = 100,
                           cond = 0,
                           verbose = TRUE,
                           files_to_write = 0,
                           num_observed_states,
                           num_hidden_states) {
  island_replicates <- list()

  for (rep in 1:replicates) {
    island_replicates[[rep]] <- list()
    full_list <- list()
    if (cond == 0) {
      number_present <- -1
    } else {
      number_present <- 0
    }

    while (number_present < cond) {

      # counts per mainland group (e.g., c(M1, M2, M3, ...))
      counts <- unlist(mainland)
      G      <- length(counts)
      cum    <- cumsum(counts)
      n_tot  <- sum(counts)

      full_list <- vector("list", n_tot)

      for (m_spec in seq_len(n_tot)) {

        # which group does this mainland species belong to?
        g <- which(m_spec <= cum)[1]

        # one-hot root state
        root <- rep(0L, num_observed_states)

        n <- num_observed_states * num_hidden_states
        if (g %in% 1:n) {
          gg <- ceiling(g / num_hidden_states)
        }

        root[gg] <- 1L
        mainland_root <- rep(0L, n)
        mainland_root[g] <- 1L

        # run model with group-specific mainland vector
        full_list[[m_spec]] <- TRAISIE_sim_core(
          time = total_time,
          mainland = as.list(mainland_root),          # e.g., list(1,0,0,...)
          trait_pars = trait_pars,
          num_observed_states = num_observed_states,
          num_hidden_states = num_hidden_states
        )

        if (!is.null(full_list[[m_spec]])) {
          full_list[[m_spec]]$root_state <- root
        }
      }
      stac_vec <- unlist(full_list)[which(names(unlist(full_list)) == "stac")]
      present <- which(stac_vec != 0)
      number_present <- length(present)
    }
    island_replicates[[rep]] <- full_list
    if (verbose == TRUE) {
      message("Island replicate ", rep)
    }
  }

  if (files_to_write == 0) {
    island_replicates <- DAISIE_format_CS(island_replicates =
                                            island_replicates,
                                          time = total_time,
                                          M = sum(unlist(mainland)),
                                          sample_freq = sample_freq,
                                          verbose = verbose)
  }

  if (files_to_write > 0) {
    rm(island_replicates)
    for (filenum in 1:files_to_write) {
      chunks <- ceiling(seq_along(1:replicates) / files_to_write)
      start <- min(which(chunks == filenum))
      end <- max(which(chunks == filenum))
      load(paste0("DAISIE_sims", start, "-", end, ".Rdata"))
      island_replicates <- DAISIE_format_CS(island_replicates =
                                              island_replicates,
                                            time = total_time,
                                            M = sum(unlist(mainland)),
                                            sample_freq = sample_freq,
                                            verbose = verbose)
      save(start, end, island_replicates,
           file = paste0("DAISIE_sims_formatted", start, "-", end, ".Rdata"))
    }
  }

  return(island_replicates)
}
