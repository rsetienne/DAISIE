#' Simulate trait-dependent island replicates
#'
#' @description
#' Simulates island communities under a **trait-dependent** diversification
#' process with optional hidden states. This function runs multiple replicates
#' by calling \code{TRAISIE_sim_core()} once per replicate.
#' If species share a single macro evolutionary process, supply the c
#' orresponding entries in \code{trait_pars} for one observed state and set
#' \code{num_observed_states = 1}. If multiple observed states (and optionally
#' hidden states) are modeled, supply rate/transition parameters per state in
#' \code{trait_pars} and set \code{num_observed_states} and
#' \code{num_hidden_states} accordingly.
#'
#' @usage
#' TRAISIE_sim(
#'   time,
#'   mainland,
#'   trait_pars,
#'   replicates,
#'   num_observed_states,
#'   num_hidden_states
#' )
#'
#' @param time Island age in Myr.
#'   \code{time = 4} simulates the full 4 Myr history;
#' @param mainland List or numeric vector describing the mainland source pool
#' per trait state.
#' @param trait_pars List of trait-dependent parameters consumed by
#'   \code{TRAISIE_sim_core()}, which includes parameters for
#'   cladogenesis, extinction, colonization, anagenesis, and state transitions.
#'   These parameters define the evolutionary processes within the model for
#'   each observed trait state and its possible transitions. The components of
#'   \code{trait_pars} are as follows:
#'
#'   \describe{
#'     \item{immig_rateX}{Immigration rate for trait X (numeric).}
#'     \item{ext_rateX}{Extinction rate for trait X (numeric).}
#'     \item{ana_rateX}{Anagenesis rate for trait X (numeric).}
#'     \item{clado_rateX}{Cladogenesis rate for trait X (numeric).}
#'     \item{\code{trans_rateX}}{A square matrix of transition rates for state changes.
#'   This matrix defines the rates at which species transition between different
#'   observed trait states. Each element in the matrix \( \code{trans_rate}[i,j] \)
#'   represents the rate of transition from state \(i\) to state \(j\).
#'       \describe{
#'         \item{\code{trans_rateX[i, j]}}{The rate at which species in state \(i\) transition to state \(j\).}
#'         \item{Diagonal elements}{\code{trans_rate[i, i]} represent the self-transition rate, and is equal to 0.}
#'       }
#'     }
#'     \item{\code{KX}}{A numeric vector specifying the carrying capacity (\( K \))
#'         for each observed state. This defines the maximum number of species that
#'         can exist in each observed trait state due to ecological or environmental
#'         constraints.}
#'     \item{\code{p}}{A scalar value specifying the probability that a trait
#'         transition between states is accompanied by anagenesis. If \( p = 1 \),
#'         every transition will result in a new species. If \( p = 0 \), the
#'         transition does not lead to the creation of a new species.}
#'   }
#' @param replicates Integer. Number of independent island replicates to simulate.
#' @param num_observed_states Integer (>= 1). Number of **observed** trait states.
#' @param num_hidden_states Integer (>= 1). Number of **hidden** trait states;
#' set to \code{1} if no hidden state is used.
#'
#'
#' @returns
#' A list of length \code{replicates}. Each element is the return value from
#' \code{TRAISIE_sim_core()} for that replicate (or \code{NULL} if
#' the replicate failed).
#'
#' @seealso
#' \code{\link{TRAISIE_sim_core}}
#'
#'
#' @examples
#' \dontrun{
#'set.seed(21)
#'trait_pars = list(immig_rate1 = 0.09,
#'                  ext_rate1 = 0.95,
#'                  ana_rate1 = 1.4,
#'                  clado_rate1 = 0.64,
#'                  immig_rate2 = 0.09,
#'                 ext_rate2 = 0.35,
#'                 ana_rate2 = 0.4,
#'                 clado_rate2 = 0.32,
#'                 trans_rate1 = 0.0,
#'                  trans_rate2 = 1.6,
#'                 trans_rate3 = 2.1,
#'                 trans_rate4 = 0.,
#'                  K1 = Inf,
#'                  K2 = Inf,
#'                  p = 0)
#' data <- TRAISIE_sim (  time = 4,
#'                         mainland = list(M1 = 100, M2 = 150),
#'                         trait_pars = trait_pars,
#'                         replicates = 1,
#'                         num_observed_states = 2,
#'                         num_hidden_states = 1)
#'
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
                                          M = mainland[[1]],
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
                                            M = mainland[[1]],
                                            sample_freq = sample_freq,
                                            verbose = verbose)
      save(start, end, island_replicates,
           file = paste0("DAISIE_sims_formatted", start, "-", end, ".Rdata"))
    }
  }

  return(island_replicates)
}
