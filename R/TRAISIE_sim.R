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
TRAISIE_sim <- function(time,
                                 mainland,
                                 trait_pars,
                                 replicates,
                                 num_observed_states,
                                 num_hidden_states) {

  island_replicates <- vector("list", replicates)
  fail_idx  <- integer(0)
  fail_msg  <- character(0)

  for (rep in seq_len(replicates)) {
    cat(sprintf("replicate %d/%d ... ", rep, replicates))
    utils::flush.console()

    res <- tryCatch(
      {
        TRAISIE_sim_core(
          time = time,
          mainland = mainland,
          trait_pars = trait_pars,
          num_observed_states = num_observed_states,
          num_hidden_states = num_hidden_states
        )
      },
      error = function(e) {
        fail_idx <<- c(fail_idx, rep)
        fail_msg <<- c(fail_msg, conditionMessage(e))
        cat("failed\n")
        NULL
      }
    )

    if (!is.null(res)) cat("done\n")
    island_replicates[[rep]] <- res
  }

  if (length(fail_idx) > 0) {
    attr(island_replicates, "failed") <- list(indices = fail_idx, messages = fail_msg)
    message(sprintf("Completed with %d/%d failures: %s",
                    length(fail_idx), replicates, paste(fail_idx, collapse = ", ")))
  } else {
    message("All replicates completed successfully.")
  }

  island_replicates
}
