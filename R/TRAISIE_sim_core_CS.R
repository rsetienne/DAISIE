#' Core simulation for multi-trait, diversity-dependent DAISIE model
#'
#' Simulates island biodiversity dynamics under a clade-specific,
#' multi-trait, diversity-dependent constant-rate process.
#' This function is an extension of the DAISIE core simulation framework,
#' supporting multiple observed and hidden trait states.
#'
#' @param time Numeric scalar. Total simulation time (island age).
#' @param mainland Named list of mainland abundances (e.g., \code{list(M1 = 100, M2 = 150)}).
#' @param extcutoff Numeric. Extinction cutoff threshold. Defaults to \code{1000}.
#' @param trait_pars List of trait-dependent rates and parameters.
#' @param num_observed_states Integer. Number of observed trait states.
#' @param num_hidden_states Integer. Number of hidden trait states.
#'
#' @return A list representing the simulated island, including:
#'   \item{island_age}{Numeric island age.}
#'   \item{stt_table}{State-through-time table.}
#'   \item{island_spec}{Island species data.}
#'   and other DAISIE-like components.
#' @details
#' The simulation proceeds as a continuous-time Markov process,
#' with event rates (immigration, extinction, anagenesis, cladogenesis,
#' trait transitions) recalculated after each event. The function
#' supports multiple observed and hidden states, with final output
#' collapsing hidden states into observed states.
#'
#' @export
TRAISIE_sim_core_CS <- function(
    time,
    mainland,
    extcutoff = 1000,
    trait_pars,
    num_observed_states,
    num_hidden_states
) {

  #### Initialization ####
  timeval <- 0
  total_time <- time

  testit::assert(length(trait_pars) > 5)

  mainland_total <- sum(unlist(mainland))

  testit::assert(mainland_total > 0)
  if (mainland[[1]] != 0) {
    mainland_spec <- seq(1, mainland[[1]], 1)
  } else {
    mainland_spec <- c()
  }
  maxspecID <- mainland_total

  island_spec <- c()
  n <- num_observed_states * num_hidden_states
  stt_table <- matrix(ncol = 3 * (n) + 1)  # 3 for each state (nI, nA, nC) and 1 for Time
  colnames(stt_table) <- c("Time",
                           paste0("nI", 1:n),
                           paste0("nA", 1:n),
                           paste0("nC", 1:n))

  # Initialize the first row
  stt_table[1, ] <- c(total_time, rep(0, 3 * n))


  num_spec <- length(island_spec[, 1])


  #### Start Monte Carlo iterations ####
  while (timeval < total_time) {
    rates <- TRAISIE_update_rates(timeval = timeval,
                                     total_time = total_time,
                                     num_spec = num_spec,
                                     mainland = mainland,
                                     trait_pars = trait_pars,
                                     island_spec = island_spec,
                                     num_observed_states = num_observed_states,
                                     num_hidden_states = num_hidden_states)

    timeval_and_dt <- TRAISIE_calc_next_timeval(
      max_rates = rates,
      timeval = timeval,
      total_time = total_time,
      num_observed_states = num_observed_states,
      num_hidden_states = num_hidden_states)

    timeval <- timeval_and_dt$timeval

    if (timeval < total_time) {
      rates <- TRAISIE_update_rates(timeval = timeval,
                                       total_time = total_time,
                                       num_spec = num_spec,
                                       mainland = mainland,
                                       trait_pars = trait_pars,
                                       island_spec = island_spec,
                                       num_observed_states = num_observed_states,
                                       num_hidden_states = num_hidden_states)


      possible_event <- TRAISIE_sample_event(
        rates = rates,
        num_observed_states = num_observed_states,
        num_hidden_states = num_hidden_states
      )

      #print(possible_event)
      updated_state <- TRAISIE_update_state_CS(
        timeval = timeval,
        total_time = total_time,
        possible_event = possible_event,
        maxspecID = maxspecID,
        island_spec = island_spec,
        stt_table = stt_table,
        trait_pars = trait_pars,
        num_observed_states = num_observed_states,
        num_hidden_states = num_hidden_states
      )

      island_spec <- updated_state$island_spec
      maxspecID   <- updated_state$maxspecID
      stt_table   <- updated_state$stt_table
      num_spec    <- length(island_spec[, 1])

    }
  }

  # Loop through all rows of island_spec
  if (length(island_spec) > 0) {
    for (i in 1:nrow(island_spec)) {

      # Get the current state (convert to numeric)
      state <- as.numeric(island_spec[i, 8])

      # Determine the new block: divide by num_hidden_states, round up, subtract 1
      new_state <- ceiling(state / num_hidden_states) - 1

      # Assign it back as a character
      island_spec[i, 8] <- as.character(new_state)
    }
  }

  #### Finalize STT ####
  stt_table <- rbind(
    stt_table,
    c(
      0,
      stt_table[nrow(stt_table), 2],
      stt_table[nrow(stt_table), 3],
      stt_table[nrow(stt_table), 4]
    )
  )
  island <- TRAISIE_create_island(
    stt_table = stt_table,
    total_time = total_time,
    island_spec = island_spec,
    mainland = mainland,
    trait_pars = trait_pars,
    num_observed_states = num_observed_states,
    num_hidden_states = num_hidden_states)
  ordered_stt_times <- sort(island$stt_table[, 1], decreasing = TRUE)
  testit::assert(all(ordered_stt_times == island$stt_table[, 1]))
  return(island)
}
