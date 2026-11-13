#' @keywords internal
TRAISIE_update_rates <- function(timeval,
                                 total_time,
                                 num_spec,
                                 mainland,
                                 trait_pars,
                                 island_spec = NULL,
                                 num_observed_states,
                                 num_hidden_states) {

  # Function to calculate rates at time = timeval. Returns a list with each rate

  # Get immigration rate
  immig_rate <- TRAISIE_get_immig_rate(
    num_spec = num_spec,
    mainland = mainland,
    trait_pars = trait_pars,
    island_spec = island_spec,
    num_observed_states = num_observed_states,
    num_hidden_states = num_hidden_states
  )

  # Get extinction rate
  ext_rate <- TRAISIE_get_ext_rate(
    num_spec = num_spec,
    trait_pars = trait_pars,
    island_spec = island_spec,
    num_observed_states = num_observed_states,
    num_hidden_states = num_hidden_states
  )

  # Get anagenesis rate
  ana_rate <- TRAISIE_get_ana_rate(
    trait_pars = trait_pars,
    island_spec = island_spec,
    num_observed_states = num_observed_states,
    num_hidden_states = num_hidden_states
  )

  # Get cladogenesis rate
  clado_rate <- TRAISIE_get_clado_rate(
    num_spec = num_spec,
    trait_pars = trait_pars,
    island_spec = island_spec,
    num_observed_states = num_observed_states,
    num_hidden_states = num_hidden_states
  )

  # Get transition rate
  trans_rate <- TRAISIE_get_trans_rate(
    trait_pars = trait_pars,
    island_spec = island_spec,
    num_observed_states = num_observed_states,
    num_hidden_states = num_hidden_states
  )

  # Initialize rates list
  rates <- list()

  # Dynamically generate the rates for each state
  for (i in 1:(num_observed_states * num_hidden_states)) {
    rates[[paste("immig_rate", i, sep = "")]] <- immig_rate[[paste("immig_rate", i, sep = "")]]
    rates[[paste("ext_rate", i, sep = "")]]   <- ext_rate[[paste("ext_rate", i, sep = "")]]
    rates[[paste("ana_rate", i, sep = "")]]   <- ana_rate[[paste("ana_rate", i, sep = "")]]
    rates[[paste("clado_rate", i, sep = "")]] <- clado_rate[[paste("clado_rate", i, sep = "")]]
  }

  # Dynamically add transition rates
  for (i in 1:((num_observed_states * num_hidden_states) ^ 2)) {
    rates[[paste("trans_rate", i, sep = "")]]  <-
      trans_rate[[paste("trans_rate", i, sep = "")]]

  }

  return(rates)
}
