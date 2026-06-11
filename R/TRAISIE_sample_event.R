#' Sample a Speciation, Extinction, or Transition Event in Trait-Dependent DAISIE Model
#'
#' This function samples a single evolutionary event (e.g., immigration, extinction,
#' anagenesis, cladogenesis, or trait transition) for the DAISIE model extended to include
#' both observed and hidden trait states. Each state has its own set of rates, and
#' trait transitions are allowed between all state combinations.
#'
#' @param rates A named list of numeric rates. This list must include:
#' \describe{
#'   \item{immig_rateX}{Immigration rate for state X}
#'   \item{ext_rateX}{Extinction rate for state X}
#'   \item{ana_rateX}{Anagenesis rate for state X}
#'   \item{clado_rateX}{Cladogenesis rate for state X}
#'   \item{trans_rateY}{Transition rate for event Y, where Y ranges from 1 to
#'   \code{(num_observed_states * num_hidden_states)^2}}
#' }
#' Here, X ranges from 1 to \code{num_observed_states * num_hidden_states}.
#'
#' @param num_observed_states Integer. Number of observed trait states.
#' @param num_hidden_states Integer. Number of hidden trait states per observed state.
#'
#' @return An integer between 1 and the total number of possible events, indicating the
#' sampled event type and its associated state.
#'
#' @details
#' The total number of possible events is:
#' \deqn{4 \times (\text{num_observed_states} \times \text{num_hidden_states}) +
#'       (\text{num_observed_states} \times \text{num_hidden_states})^2}
#'
#' This includes:
#' This includes:
#' \describe{
#'   \item{State-specific events}{
#'     Immigration, extinction, anagenesis, and cladogenesis events for each state.
#'   }
#'   \item{Transition events}{
#'     Transition events between all combinations of states.
#'   }
#' }
#' The function returns a sampled index (not the event name), so the caller is responsible
#' for mapping the index back to the actual event type if needed.
#'
#' @examples
#' # Example usage with 2 observed and 2 hidden states (i.e., 4 combined states)
#' n_obs <- 2
#' n_hid <- 2
#' num_states <- n_obs * n_hid
#'
#' rates <- list()
#' for (i in 1:num_states) {
#'   rates[[paste0("immig_rate", i)]] <- runif(1, 0.01, 0.1)
#'   rates[[paste0("ext_rate", i)]] <- runif(1, 0.01, 0.1)
#'   rates[[paste0("ana_rate", i)]] <- runif(1, 0.01, 0.1)
#'   rates[[paste0("clado_rate", i)]] <- runif(1, 0.01, 0.1)
#' }
#' for (i in 1:(num_states^2)) {
#'   rates[[paste0("trans_rate", i)]] <- runif(1, 0.01, 0.1)
#' }
#'
#' event <- TRAISIE_sample_event(
#'                   rates,
#'                   num_observed_states = 2,
#'                   num_hidden_states = 2)
#' print(event)
#' @keywords internal
TRAISIE_sample_event <- function(rates,
                                 num_observed_states,
                                 num_hidden_states) {
  # Initialize the probability vector

  # Add the rates for each state
  number_of_possible_events <- 4 * (num_observed_states * num_hidden_states) +
    ((num_observed_states * num_hidden_states) ^ 2)

  prob <- unlist(rates[1:number_of_possible_events])


  # Check for invalid probabilities
  if (any(is.na(prob)) || any(prob < 0)) {
    stop("Invalid probabilities detected in prob vector.")
  }

  # Sample a possible event based on the probabilities
  possible_event <- sample(x = seq_along(prob),
                           size = 1,
                           replace = FALSE,
                           prob = prob)

  return(possible_event)
}
