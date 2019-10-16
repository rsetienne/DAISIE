#' Samples what event to happen next
#'
#' @param rates numeric list with probability rates for each event. In the
#' ontogeny case it also contains the maximum possible probability for the
#' event at each timestep.
#' @param island_ontogeny a string describing the type of island ontogeny. Can be \code{NULL},
#' \code{"beta"} for a beta function describing area through time,
#'  or \code{"linear"} for a linear function
#' @param Tpars A named list containing diversification rates considering two trait states:
#' \itemize{
#'   \item{[1]:A numeric with the per capita transition rate with state1}
#'   \item{[2]:A numeric with the per capita immigration rate with state2}
#'   \item{[3]:A numeric with the per capita extinction rate with state2}
#'   \item{[4]:A numeric with the per capita anagenesis rate with state2}
#'   \item{[5]:A numeric with the per capita cladogenesis rate with state2}
#'   \item{[6]:A numeric with the per capita transition rate with state2}
#'   \item{[7]:A numeric with the number of species with trait state 2 on mainland}
#' }
#' @return numeric indicating what event will happen, or a supposed event that would
#' happen in some timesteps of the ontogeny algorithm.
#' \itemize{
#'   \item{[1]: immigration event}
#'   \item{[2]: extinction event}
#'   \item{[3]: cladogenesis event}
#'   \item{[4]: anagenesis event}
#'   \item{[5]: proposed extinction that will not happen or transition event}
#'   \item{[6]: proposed immigration that will not happen or immigration event with state2}
#'   \item{[7]: proposed cladogenesis that will not happen or extinction event with state2}
#'   \item{[8]: cladogenesis event with state2}
#'   \item{[9]: anagenesis event with state2}
#'   \item{[10]: transition event with state2}
#' }
#' @author Pedro Neves
DAISIE_sample_event <- function(rates, island_ontogeny = NULL, Tpars = NULL) {
  testit::assert(are_rates(rates))

  #testit::assert(DAISIE::is_island_ontogeny_runtime(island_ontogeny))

  if(is.null(Tpars)){
    # If statement prevents odd behaviour of sample when rates are 0
    if (island_ontogeny == 0) {
      possible_results <- 1:4 # Each number is a different type of event
      n_events <- 1 # We only need 1 event
      testit::assert(!is.null(rates$immig_rate))
      testit::assert(!is.null(rates$ext_rate))
      testit::assert(!is.null(rates$ana_rate))
      testit::assert(!is.null(rates$clado_rate))
      event_probabilities <- c(
        rates$immig_rate,
        rates$ext_rate,
        rates$ana_rate,
        rates$clado_rate
      )
      testit::assert(all(event_probabilities >= 0.0))
      testit::assert(length(possible_results) == length(event_probabilities))
      possible_event <- DDD::rng_respecting_sample(
        x = possible_results,
        size = n_events,
        replace = TRUE, # irrelevant for 1 draw
        prob = event_probabilities
      )
    } else {
      event_probabilities <- c(
        rates$immig_rate,
        rates$ext_rate,
        rates$ana_rate,
        rates$clado_rate,
        (rates$ext_rate_max - rates$ext_rate),
        (rates$immig_rate_max - rates$immig_rate),
        (rates$clado_rate_max - rates$clado_rate)
      )
      testit::assert(length(event_probabilities) == 7)
      possible_event <- DDD::rng_respecting_sample(
        1:7, 1,
        prob = event_probabilities,
        replace = FALSE
      )
    }
    testit::assert(is.numeric(possible_event))
    testit::assert(possible_event >= 1)
    testit::assert(possible_event <= (island_ontogeny == 0) * 4 + (island_ontogeny > 0) * 7)
  }else{
    possible_results <- 1:10 # Each number is a different type of event
    n_events <- 1 # We only need 1 event
    testit::assert(!is.null(rates$immig_rate))
    testit::assert(!is.null(rates$ext_rate))
    testit::assert(!is.null(rates$ana_rate))
    testit::assert(!is.null(rates$clado_rate))
    # event_probabilities <- c(
    #   rates$immig_rate,
    #   rates$ext_rate,
    #   rates$ana_rate,
    #   rates$clado_rate,
    #   rates$trans_rate,
    #   rates$immig_rate2,
    #   rates$ext_rate2,
    #   rates$ana_rate2,
    #   rates$clado_rate2,
    #   rates$trans_rate2
    # )
    event_probabilities <- c(
      rates$immig_rate,
      rates$immig_rate2,
      rates$ext_rate,
      rates$ext_rate2,
      rates$ana_rate,
      rates$ana_rate2,
      rates$clado_rate,
      rates$clado_rate2,
      rates$trans_rate,
      rates$trans_rate2
    )
    testit::assert(length(event_probabilities) == 10)
    # possible_event <- DDD::rng_respecting_sample(1:10, 1, prob = event_probabilities,
    #                                              replace = FALSE)
    possible_event <- DDD::rng_respecting_sample(c(1,6,2,7,3,8,4,9,5,10), 1, prob = event_probabilities,
                                                 replace = FALSE)
    testit::assert(is.numeric(possible_event))
    testit::assert(possible_event >= 1)
    testit::assert(possible_event <= 10)
  }
  return(possible_event)
}

