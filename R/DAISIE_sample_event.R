#' Samples what event to happen next
#'
#' @param rates numeric list with probability rates for each event. In the 
#' ontogeny case it also contains the maximum possible probability for the 
#' event at each timestep.
#' @param island_ontogeny a string describing the type of island ontogeny. Can be \code{NULL},
#' \code{"quadratic"} for a beta function describing area through time,
#'  or \code{"linear"} for a linear function
#'
#' @return numeric indicating what event will happen, or a supposed event that would
#' happen in some timesteps of the ontogeny algorithm.
#' \itemize{
#'   \item{[1]: immigration event}
#'   \item{[2]: extinction event}
#'   \item{[3]: cladogenesis event}
#'   \item{[4]: anagenesis event}
#'   \item{[5]: proposed extinction that will not happen}
#'   \item{[6]: proposed immigration that will not happen}
#'   \item{[7]: proposed cladogenesis that will not happen}
#' }
#' @author Pedro Neves
DAISIE_sample_event <- function(rates, island_ontogeny = NULL) {
  testit::assert(are_rates(rates))
  
  # If statement prevents odd behaviour of sample when rates are 0
  if (is.null(island_ontogeny)) {
    possible_event <- sample(1:4, 1, prob = c(rates$immig_rate,
                                              rates$ext_rate,
                                              rates$ana_rate,
                                              rates$clado_rate), 
                             replace = FALSE)
  } else {

    possible_event <- sample(1:7, 1, prob = c(
      rates$immig_rate,
      rates$ext_rate,
      rates$ana_rate,
      rates$clado_rate,
      (rates$ext_rate_max - rates$ext_rate),
      (rates$immig_rate_max - rates$immig_rate),
      (rates$clado_rate_max - rates$clado_rate)),
      replace = FALSE)
    
  }
  testit::assert(is.numeric(possible_event))
  testit::assert(possible_event > 0)
  testit::assert(possible_event < 8)
  
  possible_event
}

