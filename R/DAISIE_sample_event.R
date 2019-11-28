#' Samples what event to happen next
#'
#' @param rates numeric list with probability rates for each event. In the
#' ontogeny case it also contains the maximum possible probability for the
#' event at each timestep.
#' @param island_ontogeny a string describing the type of island ontogeny. Can be \code{NULL},
#' \code{"beta"} for a beta function describing area through time.
#' @param sea_level a numeric describing sea_level.
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
DAISIE_sample_event <- function(rates, max_rates) {
  testit::assert(are_rates(rates))
  testit::assert(are_max_rates(max_rates))
  testit::assert(are_max_rates_gt_rates(max_rates = max_rates, rates = rates))
  possible_event <- rng_respecting_sample(1:7, 1, prob = c(
    rates$immig_rate,
    rates$ext_rate,
    rates$ana_rate,
    rates$clado_rate,
    (max_rates$ext_max_rate - rates$ext_rate),
    (max_rates$immig_max_rate - rates$immig_rate),
    (max_rates$clado_max_rate - rates$clado_rate)),
    replace = FALSE)

  testit::assert(is.numeric(possible_event))
  testit::assert(possible_event >= 1)
  return(possible_event)
}
