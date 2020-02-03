#' Samples what event to happen next
#'
#' @param rates numeric list with probability rates for each event. In the
#' ontogeny case it also contains the maximum possible probability for the
#' event at each timestep.
#' \code{\link{update_rates}}.
#' @return numeric indicating what event will happen, or a supposed event that
#' would happen in some timesteps of the ontogeny algorithm.
#' \itemize{
#'   \item{[1]: immigration event}
#'   \item{[2]: extinction event}
#'   \item{[3]: cladogenesis event}
#'   \item{[4]: anagenesis event}
#' }
#' @author Pedro Neves
DAISIE_sample_event_constant_rate <- function(rates, max_rates) {
  testit::assert(are_rates(rates))
  possible_event <- sample(
    x = 1:4,
    size = 1,
    replace = FALSE,
    prob = c(rates$immig_rate,
             rates$ext_rate,
             rates$ana_rate,
             rates$clado_rate)
  )

  testit::assert(is.numeric(possible_event))
  testit::assert(possible_event >= 1)
  return(possible_event)
}
