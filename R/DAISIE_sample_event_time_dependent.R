#' Samples what event to happen next
#'
#' @inheritParams default_params_doc
#'
#' @return numeric indicating what event will happen, or a supposed event that
#' would happen in some timesteps of the ontogeny algorithm.
#' \itemize{
#'   \item{[1]: immigration event}
#'   \item{[2]: extinction event}
#'   \item{[3]: cladogenesis event}
#'   \item{[4]: anagenesis event}
#' }
#' @author Pedro Neves
DAISIE_sample_event_time_dependent <- function(max_rates) {
  testit::assert(are_max_rates(max_rates))

  possible_event <- sample(1:4,
                           1,
                           replace = FALSE,
                           prob = c(max_rates$immig_max_rate,
                                    max_rates$ext_max_rate,
                                    max_rates$ana_max_rate,
                                    max_rates$clado_max_rate)
  )
  testit::assert(is.numeric(possible_event))
  testit::assert(possible_event >= 1)
  return(possible_event)
}
