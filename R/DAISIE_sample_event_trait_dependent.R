#' Samples what event to happen next
#'
#' @inheritParams default_params_doc
#'
#' @return numeric indicating what event will happen, or a supposed event that
#' would happen in some timesteps of the ontogeny algorithm.
#' \itemize{
#'   \item{[1]: immigration event with trait1}
#'   \item{[2]: extinction event with trait1}
#'   \item{[3]: cladogenesis event with trait1}
#'   \item{[4]: anagenesis event with trait1}
#'   \item{[5]: transition event with trait1}
#'   \item{[6]: immigration event with trait2}
#'   \item{[7]: extinction event with trait2}
#'   \item{[8]: cladogenesis event with trait2}
#'   \item{[9]: anagenesis event with trait2}
#'   \item{[10]: transition event with trait2}
#' }
#' @author Pedro Neves
DAISIE_sample_event_trait_dependent <- function(rates) {
  testit::assert(are_rates(rates))
  
  possible_event <- sample(1:10,
                           1,
                           replace = FALSE,
                           prob = c(rates$immig_rate,
                                    rates$ext_rate,
                                    rates$ana_rate,
                                    rates$clado_rate,
                                    rates$trans_rate,
                                    rates$immig_rate2,
                                    rates$ext_rate2,
                                    rates$ana_rate2,
                                    rates$clado_rate2,
                                    rates$trans_rate2)
  )
  testit::assert(is.numeric(possible_event))
  testit::assert(possible_event >= 1)
  return(possible_event)
}

# rng_respecting_sample <- function (x, size, replace, prob)
# {
#   which_non_zero <- prob > 0
#   non_zero_prob <- prob[which_non_zero]
#   non_zero_x <- x[which_non_zero]
#   return(DDD::sample2(x = non_zero_x, size = size, replace = replace,
#                       prob = non_zero_prob))
# }
