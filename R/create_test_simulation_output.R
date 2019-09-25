#' Creates the output of a as-short-as-possible
#' simulation run
#'
#' @param island_ontogeny Boolean indicating if an ontogeny or no ontogeny scenario
#' should be generated. Default is \code{FALSE}.
#'
#' @return simulation_outputs
#' @examples
#'
#' standard_simulation <- create_test_simulation_outputs()
#'
#' island_ontogeny_simulation <- create_test_simulation_outputs(
#'   island_ontogeny = "beta"
#' )
#' @author Pedro Neves, Richel J.C. Bilderbeek
#' @export
create_test_simulation_outputs <- function(island_ontogeny = NULL) {
  set.seed(42)
  if (is.null(island_ontogeny)) {
    pars_equal <- c(2.550687345, 2.683454548, Inf, 0.00933207, 1.010073119)
    DAISIE::DAISIE_sim(
      time = 0.4,
      M = 10,
      pars = pars_equal,
      replicates = 2,
      plot_sims = FALSE,
      verbose = FALSE
    )
  } else {
    testit::assert(is_island_ontogeny_input(island_ontogeny))
    pars_ontogoney_run <- c(7.48223e-05, 1, 0.05, 0.001, 1)
    area_params <- create_area_params(max_area = 10000,
                                      proportional_peak_t = 0.1,
                                      peak_sharpness = 1,
                                      total_island_age = 3)
    DAISIE::DAISIE_sim(
      time = 2,
      M = 500,
      pars = pars_ontogoney_run,
      replicates = 1,
      Apars = area_params,
      Epars = c(0.1, 15),
      island_ontogeny = island_ontogeny,
      plot_sims = FALSE,
      verbose = FALSE
    )
  }
}
