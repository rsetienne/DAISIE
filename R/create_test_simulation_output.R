#' Creates the output of a as-short-as-possible
#' simulation run
#' @return simulation_outputs
#' @author Pedro Neves, Richel J.C. Bilderbeek
#' @export
create_test_simulation_outputs <- function() {
  pars_equal = c(2.550687345,2.683454548,Inf,0.00933207,1.010073119)
  DAISIE::DAISIE_sim(
    time = 0.4,
    M = 1000,
    pars = pars_equal,
    replicates = 2
  )
}
