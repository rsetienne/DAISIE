#' Measures if the input is a valid collection of simulation
#' outputs.
#' @param simulation_outputs A list with matrices? of simulation produced by
#' DAISIE_sim.  
#' @return TRUE of the input is a valid collection of simulation
#' outputs.
#' @author Richel J.C Bilderbeek, Pedro Neves
#' @examples
#' library(testthat)
#'  
#' expect_false(is_simulation_outputs("nonsense"))   
#' @export
is_simulation_outputs <- function(simulation_outputs) {
  is.list(simulation_outputs)
}