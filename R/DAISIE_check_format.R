#' Determines if object x are rates
#' @param x object to be determined to be rates
#' @return TRUE if object x is a list of rates
are_rates <- function(x) {
  if (!"immig_rate" %in% names(x)) return(FALSE)
  if (!"ext_rate" %in% names(x)) return(FALSE)
  if (!"ana_rate" %in% names(x)) return(FALSE)
  if (!"clado_rate" %in% names(x)) return(FALSE)
  if (!"ext_rate_max" %in% names(x)) return(FALSE)
  if (!"immig_rate_max" %in% names(x)) return(FALSE)
  if (!"clado_rate_max" %in% names(x)) return(FALSE)
  if (x$immig_rate < 0.0) return(FALSE)
  if (x$ext_rate < 0.0) return(FALSE)
  if (x$ana_rate < 0.0) return(FALSE)
  if (x$clado_rate < 0.0) return(FALSE)
  if (x$ext_rate_max < 0.0 || x$ext_rate_max < x$ext_rate) return(FALSE)
  if (x$immig_rate_max < 0.0 || x$immig_rate_max < x$immig_rate) return(FALSE)
  if (x$clado_rate_max < 0.0 || x$clado_rate_max < x$clado_rate) return(FALSE)
  TRUE
}

#' Check if island_ontogeny is correct after user input
#'
#' @param island_ontogeny Character string that can be \code{"const"},
#' \code{"linear"} or \code{"beta"} depending on type of island ontogeny desired
#' @seealso is_island_ontogeny_runtime
#' @return Boolean stating if island_ontogeny is correct.
#' @export
is_island_ontogeny_input <- function(island_ontogeny) {
  if (class(island_ontogeny) != class(character())) return(FALSE)
  if (island_ontogeny != "const" && island_ontogeny != "linear" && island_ontogeny != "beta") return(FALSE)
  TRUE
}

#' Check if island_ontogeny is correct during runtime (i.e. numeric)
#'
#' @param island_ontogeny Character string that can be \code{"const"},
#' \code{"linear"} or \code{"beta"} depending on type of island ontogeny desired
#' @seealso is_island_ontogeny_runtime
#' @return Boolean stating if island_ontogeny is correct.
#' @export
is_island_ontogeny_runtime <- function(island_ontogeny) {
  if (class(island_ontogeny) != class(numeric())) return(FALSE)
  if (island_ontogeny != 0 && island_ontogeny != 1 && island_ontogeny != 2) return(FALSE)
  TRUE
}

#' Measures if the input is a valid collection of simulation
#' outputs.
#' @param simulation_outputs A list with matrices? of simulation produced by
#' DAISIE_sim.  
#' @return TRUE if the input is a valid collection of simulation
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