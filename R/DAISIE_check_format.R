#' Determines if object x are rates
#' @param x object to be determined to be rates
#' @return TRUE if object x is a list of rates
#' @noRd
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
  if(length(x) == 7){
    if (x$ext_rate_max < 0.0 || x$ext_rate_max < x$ext_rate) return(FALSE)
    if (x$immig_rate_max < 0.0 || x$immig_rate_max < x$immig_rate) return(FALSE)
    if (x$clado_rate_max < 0.0 || x$clado_rate_max < x$clado_rate) return(FALSE)
  }else{
    if (!"immig_rate2" %in% names(x)) return(FALSE)
    if (!"ext_rate2" %in% names(x)) return(FALSE)
    if (!"ana_rate2" %in% names(x)) return(FALSE)
    if (!"clado_rate2" %in% names(x)) return(FALSE)
    if (!"trans_rate" %in% names(x)) return(FALSE)
    if (!"trans_rate2" %in% names(x)) return(FALSE)
    if (x$immig_rate2 < 0.0) return(FALSE)
    if (x$ext_rate2 < 0.0) return(FALSE)
    if (x$ana_rate2 < 0.0) return(FALSE)
    if (x$clado_rate2 < 0.0) return(FALSE)
    if (x$trans_rate < 0.0) return(FALSE)
    if (x$trans_rate2 < 0.0) return(FALSE)
    if (x$ext_rate_max < 0.0 || x$ext_rate_max < x$ext_rate || x$ext_rate_max < x$ext_rate2 ) return(FALSE)
    if (x$immig_rate_max < 0.0 || x$immig_rate_max < x$immig_rate || x$immig_rate_max < x$immig_rate2) return(FALSE)
    if (x$clado_rate_max < 0.0 || x$clado_rate_max < x$clado_rate || x$clado_rate_max < x$clado_rate2) return(FALSE)}
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
#' @param simulation_outputs A list with matrices and vectors of simulation 
#' produced by DAISIE_sim.  
#' @return TRUE if the input is a valid collection of simulation
#' outputs.
#' @author Richel J.C Bilderbeek, Pedro Neves
#' @examples
#' library(testthat)
#'  
#' expect_false(is_simulation_outputs("nonsense")) 
#' 
#' simulation_outputs <- create_test_simulation_outputs() 
#' expect_true(is_simulation_outputs(simulation_outputs))
#' @export
is_simulation_outputs <- function(simulation_outputs) {
  for (n_replicate in seq_along(simulation_outputs)) {
    
    if (!"island_age" %in% names(simulation_outputs[[n_replicate]][[1]])) return(FALSE)
    if (!(!"not_present" %in% names(simulation_outputs[[n_replicate]][[1]]) ||
        !"not_present_type1" %in% names(simulation_outputs[[n_replicate]][[1]]))) {
      return(FALSE)
    }
    if (!"stt_all" %in% names(simulation_outputs[[n_replicate]][[1]])) return(FALSE)
    # TODO: Figure out how to test this?
    # if (!"branching_times" %in% names(simulation_outputs)) return(FALSE)
    # if (!"stac" %in% names(simulation_outputs)) return(FALSE)
    # if (!"missing_species" %in% names(simulation_outputs)) return(FALSE)
  }
  if (is.list(simulation_outputs) && length(simulation_outputs) >= 1) {
    return(TRUE)
  }
  
}
