#' Determines if object rates are rates
#'
#' @inheritParams default_params_doc
#'
#' @export
#'
#' @return TRUE if object rates is a list of rates
are_rates <- function(rates) {
  if (!all(sapply(rates, is.numeric))) return(FALSE)
  if (!"immig_rate" %in% names(rates)) return(FALSE)
  if (!"ext_rate" %in% names(rates)) return(FALSE)
  if (!"ana_rate" %in% names(rates)) return(FALSE)
  if (!"clado_rate" %in% names(rates)) return(FALSE)
  if (rates$immig_rate < 0.0) return(FALSE)
  if (rates$ext_rate < 0.0) return(FALSE)
  if (rates$ana_rate < 0.0) return(FALSE)
  if (rates$clado_rate < 0.0) return(FALSE)
  if(length(rates) > 4){
    if (!"immig_rate2" %in% names(rates)) return(FALSE)
    if (!"ext_rate2" %in% names(rates)) return(FALSE)
    if (!"ana_rate2" %in% names(rates)) return(FALSE)
    if (!"clado_rate2" %in% names(rates)) return(FALSE)
    if (!"trans_rate" %in% names(rates)) return(FALSE)
    if (!"trans_rate2" %in% names(rates)) return(FALSE)
    if (rates$immig_rate2 < 0.0) return(FALSE)
    if (rates$ext_rate2 < 0.0) return(FALSE)
    if (rates$ana_rate2 < 0.0) return(FALSE)
    if (rates$clado_rate2 < 0.0) return(FALSE)
    if (rates$trans_rate < 0.0) return(FALSE)
    if (rates$trans_rate2 < 0.0) return(FALSE)
    }
  TRUE
}

#' Determines if object max_rates are max_rates
#'
#' @inheritParams default_params_doc
#'
#' @export
#'
#' @return \code{TRUE} if object max_rates is a list of rates,
#' \code{FALSE} otherwise.
are_max_rates <- function(max_rates) {
  if (!all(sapply(max_rates, is.numeric))) return(FALSE)
  if (!"ana_max_rate" %in% names(max_rates)) return(FALSE)
  if (!"ext_max_rate" %in% names(max_rates)) return(FALSE)
  if (!"immig_max_rate" %in% names(max_rates)) return(FALSE)
  if (!"clado_max_rate" %in% names(max_rates)) return(FALSE)
  if (max_rates$ext_max_rate < 0.0) return(FALSE)
  if (max_rates$immig_max_rate < 0.0) return(FALSE)
  if (max_rates$ana_max_rate < 0.0) return(FALSE)
  if (max_rates$clado_max_rate < 0.0) return(FALSE)
  TRUE
}

#' Check if maximum rates are greater or equal to rates
#'
#' @inheritParams default_params_doc
#'
#' @return \code{TRUE} if maximum rates are greater or equal than rates,
#' \code{FALSE} otherwise.
#' @author Joshua Lambert, Pedro Neves
#' @seealso \code{\link{are_rates}}, \code{\link{are_max_rates}}
#'
#' @examples
#' rates <- list(
#'   ext_rate = 0.1,
#'   immig_rate = 0.1,
#'   ana_rate = 0.1,
#'   clado_rate = 0.1
#' )
#' max_rates <- list(
#'   ext_max_rate = 1,
#'   immig_max_rate = 1,
#'   ana_max_rate = 1,
#'   clado_max_rate = 1
#' )
#' testthat::expect_true(
#'   DAISIE:::are_max_rates_gt_rates(
#'     rates = rates,
#'     max_rates = max_rates
#'   )
#' )
are_max_rates_gt_rates <- function(rates, max_rates) {
  if (!all(sapply(rates, is.numeric))) return(FALSE)
  if (!all(sapply(max_rates, is.numeric))) return(FALSE)
  if (max_rates$ext_max_rate < rates$ext_rate)  return(FALSE)
  if (max_rates$immig_max_rate < rates$immig_rate)  return(FALSE)
  if (max_rates$ana_max_rate < rates$ana_rate)  return(FALSE)
  if (max_rates$clado_max_rate < rates$clado_rate)  return(FALSE)
  TRUE
}

#' Check if island_ontogeny is correct after user input
#'
#' @inheritParams default_params_doc
#'
#' @seealso is_island_ontogeny_runtime
#' @return Boolean stating if island_ontogeny is correct.
#' @export
is_island_ontogeny_input <- function(island_ontogeny) {
  if (class(island_ontogeny) != class(character())) return(FALSE)
  if (island_ontogeny != "const" && island_ontogeny != "beta") return(FALSE)
  TRUE
}

#' Check if sea_level is correct after user input
#'
#' @inheritParams default_params_doc
#'
#' @seealso is_sea_level_runtime
#' @return Boolean stating if sea_level is correct.
#' @export
is_sea_level_input <- function(sea_level) {
  if (class(sea_level) != class(character())) return(FALSE)
  if (sea_level != "const" && sea_level != "sine") return(FALSE)
  TRUE
}


#' Measures if the input is a valid collection of simulation
#' outputs.
#'
#' @inheritParams default_params_doc
#'
#' @return TRUE if the input is a valid collection of simulation
#' outputs.
#' @author Richel J.C Bilderbeek, Pedro Neves
#' @examples
#' testthat::expect_false(is_simulation_outputs("nonsense"))
#'
#' simulation_outputs <- DAISIE_sim_constant_rate(
#'   time = 2,
#'   M = 1000,
#'   pars = c(2, 1, 20, 0.0001, 1),
#'   replicates = 1,
#'   plot_sims = FALSE
#'  )
#' testthat::expect_true(is_simulation_outputs(simulation_outputs))
#' @export
is_simulation_outputs <- function(simulation_outputs) {
  for (n_replicate in seq_along(simulation_outputs)) {
    if (!"island_age" %in% names(simulation_outputs[[n_replicate]][[1]]))
      return(FALSE)
    if (!(!"not_present" %in% names(simulation_outputs[[n_replicate]][[1]]) ||
        !"not_present_type1" %in%
        names(simulation_outputs[[n_replicate]][[1]]))) {
      return(FALSE)
    }
    if (!"stt_all" %in% names(simulation_outputs[[n_replicate]][[1]]))
      return(FALSE)
    # TODO: Figure out how to test this?
    # if (!"branching_times" %in% names(simulation_outputs)) return(FALSE)
    # if (!"stac" %in% names(simulation_outputs)) return(FALSE)
    # if (!"missing_species" %in% names(simulation_outputs)) return(FALSE)
  }
  if (is.list(simulation_outputs) && length(simulation_outputs) >= 1) {
    return(TRUE)
  }
}
