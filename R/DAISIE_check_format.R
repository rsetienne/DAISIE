#' Determines if object rates are rates
#'
#' @inheritParams default_params_doc
#'
#' @return TRUE if object rates is a list of rates
#' @keywords internal
are_rates <- function(rates) {
  # TODO: check if fails on regular 2type model
  if (!all(sapply(rates, is.numeric))) return(FALSE)
  if (!"immig_rate" %in% names(rates)) return(FALSE)
  if (!"ext_rate" %in% names(rates)) return(FALSE)
  if (!"ana_rate" %in% names(rates)) return(FALSE)
  if (!"clado_rate" %in% names(rates)) return(FALSE)
  if (rates$immig_rate < 0.0) return(FALSE)
  if (rates$ext_rate < 0.0) return(FALSE)
  if (rates$ana_rate < 0.0) return(FALSE)
  if (rates$clado_rate < 0.0) return(FALSE)
  if(length(rates) > 4) {
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
#' @return \code{TRUE} if object max_rates is a list of rates,
#' \code{FALSE} otherwise.
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
is_simulation_outputs <- function(simulation_outputs) {
  for (n_replicate in seq_along(simulation_outputs)) {
    if (!"island_age" %in% names(simulation_outputs[[n_replicate]][[1]]))
      return(FALSE)
    if (!(names(simulation_outputs[[n_replicate]][[1]])[2] %in%
        c("not_present","not_present_type1"))) {
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

#' Test if list has area parameters
#'
#' @inheritParams default_params_doc
#'
#' @return Boolean that indicates if list conforms to expected area parameters
#' as created by \link{create_area_pars}
#' @author Richel J.C Bilderbeek, Joshua Lambert, Pedro Neves
#' @noRd
are_area_pars <- function(area_pars) {
  if (is.null(area_pars) == TRUE) return(TRUE)
  if (class(area_pars) != class(list())) return(FALSE)
  if (!"max_area" %in% names(area_pars)) return(FALSE)
  if (!"current_area" %in% names(area_pars)) return(FALSE)
  if (!"proportional_peak_t" %in% names(area_pars)) return(FALSE)
  if (!"total_island_age" %in% names(area_pars)) return(FALSE)
  if (!"sea_level_amplitude" %in% names(area_pars)) return(FALSE)
  if (!"sea_level_frequency" %in% names(area_pars)) return(FALSE)
  if (!"island_gradient_angle" %in% names(area_pars)) return(FALSE)
  if (area_pars$max_area < 0.0) return(FALSE)
  if (area_pars$proportional_peak_t < 0.0) return(FALSE)
  if (area_pars$proportional_peak_t >= 1.0) return(FALSE)
  if (area_pars$total_island_age < 0.0) return(FALSE)
  if (area_pars$sea_level_amplitude < 0.0) return(FALSE)
  if (area_pars$sea_level_frequency < 0.0) return(FALSE)
  if (area_pars$island_gradient_angle < 0.0) return(FALSE)
  if (area_pars$island_gradient_angle > 90) return(FALSE)
  TRUE
}

#' Test if a list has hyperparameters
#'
#' @inheritParams default_params_doc
#'
#' @return \code{TRUE} if list contains hyperparameters, \code{FALSE} otherwise.
#' @author Pedro Neves, Joshua Lambert
#'
#' @noRd
are_hyper_pars <- function(hyper_pars) {
  if (!is.list(hyper_pars)) return(FALSE)
  if (!is.numeric(unlist(hyper_pars))) return(FALSE)
  if (!"d" %in% names(hyper_pars)) return(FALSE)
  if (!"x" %in% names(hyper_pars)) return(FALSE)
  if (!is.numeric(hyper_pars$x)) return(FALSE)
  if (!is.numeric(hyper_pars$d)) return(FALSE)
  TRUE
}

#' Test if list has trait state parameters
#'
#' @inheritParams default_params_doc
#'
#' @return Boolean that indicates if list conforms to expected area parameters
#' as created by \link{create_trait_pars}
#'
#' @noRd
are_trait_pars <- function(trait_pars) {
  if (is.null(trait_pars) == TRUE) return(TRUE)
  if (class(trait_pars) != class(list())) return(FALSE)
  if (!"trans_rate" %in% names(trait_pars)) return(FALSE)
  if (!"immig_rate2" %in% names(trait_pars)) return(FALSE)
  if (!"ext_rate2" %in% names(trait_pars)) return(FALSE)
  if (!"ana_rate2" %in% names(trait_pars)) return(FALSE)
  if (!"clado_rate2" %in% names(trait_pars)) return(FALSE)
  if (!"trans_rate2" %in% names(trait_pars)) return(FALSE)
  if (!"M2" %in% names(trait_pars)) return(FALSE)
  if (trait_pars$trans_rate < 0.0) return(FALSE)
  if (trait_pars$immig_rate2 < 0.0) return(FALSE)
  if (trait_pars$ext_rate2 < 0.0) return(FALSE)
  if (trait_pars$ana_rate2 < 0.0) return(FALSE)
  if (trait_pars$clado_rate2 < 0.0) return(FALSE)
  if (trait_pars$trans_rate2 < 0.0) return(FALSE)
  if (trait_pars$M2 < 0.0) return(FALSE)
  TRUE
}
