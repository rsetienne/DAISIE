#' Determines if object x are rates
#' @param x object to be determined to be rates
#' @export
#'
#' @return TRUE if object x is a list of rates
are_rates <- function(x) {
  if (!all(sapply(x, is.numeric))) return(FALSE)
  if (!"immig_rate" %in% names(x)) return(FALSE)
  if (!"ext_rate" %in% names(x)) return(FALSE)
  if (!"ana_rate" %in% names(x)) return(FALSE)
  if (!"clado_rate" %in% names(x)) return(FALSE)
  if (x$immig_rate < 0.0) return(FALSE)
  if (x$ext_rate < 0.0) return(FALSE)
  if (x$ana_rate < 0.0) return(FALSE)
  if (x$clado_rate < 0.0) return(FALSE)
  TRUE
}

#' Determines if object x are rates
#' @param x object to be determined to be rates
#' @export
#'
#' @return \code{TRUE} if object x is a list of rates, \code{FALSE} otherwise.
are_max_rates <- function(x) {
  if (!all(sapply(x, is.numeric))) return(FALSE)
  if (!"ana_max_rate" %in% names(x)) return(FALSE)
  if (!"ext_max_rate" %in% names(x)) return(FALSE)
  if (!"immig_max_rate" %in% names(x)) return(FALSE)
  if (!"clado_max_rate" %in% names(x)) return(FALSE)
  if (x$ext_max_rate < 0.0) return(FALSE)
  if (x$immig_max_rate < 0.0) return(FALSE)
  if (x$ana_max_rate < 0.0) return(FALSE)
  if (x$clado_max_rate < 0.0) return(FALSE)
  TRUE
}

#' Check if maximum rates are greater or equal to rates
#'
#' @param rates named list of rates as returned by \code{\link{update_rates}}.
#' @param max_rates named list of max rates as returned by
#' \code{\link{update_rates}}.
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
#' @param island_ontogeny Character string that can be \code{"const"},
#' or \code{"beta"} depending on type of island ontogeny desired
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
#' @param sea_level Character string that can be \code{"const"} or
#' \code{"sine"} depending on if sea-level is speified
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


#' Checks if parameters are valid
#'
#' @param pars DAISIE simulation parameters in the form of a named list.
#'
#' @return A boolean stating whether checks are TRUE
#' @export
are_DAISIE_create_sim_pars <- function(pars) {
  if (!"time" %in% names(pars)) return(FALSE)
  if (!"M" %in% names(pars)) return(FALSE)
  if (!"pars" %in% names(pars)) return(FALSE)
  if (!"replicates" %in% names(pars)) return(FALSE)
  if (!"divdepmodel" %in% names(pars)) return(FALSE)
  if (!"ddmodel_sim" %in% names(pars)) return(FALSE)
  if (!"island_type" %in% names(pars)) return(FALSE)
  if (!"nonoceanic_pars" %in% names(pars)) return(FALSE)
  if (!"prop_type2_pool" %in% names(pars)) return(FALSE)
  if (!"replicates_apply_type2" %in% names(pars)) return(FALSE)
  if (!"sample_freq" %in% names(pars)) return(FALSE)
  if (!"plot_sims" %in% names(pars)) return(FALSE)
  if (!"island_ontogeny" %in% names(pars)) return(FALSE)
  if (!"area_pars" %in% names(pars)) return(FALSE)
  if (!"ext_pars" %in% names(pars)) return(FALSE)
  if (!"verbose" %in% names(pars)) return(FALSE)
  if (!pars$time > 0) return(FALSE)
  if (!is.numeric(pars$time)) return(FALSE)
  if (!pars$M > 0) return(FALSE)
  if (!is.numeric(pars$M)) return(FALSE)
  if (!length(pars$pars) == 5 || length(pars$pars) == 10) return(FALSE)
  if (!is.numeric(pars$pars)) return(FALSE)
  if (!pars$replicates >= 1) return(FALSE)
  if (!is.numeric(pars$replicates)) return(FALSE)
  testit::assert(pars$divdepmodel == "CS" || pars$divdepmodel == "IW")
  if (!is.numeric(pars$ddmodel_sim)) return(FALSE)
  testit::assert(pars$island_type == "oceanic" ||
      pars$island_type == "nonoceanic")
  testit::assert(length(pars$nonoceanic_pars) == 2 ||
      is.null(pars$nonoceanic_pars))
  #testit::assert(pars$prop_type2_pool) Pedro write test
  if (!pars$replicates_apply_type2 == TRUE ||
      pars$replicates_apply_type2 == FALSE) return(FALSE)
  if (!is.numeric(pars$sample_freq)) return(FALSE)
  if (!pars$sample_freq > 0) return(FALSE)
  if (!pars$plot_sims == TRUE || pars$plot_sims == FALSE) return(FALSE)
  if (!pars$island_ontogeny == "const" ||
      pars$island_ontogeny == "beta") return(FALSE)
  testit::assert(length(pars$area_pars) == 3 || is.null(pars$area_pars))
  testit::assert(length(pars$ext_pars) == 2 || is.null(pars$ext_pars))
  testit::assert(pars$verbose == TRUE || pars$verbose == FALSE)
  return(TRUE)
}
