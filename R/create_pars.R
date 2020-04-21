#' Test if list has area parameters
#'
#' @inheritParams default_params_doc
#'
#' @export
#' @return Boolean that indicates if list conforms to expected area parameters
#' as created by \link{create_area_pars}
#' @author Richel J.C Bilderbeek, Joshua Lambert, Pedro Neves
#'
#' @examples
#' testit::assert(DAISIE::are_area_pars(
#'   create_area_pars(
#'     max_area = 10,
#'     current_area = 4,
#'     proportional_peak_t = 0.5,
#'     total_island_age = 5,
#'     sea_level_amplitude = 5,
#'     sea_level_frequency = 10,
#'     island_gradient_angle = 0)))
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

#' Create named list of area parameters
#'
#' @inheritParams default_params_doc
#'
#' @return list of numerical values containing area and sea level parameters
#' for island ontogeny simulation
#' @export
#' @author Richel J.C Bilderbeek, Joshua Lambert, Pedro Neves
#'
#'
#' @examples
#' testit::assert(DAISIE::are_area_pars(
#' create_area_pars(
#'   max_area = 10,
#'   current_area = 1,
#'   proportional_peak_t = 0.5,
#'   total_island_age = 5,
#'   sea_level_amplitude = 5,
#'   sea_level_frequency = 10,
#'   island_gradient_angle = 0
#'       )
#'     )
#'   )
create_area_pars <- function(max_area,
                             current_area,
                             proportional_peak_t,
                             total_island_age,
                             sea_level_amplitude,
                             sea_level_frequency,
                             island_gradient_angle) {
  testit::assert(max_area > 0.0)
  testit::assert(current_area > 0.0)
  testit::assert(proportional_peak_t >= 0.0)
  testit::assert(proportional_peak_t <= 1.0)
  testit::assert(total_island_age >= 0.0)
  testit::assert(sea_level_amplitude >= 0.0)
  testit::assert(sea_level_frequency >= 0.0)
  testit::assert(island_gradient_angle >= 0)
  testit::assert(island_gradient_angle <= 90)
  list(max_area = max_area,
       current_area = current_area,
       proportional_peak_t = proportional_peak_t,
       total_island_age = total_island_age,
       sea_level_amplitude = sea_level_amplitude,
       sea_level_frequency = sea_level_frequency,
       island_gradient_angle = island_gradient_angle)
}

#' Test if a list has hyperparameters
#'
#' @inheritParams default_params_doc
#'
#' @return \code{TRUE} if list contains hyperparameters, \code{FALSE} otherwise.
#' @export
#' @author Pedro Neves, Joshua Lambert
#'
#' @examples
#' testit::assert(
#'   are_hyper_pars(create_hyper_pars(d = 0.027, x = 0.15))
#' )
are_hyper_pars <- function(hyper_pars) {
  if (!is.list(hyper_pars)) return(FALSE)
  if (!is.numeric(unlist(hyper_pars))) return(FALSE)
  if (!"d" %in% names(hyper_pars)) return(FALSE)
  if (!"x" %in% names(hyper_pars)) return(FALSE)
  if (!is.numeric(hyper_pars$x)) return(FALSE)
  if (!is.numeric(hyper_pars$d)) return(FALSE)
  TRUE
}

#' Create list of hyperparameters
#'
#' @inheritParams default_params_doc
#'
#' @return Named list with hyperparameters
#' @export
#' @author Pedro Neves, Joshua Lambert
#'
#' @examples
#' hyper_pars <- create_hyper_pars(d = 0.027, x = 0.15)
create_hyper_pars <- function(d, x) {
  testit::assert(d >= 0.0)
  testit::assert(is.numeric(x))
  list(
    d = d,
    x = x
  )
}

