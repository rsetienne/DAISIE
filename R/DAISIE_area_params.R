#' Test if list has area parameters
#'
#' @param area_pars object to be tested if conforms to area parameters
#'
#' @export
#' @return Boolean that indicates if list conforms to expected area parameters
#' as created by \link{create_area_pars}
#'
#' @examples
#' testit::assert(DAISIE::are_area_pars(
#'   create_area_pars(
#'     max_area = 10,
#'     proportional_peak_t = 0.5,
#'     peak_sharpness = 1,
#'     total_island_age = 5,
#'     sea_level_amplitude = 5,
#'     sea_level_frequency = 10)) == TRUE)
are_area_pars <- function(area_pars) {
  if (is.null(area_pars) == TRUE) return(TRUE)
  if (class(area_pars) != class(list())) return(FALSE)
  if (!"max_area" %in% names(area_pars)) return(FALSE)
  if (!"proportional_peak_t" %in% names(area_pars)) return(FALSE)
  if (!"peak_sharpness" %in% names(area_pars)) return(FALSE)
  if (!"total_island_age" %in% names(area_pars)) return(FALSE)
  if (!"sea_level_amplitude" %in% names(area_pars)) return(FALSE)
  if (!"sea_level_frequency" %in% names(area_pars)) return(FALSE)
  if (area_pars$max_area < 0.0) return(FALSE)
  if (area_pars$proportional_peak_t <= 0.0) return(FALSE)
  if (area_pars$proportional_peak_t >= 1.0) return(FALSE)
  if (area_pars$peak_sharpness <= 0) return(FALSE)
  if (area_pars$total_island_age <= 0.0) return(FALSE)
  if (area_pars$sea_level_amplitude < 0.0) return(FALSE)
  if (area_pars$sea_level_frequency < 0.0) return(FALSE)
  TRUE
}

#' Create named list of area parameters
#'
#' @param max_area maximum area
#' @param proportional_peak_t value from 0 to 1 indicating
#' where in the island's history the peak area is achieved
#' @param peak_sharpness sharpness of peak
#' @param total_island_age total island age
#' @param sea_level_amplitude amplitude of area fluctuation
#' from sea level
#' @param sea_level_frequency frequency of sine wave of
#' area change from sea level
#'
#' @return list of numerical values containing area parameters
#' for island ontogeny simulation
#' @export
#'
#' @examples
#' testit::assert(DAISIE::are_area_pars(
#' create_area_pars(
#'   max_area = 10,
#'   proportional_peak_t = 0.5,
#'   peak_sharpness = 1,
#'   total_island_age = 5,
#'   sea_level_amplitude = 5,
#'   sea_level_frequency = 10
#'       )
#'     )
#'   )
create_area_pars <- function(max_area,
                             proportional_peak_t,
                             peak_sharpness,
                             total_island_age,
                             sea_level_amplitude,
                             sea_level_frequency){
  testit::assert(max_area > 0.0)
  testit::assert(proportional_peak_t >= 0.0)
  testit::assert(proportional_peak_t <= 1.0)
  testit::assert(peak_sharpness >= 0)
  testit::assert(total_island_age >= 0.0)
  testit::assert(sea_level_amplitude >= 0.0)
  testit::assert(sea_level_frequency >= 0.0)
  list(max_area = max_area,
       proportional_peak_t = proportional_peak_t,
       peak_sharpness = peak_sharpness,
       total_island_age = total_island_age,
       sea_level_amplitude = sea_level_amplitude,
       sea_level_frequency = sea_level_frequency)
}
