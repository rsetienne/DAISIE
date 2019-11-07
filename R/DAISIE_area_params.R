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
#'     total_island_age = 5)) == TRUE)
are_area_pars <- function(area_pars) {
  if (is.null(area_pars) == TRUE) return(TRUE)
  if (class(area_pars) != class(list())) return(FALSE)
  if (!"max_area" %in% names(area_pars)) return(FALSE)
  if (!"proportional_peak_t" %in% names(area_pars)) return(FALSE)
  if (!"peak_sharpness" %in% names(area_pars)) return(FALSE)
  if (!"total_island_age" %in% names(area_pars)) return(FALSE)
  if (area_pars$max_area < 0.0) return(FALSE)
  if (area_pars$proportional_peak_t <= 0.0) return(FALSE)
  if (area_pars$proportional_peak_t >= 1.0) return(FALSE)
  if (area_pars$peak_sharpness <= 0) return(FALSE)
  if (area_pars$total_island_age <= 0.0) return(FALSE)
  TRUE
}

#' Create named list of area parameters
#'
#' @param max_area maximum area
#' @param proportional_peak_t value from 0 to 1 indicating
#' where in the island's history the peak area is achieved
#' @param peak_sharpness sharpness of peak
#' @param total_island_age total island age
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
#'   total_island_age = 5
#'       )
#'     )
#'   )
create_area_pars <- function(max_area,
                               proportional_peak_t,
                               peak_sharpness,
                               total_island_age) {
  testit::assert(max_area > 0.0)
  testit::assert(proportional_peak_t >= 0.0)
  testit::assert(proportional_peak_t <= 1.0)
  testit::assert(peak_sharpness >= 0)
  testit::assert(total_island_age >= 0.0)
  list(max_area = max_area,
       proportional_peak_t = proportional_peak_t,
       peak_sharpness = peak_sharpness,
       total_island_age = total_island_age)
}

#' Create named list of sea-level parameters
#'
#' @param amplitude the amplitude of the sine function
#' @param frequency the frequency of the sine function
#' @param phase the phase of the sine function
#'
#' @return list of numerical values containing parameters
#' for sea-level simulation
#' @export
#'
#' @examples
#' create_sea_level_pars(
#'   amplitude = 50,
#'   frequency = 2,
#'   phase = 0)
create_sea_level_pars <- function(amplitude,
                                  frequency,
                                  phase) {
  testit::assert(amplitude > 0.0)
  testit::assert(frequency > 0.0)
  testit::assert(phase >= 0)
  list(amplitude = amplitude,
       frequency = frequency,
       phase = phase)
}
