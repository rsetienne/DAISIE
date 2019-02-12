#' Test if list has area parameters
#'
#' @param area_params object to be tested if conforms to area parameters
#'
#' @return Boolean that indicates if list conforms to expected area parameters
#' as created by \link{create_area_params}
#'
#' @examples
#' testit::assert(are_area_params(
#'   create_area_params(
#'     max_area = 10,
#'     proportional_peak_t = 0.5,
#'     peak_sharpness = 1,
#'     total_island_age = 5)) == TRUE)
are_area_params <- function(area_params) {
  if (is.null(area_params) == TRUE) return(TRUE)
  if (class(area_params) != class(list())) return(FALSE)
  if (!"max_area" %in% names(area_params)) return(FALSE)
  if (!"proportional_peak_t" %in% names(area_params)) return(FALSE)
  if (!"peak_sharpness" %in% names(area_params)) return(FALSE)
  if (!"total_island_age" %in% names(area_params)) return(FALSE)
  if (area_params$max_area < 0.0) return(FALSE)
  if (area_params$proportional_peak_t <= 0.0) return(FALSE)
  if (area_params$proportional_peak_t >= 1.0) return(FALSE)
  if (area_params$peak_sharpness <= 0) return(FALSE)
  if (area_params$total_island_age <= 0.0) return(FALSE)
  TRUE
}
#' Create named list of area parameters
#'
#' @param max_area maximum area
#' @param proportional_peak_t value from 0 to 1 indicating where in the island's history the 
#' peak area is achieved
#' @param peak_sharpness sharpness of peak
#' @param total_island_age total island age
#'
#' @return list of numerical values containing area parameters for island ontogeny 
#' simulation
#' @export
#'
#' @examples 
#' testit::assert(are_area_params(
#' create_area_params(
#'   max_area = 10,
#'   proportional_peak_t = 0.5,
#'   peak_sharpness = 1,
#'   total_island_age = 5
#'       )
#'     )
#'   )
create_area_params <- function(max_area,
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