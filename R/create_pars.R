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
  if (area_pars$proportional_peak_t < 0.0) return(FALSE)
  if (area_pars$proportional_peak_t >= 1.0) return(FALSE)
  if (area_pars$peak_sharpness < 0) return(FALSE)
  if (area_pars$total_island_age < 0.0) return(FALSE)
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
#' @return list of numerical values containing area and sea level parameters
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
                             sea_level_frequency) {
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

are_hyper_pars <- function(variables) {
  #stubby stub stub
}


#' Create list of hyperparameters
#'
#' @param hyper_pars A numeric vector for hyperparameters for the rate
#' calculations:
#' \itemize{
#' \item{[1]: is d_0 the scaling parameter for exponent for calculating
#' cladogenesis rate}
#' \item{[2]: is x the exponent for calculating extinction rate}
#' \item{[3]: is alpha, the exponent for calculating the immigration rate}
#' \item{[4]: is beta the exponent for calculating the anagenesis rate.}
#' }
#'
#' @return Named list with hyperparameters
#' @export
#'
#' @examples
#' hyper_pars <- create_hyper_pars(d_0 = 0.027, x = 0.15, 0.294, 0.383)
create_hyper_pars <- function(d_0, x, alpha, beta) {
  testit::assert(d_0 >= 0.0)
  testit::assert(x >= 0.0)
  testit::assert(alpha >= 0.0)
  testit::assert(beta >= 0.0)
  list(
    d_0 = d_0,
    x = x,
    alpha = alpha,
    beta = beta
  )
}

#' Create list of distance pars
#'
#' @param dist_pars A numeric vector for distance parameters for the rate
#' calculations:
#' \itemize{
#' \item{[1]: is D distance from the mainland}
#' }
#'
#' @return Named list with distance parameters
#' @export
#'
#' @examples
#' dist_pars <- create_dist_pars(D = 2500)
create_dist_pars <- function(D) {
  testit::assert(D > 0)
  list(
    D = D
  )
}
