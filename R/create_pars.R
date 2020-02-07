#' Test if list has area parameters
#'
#' @param area_pars object to be tested if conforms to area parameters
#'
#' @export
#' @return Boolean that indicates if list conforms to expected area parameters
#' as created by \link{create_area_pars}
#' @author Richel J.C Bilderbeek, Joshua Lambert, Pedro Neves
#'
#' @examples
#' testit::assert(DAISIE::are_constant_rate_area_pars(area_pars = 1))
are_constant_rate_area_pars <- function(area_pars) {
  if (!is.numeric(area_pars)) return(FALSE)
  if (is.list(area_pars)) return(FALSE)
  if (is.null(area_pars)) return(FALSE)
  if (area_pars < 0.0) return(FALSE)
  TRUE
}

#' Test if list has area parameters
#'
#' @param area_pars object to be tested if conforms to area parameters
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
#'     proportional_peak_t = 0.5,
#'     peak_sharpness = 1,
#'     total_island_age = 5,
#'     sea_level_amplitude = 5,
#'     sea_level_frequency = 10,
#'     island_gradient_angle = 0)))
are_area_pars <- function(area_pars) {
  if (is.null(area_pars) == TRUE) return(TRUE)
  if (class(area_pars) != class(list())) return(FALSE)
  if (!"max_area" %in% names(area_pars)) return(FALSE)
  if (!"proportional_peak_t" %in% names(area_pars)) return(FALSE)
  if (!"peak_sharpness" %in% names(area_pars)) return(FALSE)
  if (!"total_island_age" %in% names(area_pars)) return(FALSE)
  if (!"sea_level_amplitude" %in% names(area_pars)) return(FALSE)
  if (!"sea_level_frequency" %in% names(area_pars)) return(FALSE)
  if (!"island_gradient_angle" %in% names(area_pars)) return(FALSE)
  if (area_pars$max_area < 0.0) return(FALSE)
  if (area_pars$proportional_peak_t < 0.0) return(FALSE)
  if (area_pars$proportional_peak_t >= 1.0) return(FALSE)
  if (area_pars$peak_sharpness < 0) return(FALSE)
  if (area_pars$total_island_age < 0.0) return(FALSE)
  if (area_pars$sea_level_amplitude < 0.0) return(FALSE)
  if (area_pars$sea_level_frequency < 0.0) return(FALSE)
  if (area_pars$island_gradient_angle < 0.0) return(FALSE)
  if (area_pars$island_gradient_angle > 90) return(FALSE)
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
#' area change from sea level.
#' @param island_gradient_angle angle in degrees specifying the slope of the
#' island.
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
#'   proportional_peak_t = 0.5,
#'   peak_sharpness = 1,
#'   total_island_age = 5,
#'   sea_level_amplitude = 5,
#'   sea_level_frequency = 10,
#'   island_gradient_angle = 0
#'       )
#'     )
#'   )
create_area_pars <- function(max_area,
                             proportional_peak_t,
                             peak_sharpness,
                             total_island_age,
                             sea_level_amplitude,
                             sea_level_frequency,
                             island_gradient_angle) {
  testit::assert(max_area > 0.0)
  testit::assert(proportional_peak_t >= 0.0)
  testit::assert(proportional_peak_t <= 1.0)
  testit::assert(peak_sharpness >= 0)
  testit::assert(total_island_age >= 0.0)
  testit::assert(sea_level_amplitude >= 0.0)
  testit::assert(sea_level_frequency >= 0.0)
  testit::assert(island_gradient_angle >= 0)
  testit::assert(island_gradient_angle <= 90)
  list(max_area = max_area,
       proportional_peak_t = proportional_peak_t,
       peak_sharpness = peak_sharpness,
       total_island_age = total_island_age,
       sea_level_amplitude = sea_level_amplitude,
       sea_level_frequency = sea_level_frequency,
       island_gradient_angle = island_gradient_angle)
}

#' Test if a list has hyperparameters
#'
#' @param hyper_pars A named list of hyperparameters as returned
#' by \code{\link{create_hyper_pars}}:
#' \itemize{
#' \item{[1]: is d_0 the scaling parameter for exponent for calculating
#' cladogenesis rate}
#' \item{[2]: is x the exponent for calculating extinction rate}
#' \item{[3]: is alpha, the exponent for calculating the immigration rate}
#' \item{[4]: is beta the exponent for calculating the anagenesis rate.}
#' }
#'
#' @return \code{TRUE} if list contains hyperparameters, \code{FALSE} otherwise.
#' @export
#' @author Pedro Neves, Joshua Lambert
#'
#' @examples
#' testit::assert(
#'   are_hyper_pars(create_hyper_pars(d_0 = 0.027, x = 0.15, 0.294, 0.383))
#' )
are_hyper_pars <- function(hyper_pars) {
  if (!is.list(hyper_pars)) return(FALSE)
  if (!is.numeric(unlist(hyper_pars))) return(FALSE)
  if (!"d_0" %in% names(hyper_pars)) return(FALSE)
  if (!"x" %in% names(hyper_pars)) return(FALSE)
  if (!"alpha" %in% names(hyper_pars)) return(FALSE)
  if (!"beta" %in% names(hyper_pars)) return(FALSE)
  if (hyper_pars$d_0 < 0.0) return(FALSE)
  if (hyper_pars$x < 0.0) return(FALSE)
  if (hyper_pars$alpha < 0.0) return(FALSE)
  if (hyper_pars$beta < 0.0) return(FALSE)
  TRUE
}

#' Create list of hyperparameters
#'
#' @param d_0 the scaling parameter for exponent for calculating
#' cladogenesis rate.
#' @param x the exponent for calculating extinction rate.
#' @param alpha the exponent for calculating the immigration rate.
#' @param beta the exponent for calculating the anagenesis rate.
#'
#' @return Named list with hyperparameters
#' @export
#' @author Pedro Neves, Joshua Lambert
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

#' Create list of distance parameters
#'
#' @param dist_pars A named list of distance parameters
#' as created  by \code{\link{create_dist_pars}}:
#' \itemize{
#' \item{[1]: is D distance from the mainland}
#' }
#' @return \code{TRUE} if list contains distance parameters,
#' \code{FALSE} otherwise.
#' @export
#' @author Pedro Neves, Joshua Lambert
#'
#' @examples
#' testit::assert(are_dist_pars(create_dist_pars(D = 2500)))
are_dist_pars <- function(dist_pars) {
  if (!is.list(dist_pars)) return(FALSE)
  if (!is.numeric(unlist(dist_pars))) return(FALSE)
  if (!"D" %in% names(dist_pars)) return(FALSE)
  if (dist_pars$D < 0) return(FALSE)
  TRUE
}

#' Create list of distance pars
#'
#' @param D A numeric vector for distance parameters for the rate
#' calculations:
#' \itemize{
#' \item{[1]: is D distance from the mainland}
#' }
#'
#' @return Named list with distance parameters
#' @export
#' @author Pedro Neves, Joshua Lambert
#'
#' @examples
#' dist_pars <- create_dist_pars(D = 2500)
create_dist_pars <- function(D) {
  testit::assert(D > 0)
  list(
    D = D
  )
}

#' Creates standard metaparameters to defaults when NULL
#'
#' @param island_ontogeny a numeric describing the type of island ontogeny.
#' @param sea_level a numeric describing the type of sea level.
#' @param area_pars a named list containing area and sea level parameters as
#' created by \code{\link{create_area_pars}}:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#'   \item{[5]: amplitude of area fluctuation from sea level}
#'   \item{[6]: frequency of sine wave of area change from sea level}
#' }
#' @param hyper_pars A numeric vector for hyperparameters for the rate
#' calculations:
#' \itemize{
#' \item{[1]: is d_0 the scaling parameter for exponent for calculating
#' cladogenesis rate}
#' \item{[2]: x, the exponent for calculating extinction rate}
#' \item{[3]: alpha, the exponent for calculating the immigration rate}
#' \item{[4]: beta, the exponent for calculating the anagenesis rate.}
#' }
#' @param dist_pars a numeric for the distance from the mainland.
#' @param ext_pars A numeric vector:
#' \itemize{
#'   \item{[1]: minimum extinction when area is at peak}
#'   \item{[2]: extinction rate when current area is 0.10 of maximum area}
#' }
#' @param totaltime Length of the simulation in time units.
#' @param pars Contains the model parameters:
#' \itemize{
#'   \item{[1]: lambda^c (cladogenesis rate)}
#'   \item{[2]: mu (extinction rate)}
#'   \item{[3]: K (clade-level carrying capacity). Set K=Inf for non-diversity
#'   dependence}
#'   \item{[4]: gamma (immigration rate)}
#'   \item{[5]: lambda^a (anagenesis rate)}
#' }
#' @author Joshua Lambert, Pedro Neves
#' @examples
#' default_pars <- create_default_pars(
#'   island_ontogeny = 0,
#'   sea_level = 0,
#'   area_pars = create_area_pars(
#'     max_area = 1000,
#'     proportional_peak_t = 0.1,
#'     peak_sharpness = 1,
#'     total_island_age = 10,
#'     sea_level_amplitude = 1,
#'     sea_level_frequency = 10,
#'     island_gradient_angle = 0
#'   ),
#'   hyper_pars = create_hyper_pars(
#'     d_0 = 0,
#'     x = 0,
#'     alpha = 0,
#'     beta = 0
#'   ),
#'   dist_pars = create_dist_pars(D = 25000),
#'   ext_pars = c(5, 10),
#'   totaltime = 15,
#'   pars = c(2, 1, 0.05, 0.001, 1)
#' )
#'
#' @export
#' @return Named list with standard parameters and metaparameters needed by
#' \code{DAISIE}.
create_default_pars <- function(island_ontogeny = 0,
                                sea_level = 0,
                                area_pars = NULL,
                                hyper_pars = NULL,
                                dist_pars = NULL,
                                ext_pars = NULL,
                                totaltime,
                                pars) {
  if (island_ontogeny == 0 && sea_level == 0) {
    area_pars <- create_area_pars(
      max_area = 1,
      proportional_peak_t = 0,
      peak_sharpness = 0,
      total_island_age = totaltime,
      sea_level_amplitude = 0,
      sea_level_frequency = 0,
      island_gradient_angle = 0
    )
  }
  if (is.null(hyper_pars)) {
    hyper_pars <- create_hyper_pars(d_0 = 1, x = 1, alpha = 0, beta = 0)
  }
  if (is.null(dist_pars)) {
    dist_pars <- create_dist_pars(D = exp(1))
  }
  if (is.null(ext_pars)) {
    ext_pars <- pars[2]
  }
  testit::assert(is.list(area_pars))
  testit::assert(is.list(hyper_pars))
  testit::assert(is.list(dist_pars))
  testit::assert(is.numeric(ext_pars))
  out <- list(
    area_pars = area_pars,
    hyper_pars = hyper_pars,
    dist_pars = dist_pars,
    ext_pars = ext_pars
  )
  return(out)
}
