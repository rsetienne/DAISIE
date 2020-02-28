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
#' @inheritParams default_params_doc
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
#' @inheritParams default_params_doc
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
#' @inheritParams default_params_doc
#'
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
#' @inheritParams default_params_doc
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
#' @inheritParams default_params_doc
#'
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

#' Test if list has trait state parameters
#'
#' @inheritParams default_params_doc
#'
#' @return Boolean that indicates if list conforms to expected area parameters
#' as created by \link{create_trait_pars}
#'
#' @examples
#' testit::assert(DAISIE:::are_trait_pars(
#'   create_trait_pars(
#'     trans_rate = 0.5,
#'     immig_rate2 = 0.1,
#'     ext_rate2 = 0.2,
#'     ana_rate2 = 0.3,
#'     clado_rate2 = 0.4,
#'     trans_rate2 = 0.5,
#'     M2 = 1000)) == TRUE)
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

#' Create named list of trait state parameters
#'
#' @param trans_rate   A numeric with the per capita transition rate with state1
#' @param immig_rate2  A numeric with the per capita immigration rate with state2
#' @param ext_rate2    A numeric with the per capita extinction rate with state2
#' @param ana_rate2    A numeric with the per capita anagenesis rate with state2
#' @param clado_rate2  A numeric with the per capita cladogenesis rate with state2
#' @param trans_rate2  A numeric with the per capita transition rate with state2
#' @param M2           A numeric with the number of species with trait state 2 on mainland
#'
#' @return list of numerical values containing trait state parameters 
#' @export
#'
#' @examples 
#' testit::assert(DAISIE:::are_trait_pars(
#'   create_trait_pars(
#'     trans_rate = 0.5,
#'     immig_rate2 = 0.1,
#'     ext_rate2 = 0.2,
#'     ana_rate2 = 0.3,
#'     clado_rate2 = 0.4,
#'     trans_rate2 = 0.5,
#'     M2 = 1000)) == TRUE)
create_trait_pars <- function(trans_rate,
                              immig_rate2,
                              ext_rate2,
                              ana_rate2,
                              clado_rate2,
                              trans_rate2,
                              M2) {
  assertive::assert_is_a_double(trans_rate)
  assertive::assert_is_a_double(immig_rate2)
  assertive::assert_is_a_double(ext_rate2)
  assertive::assert_is_a_double(ana_rate2)
  assertive::assert_is_a_double(clado_rate2)
  assertive::assert_is_a_double(trans_rate2)
  assertive::assert_is_a_number(M2)
  testit::assert(trans_rate >= 0.0)
  testit::assert(immig_rate2 >= 0.0)
  testit::assert(ext_rate2 >= 0.0)
  testit::assert(ana_rate2 >= 0.0)
  testit::assert(clado_rate2 >= 0.0)
  testit::assert(trans_rate2 >=0.0)
  testit::assert(M2 >=0)
  list(trans_rate = trans_rate,
       immig_rate2 = immig_rate2,
       ext_rate2 = ext_rate2,
       ana_rate2 = ana_rate2,
       clado_rate2 = clado_rate2,
       trans_rate2 = trans_rate2,
       M2 = M2)
}