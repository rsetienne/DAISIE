#' Calculates the maximum rates for a Gillespie simulation
#' @description Internal function that updates the all the max rates at time t.
#' @family rate calculations
#'
#' @param timeval A numeric with the current time of simulation
#' @param totaltime A numeric with the total time of simulation
#' @param gam A numeric with the per capita immigration rate
#' @param laa A numeric with the per capita anagenesis rate
#' @param lac A numeric with the per capita cladogenesis rate
#' @param hyper_pars A numeric vector for hyperparameters for the rate
#' calculations, \code{hyper_pars[1]} is d_0 the scaling parameter for
#' exponent for calculating cladogenesis rate, \code{hyper_pars[2]}
#' is x the exponent for calculating extinction rate,
#' \code{hyper_pars[3]} is alpha the exponent for calculating the
#' immigration rate, \code{hyper_pars[4]} is beta the exponent for
#' calculating the anagenesis rate.
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
#' @param dist_pars a numeric for the distance from the mainland.
#' @param ext_pars a numeric vector:
#' \itemize{
#'   \item{[1]: minimum extinction when area is at peak}
#'   \item{[2]: extinction rate when current area is 0.10 of maximum area}
#' }
#' @param island_ontogeny a numeric describing the type of island ontogeny.
#' Can be \code{NULL}, \code{1} for a beta function describing
#' area through time.
#' @param sea_level a numeric describing the type of sea level.
#' @param extcutoff A numeric with the cutoff for extinction rate
#' preventing it from being too large and slowing down simulation.
#' Should be big.
#' @param K a numeric with carrying capacity
#' @param num_spec a numeric with the current number of species.
#' @param num_immigrants a numeric with the current number of non-endemic
#' species (a.k.a non-endemic species).
#' @param mainland_n number of mainland species
#' @seealso \code{\link{update_rates}}
#'
#' @return a named list with the updated effective rates.
#' @export
update_max_rates <- function(timeval,
                             totaltime,
                             gam,
                             laa,
                             lac,
                             hyper_pars = NULL,
                             area_pars,
                             dist_pars = NULL,
                             ext_pars = NULL,
                             island_ontogeny = NULL,
                             sea_level = NULL,
                             extcutoff,
                             K,
                             num_spec,
                             num_immigrants,
                             mainland_n) {

  global_max_area_time <- get_global_max_area_time(
    totaltime = totaltime,
    area_pars = area_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level
  )
  global_min_area_time <- get_global_min_area_time(
    totaltime = totaltime,
    area_pars = area_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level
  )
  testit::assert(is.numeric(global_max_area_time))
  testit::assert(is.finite(global_max_area_time))
  testit::assert(is.numeric(global_min_area_time))
  testit::assert(is.finite(global_min_area_time))

  immig_max_rate <- get_immig_rate(
    timeval = global_max_area_time,
    totaltime = totaltime,
    gam = gam,
    mainland_n = mainland_n,
    hyper_pars = hyper_pars,
    area_pars = area_pars,
    dist_pars = dist_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    num_spec = num_spec,
    K = K
  )

  testit::assert(is.numeric(immig_max_rate))
  clado_max_rate <- get_clado_rate(
    timeval = global_max_area_time,
    lac = lac,
    hyper_pars = hyper_pars,
    area_pars = area_pars,
    dist_pars = dist_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    num_spec = num_spec,
    K = K
  )
  testit::assert(is.numeric(clado_max_rate))

  ext_max_rate <- get_ext_rate(
    timeval = global_min_area_time,
    hyper_pars = hyper_pars,
    area_pars = area_pars,
    ext_pars = ext_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    extcutoff = extcutoff,
    num_spec = num_spec,
    K = K
  )
  testit::assert(is.numeric(ext_max_rate) && ext_max_rate >= 0.0)

  ana_max_rate <- get_ana_rate(
    laa = laa,
    hyper_pars = hyper_pars,
    dist_pars = dist_pars,
    num_immigrants = num_immigrants
  )
  testit::assert(is.numeric(ana_max_rate) && ana_max_rate >= 0.0)

  max_rates <- list(
    ext_max_rate = ext_max_rate,
    immig_max_rate = immig_max_rate,
    ana_max_rate = ana_max_rate,
    clado_max_rate = clado_max_rate
  )
  return(max_rates)
}


#' Get the time of maximum area
#'
#' @param totaltime stub
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
#' @param island_ontogeny a numeric describing the type of island ontogeny.
#' Can be \code{NULL}, \code{1} for a beta function describing
#' area through time.
#' @param sea_level stub
#'
#' @return Numeric with time at which area is maximum during the simulation.
#'
#' @examples
#' timeval <- 1
#' totaltime <- 10
#' area_pars <- DAISIE::create_area_pars(
#'   max_area = 5000,
#'   proportional_peak_t = 0.5,
#'   peak_sharpness = 1,
#'   total_island_age = 15,
#'   sea_level_amplitude = 0,
#'   sea_level_frequency = 0
#' )
#' island_ontogeny <- 1
#' sea_level <- 0
#'
#' testthat::expect_silent(
#'   global_max_area_time <- DAISIE:::get_global_max_area_time(
#'     totaltime = totaltime,
#'     area_pars = area_pars,
#'     island_ontogeny = island_ontogeny,
#'     sea_level = sea_level
#'   )
#' )
#'
#' @author Pedro Neves, Joshua Lambert, Shu Xie
get_global_max_area_time <- function(totaltime,
                                     area_pars,
                                     island_ontogeny,
                                     sea_level) {
  # Intervals are temporarily set so the function computes only the global
  # maximum
  interval_min <- 0
  interval_max <- totaltime

  max <- stats::optimize(
    f = DAISIE::island_area,
    interval = c(interval_min, interval_max),
    area_pars = area_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level, # Fixed at no sea_level for the moment
    maximum = TRUE,
    tol = .Machine$double.eps
  )

  global_max_area_time <- max$maximum

  testit::assert(is.numeric((global_max_area_time)))
  global_max_area_time <- DDD::roundn(global_max_area_time, 14)
  return(global_max_area_time)
}

#' Get the time when area is minimum
#'
#' @param totaltime total time of simulation
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
#' @param island_ontogeny a numeric describing the type of island ontogeny.
#' Can be \code{NULL}, \code{1} for a beta function describing
#' area through time.
#' @param sea_level stub
#'
#' @return Numeric with time at which area is minimum during the simulation
#'
#' @examples
#' timeval <- 1
#' totaltime <- 10
#' area_pars <- DAISIE::create_area_pars(
#'   max_area = 5000,
#'   proportional_peak_t = 0.5,
#'   peak_sharpness = 1,
#'   total_island_age = 15,
#'   sea_level_amplitude = 0,
#'   sea_level_frequency = 0
#' )
#' island_ontogeny <- 1
#' sea_level <- 0
#'
#' testthat::expect_silent(
#'   DAISIE:::get_global_min_area_time(
#'     totaltime = totaltime,
#'     area_pars = area_pars,
#'     island_ontogeny = island_ontogeny,
#'     sea_level = sea_level
#'   )
#' )
#'
#' @author Pedro Neves, Joshua Lambert, Shu Xie
get_global_min_area_time <- function(totaltime,
                                     area_pars,
                                     island_ontogeny,
                                     sea_level) {
  fx <- function(timeval) {
    y <- island_area(
      timeval,
      area_pars = area_pars,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level
    )
    if (is.nan(y)) {
      return(Inf)
    } else {
      return(y)
    }
  }
  global_min_area_time <- subplex::subplex(par = 0, fn = fx)$par
  testit::assert(is.numeric((global_min_area_time)))
  global_min_area_time <- DDD::roundn(global_min_area_time, 14)
  return(global_min_area_time)
}
