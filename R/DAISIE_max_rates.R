#' Calculates the maximum rates for a Gillespie simulation
#' @description Internal function that updates the all the max rates at time t.
#' @family rates calculation
#'
#' @param timeval A numeric with the current time of simulation
#' @param totaltime A numeric with the total time of simulation
#' @param gam A numeric with the per capita immigration rate
#' @param mu A numeric with the per capita extinction rate
#' @param laa A numeric with the per capita anagenesis rate
#' @param lac A numeric with the per capita cladogenesis rate
#' @param ddmodel_sim A numeric determining which parameters are diversity-
#' dependent.
#' @param hyper_pars A vector of hyper parameters
#' @param area_pars stub
#' @param dist_pars stub
#' @param ext_pars stub
#' @param island_ontogeny stub
#' @param sea_level stub
#' @param extcutoff stub
#' @param K stub
#' @param num_spec stub
#' @param num_immigrants stub
#' @param mainland_n stub
#'
#' @return stub
#' @export
#'
#' @examples stub
update_max_rates <- function(timeval,
                             totaltime,
                             gam,
                             mu,
                             laa,
                             lac,
                             ddmodel_sim,
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
    ddmodel_sim = ddmodel_sim,
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
    ddmodel_sim = ddmodel_sim,
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
    mu = mu,
    ddmodel_sim = ddmodel_sim,
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

#' Dynamically update horizon time
#'
#' Update horizon time according to the dynamic maxima of an area/sea-level
#' function.
#'
#' @inheritParams get_t_hor
#' @family rates calculation
#' @return Numeric value with updated time of global maximum area
#' @note At the moment sea-level is set to 0 and only global maximum of function
#' is calculated.
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
#'     timeval = timeval,
#'     totaltime = totaltime,
#'     area_pars = area_pars,
#'     island_ontogeny = island_ontogeny,
#'     sea_level = sea_level,
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
    island_ontogeny = 1,
    sea_level = 0, # Fixed at no sea_level for the moment
    maximum = TRUE,
    tol = .Machine$double.eps
  )
  global_max_area_time <- max$maximum

  testit::assert(is.numeric((global_max_area_time)))
  global_max_area_time <- DDD::roundn(global_max_area_time, 14)
  return(global_max_area_time)
}

#' Dynamically update horizon time
#'
#' Update horizon time according to the dynamic maxima of an area/sea-level
#' function.
#'
#' @inheritParams get_t_hor
#' @family rates calculation
#' @return Numeric value with updated time of global maximum area
#' @note At the moment sea-level is set to 0 and only global maximum of function
#' is calculated.
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
#'   dynamic_t_hor <- DAISIE:::get_global_min_area_time(
#'     timeval = timeval,
#'     totaltime = totaltime,
#'     area_pars = area_pars,
#'     island_ontogeny = island_ontogeny,
#'     sea_level = sea_level,
#'   )
#' )
#'
#' @author Pedro Neves, Joshua Lambert, Shu Xie
get_global_min_area_time <- function(totaltime,
                                     area_pars,
                                     island_ontogeny,
                                     sea_level) {
  # Intervals are temporarily set so the function computes only the global
  # maximum
  interval_min <- 0
  interval_max <- totaltime

  min <- stats::optimize(
    f = DAISIE::island_area,
    interval = c(interval_min, interval_max),
    area_pars = area_pars,
    island_ontogeny = 1,
    sea_level = 0, # Fixed at no sea_level for the moment
    maximum = FALSE,
    tol = .Machine$double.eps
  )
  global_min_area_time <- min$minimum
  testit::assert(is.numeric((global_min_area_time)))
  global_min_area_time <- DDD::roundn(global_min_area_time, 14)
  return(global_min_area_time)
}

