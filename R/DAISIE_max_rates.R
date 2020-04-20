#' Calculates the maximum rates for a Gillespie simulation
#' @description Internal function that updates the all the max rates at time t.
#' @family rate calculations
#'
#' @inheritParams default_params_doc
#'
#' @seealso \code{\link{update_rates}()}
#'
#' @return a named list with the updated effective rates.
#' @export
update_max_rates <- function(gam,
                             laa,
                             lac,
                             mu,
                             hyper_pars = NULL,
                             extcutoff,
                             K,
                             num_spec,
                             num_immigrants,
                             mainland_n,
                             Amin,
                             Amax) {

  immig_max_rate <- get_immig_rate(
    gam = gam,
    A = Amax,
    num_spec = num_spec,
    K = K,
    mainland_n = mainland_n
  )

  testit::assert(is.numeric(immig_max_rate))
  clado_max_rate <- get_clado_rate(
    lac = lac,
    hyper_pars = hyper_pars,
    num_spec = num_spec,
    K = K,
    A = Amax
  )
  testit::assert(is.numeric(clado_max_rate))

  ext_max_rate <- get_ext_rate(
    mu = mu,
    hyper_pars = hyper_pars,
    extcutoff = extcutoff,
    num_spec = num_spec,
    K = K,
    A = Amin
  )
  testit::assert(is.numeric(ext_max_rate) && ext_max_rate >= 0.0)

  ana_max_rate <- get_ana_rate(
    laa = laa,
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


#' Get the maximum area
#'
#' @inheritParams default_params_doc
#'
#' @return Numeric maximum area during the simulation.
#'
#' @examples
#' timeval <- 1
#' totaltime <- 10
#' area_pars <- DAISIE::create_area_pars(
#'   max_area = 5000,
#'   current_area = 4000,
#'   proportional_peak_t = 0.5,
#'   total_island_age = 15,
#'   sea_level_amplitude = 0,
#'   sea_level_frequency = 0,
#'   island_gradient_angle = 0
#' )
#' island_ontogeny <- 1
#' sea_level <- 0
#' peak <- DAISIE:::calc_peak(totaltime = totaltime, area_pars = area_pars)
#' testthat::expect_silent(
#'   global_max_area_time <- DAISIE:::get_global_max_area(
#'     totaltime = totaltime,
#'     peak = peak,
#'     area_pars = area_pars,
#'     island_ontogeny = island_ontogeny,
#'     sea_level = sea_level
#'   )
#' )
#'
#' @author Pedro Neves, Joshua Lambert, Shu Xie
get_global_max_area <- function(totaltime,
                                area_pars,
                                peak,
                                island_ontogeny,
                                sea_level) {

  max <- stats::optimize(
    f = DAISIE::island_area,
    interval = c(0, totaltime),
    area_pars = area_pars,
    peak = peak,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level, # Fixed at no sea_level for the moment
    maximum = TRUE,
    tol = .Machine$double.eps
  )

  global_max_area_time <- max$maximum

  testit::assert(is.numeric((global_max_area_time)))
  global_max_area_time <- DDD::roundn(global_max_area_time, 14)

  Amax <- DAISIE::island_area(
    timeval = global_max_area_time,
    area_pars = area_pars,
    peak = peak,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level
  )
  return(Amax)
}

#' Get the minimum area
#'
#' @inheritParams default_params_doc
#'
#' @return Numeric with time at which area is minimum during the simulation
#'
#' @examples
#' timeval <- 1
#' totaltime <- 10
#' area_pars <- DAISIE::create_area_pars(
#'   max_area = 5000,
#'   current_area = 4000,
#'   proportional_peak_t = 0.5,
#'   total_island_age = 15,
#'   sea_level_amplitude = 0,
#'   sea_level_frequency = 0,
#'   island_gradient_angle = 0
#' )
#' island_ontogeny <- 1
#' sea_level <- 0
#' peak <- DAISIE:::calc_peak(totaltime = totaltime, area_pars = area_pars)
#'
#' testthat::expect_silent(
#'   DAISIE:::get_global_min_area(
#'     totaltime = totaltime,
#'     area_pars = area_pars,
#'     peak = peak,
#'     island_ontogeny = island_ontogeny,
#'     sea_level = sea_level
#'   )
#' )
#'
#' @author Pedro Neves, Joshua Lambert, Shu Xie
get_global_min_area <- function(totaltime,
                                area_pars,
                                peak,
                                island_ontogeny,
                                sea_level) {
  fx <- function(timeval) {
    y <- island_area(
      timeval,
      area_pars = area_pars,
      peak = peak,
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

  Amin <- DAISIE::island_area(
    timeval = global_min_area_time,
    area_pars = area_pars,
    peak = peak,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level
  )
  return(Amin)
}
