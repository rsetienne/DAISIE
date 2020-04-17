#' Calculates algorithm rates
#' @description Internal function that updates the all the rates and
#' max extinction horizon at time t.
#' @family rate calculations
#'
#' @inheritParams default_params_doc
#'
#' @seealso \code{\link{update_max_rates}()}
#'
#' @return a named list with the updated effective rates.
update_rates <- function(timeval,
                         totaltime,
                         gam,
                         laa,
                         lac,
                         mu,
                         hyper_pars = hyper_pars,
                         area_pars = NULL,
                         island_ontogeny = NULL,
                         sea_level = NULL,
                         extcutoff,
                         K,
                         num_spec,
                         num_immigrants,
                         mainland_n) {
  # Function to calculate rates at time = timeval. Returns list with each rate.
  testit::assert(is.numeric(timeval))
  testit::assert(is.numeric(totaltime))
  testit::assert(is.numeric(gam))
  testit::assert(is.numeric(laa))
  testit::assert(is.numeric(lac))
  testit::assert(is.numeric(mu))
  testit::assert(are_hyper_pars(hyper_pars))
  testit::assert(are_area_pars(area_pars))
  testit::assert(is.numeric(island_ontogeny))
  testit::assert(is.numeric(extcutoff) || is.null(extcutoff))
  testit::assert(is.numeric(K))
  testit::assert(is.numeric(num_spec) || is.null(num_spec))
  testit::assert(is.numeric(num_immigrants) || is.null(num_immigrants))
  testit::assert(is.numeric(mainland_n))
  testit::assert(is.numeric(sea_level))

  A <- DAISIE::island_area(
    timeval = timeval,
    area_pars = area_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level
  )

  immig_rate <- get_immig_rate(
    gam = gam,
    A = A,
    num_spec = num_spec,
    K = K,
    mainland_n = mainland_n
  )
  testit::assert(is.numeric(immig_rate))
  ext_rate <- get_ext_rate(
    mu = mu,
    hyper_pars = hyper_pars,
    extcutoff = extcutoff,
    num_spec = num_spec,
    K = K,
    A = A
  )
  testit::assert(is.numeric(ext_rate))
  ana_rate <- get_ana_rate(
    laa = laa,
    num_immigrants = num_immigrants
  )
  testit::assert(is.numeric(ana_rate))
  clado_rate <- get_clado_rate(
    lac = lac,
    hyper_pars = hyper_pars,
    num_spec = num_spec,
    K = K,
    A = A
  )
  testit::assert(is.numeric(clado_rate))


  rates <- list(
    immig_rate = immig_rate,
    ext_rate = ext_rate,
    ana_rate = ana_rate,
    clado_rate = clado_rate
  )
  return(rates)
}

#' Function to describe changes in area through time
#'
#' @inheritParams default_params_doc
#'
#' @export
#' @family rate calculations
#' @author Pedro Neves, Joshua Lambert
#' @references
#' Valente, Luis M., Rampal S. Etienne, and Albert B. Phillimore.
#' "The effects of island ontogeny on species diversity and phylogeny."
#' Proceedings of the Royal Society of London B: Biological
#' Sciences 281.1784 (2014): 20133227.
island_area <- function(timeval,
                        area_pars,
                        island_ontogeny,
                        sea_level) {
  testit::assert(are_area_pars(area_pars))
  Tmax <- area_pars$total_island_age
  Amax <- area_pars$max_area
  Topt <- area_pars$proportional_peak_t
  peak <- area_pars$peak_sharpness
  ampl <- area_pars$sea_level_amplitude
  freq <- area_pars$sea_level_frequency
  theta <- area_pars$island_gradient_angle
  proptime <- timeval / Tmax
  theta <- theta * (pi / 180)
  # Constant ontogeny and sea-level
  if ((island_ontogeny == 0 & sea_level == 0)) {
    if (Amax != 1 || is.null(Amax)) {
      warning("Constant island area requires a maximum area of 1.")
    }
    return(1)
  }

  # Beta function ontogeny and constant sea-level
  if (island_ontogeny == 1 & sea_level == 0) {
    f <- Topt / (1 - Topt)
    a <- f * peak / (1 + f)
    b <- peak / (1 + f)
    At <-
      Amax * proptime ^ a *
      (1 - proptime) ^ b / ((a / (a + b)) ^ a * (b / (a + b)) ^ b)
    return(At)
  }

  if (island_ontogeny == 0 & sea_level == 1) {
    angular_freq <- 2 * pi * freq
    delta_sl <- ampl * sin(proptime * angular_freq)
    r_zero <- sqrt(Amax / pi)
    h_zero <- tan(theta) * r_zero
    h_delta <- max(0, h_zero - delta_sl)
    At <- pi * (h_delta / tan(theta)) ^ 2
    return(At)
  }
  if (island_ontogeny == 1 && sea_level == 1) {
    f <- Topt / (1 - Topt)
    a <- f * peak / (1 + f)
    b <- peak / (1 + f)
    A_beta <-
      Amax * proptime ^ a *
      (1 - proptime) ^ b / ((a / (a + b)) ^ a * (b / (a + b)) ^ b)
    angular_freq <- 2 * pi * freq
    delta_sl <- ampl * sin(proptime * angular_freq)
    r_zero <- sqrt(A_beta / pi)
    h_zero <- tan(theta) * r_zero
    h_delta <- max(0, h_zero - delta_sl)
    At <- pi * (h_delta / tan(theta)) ^ 2
    return(At)
  }
}

#' Function to describe changes in extinction rate through time.
#'
#' @inheritParams default_params_doc
#'
#' @export
#' @family rate calculations
#' @author Pedro Neves, Joshua Lambert
get_ext_rate <- function(mu,
                         hyper_pars,
                         extcutoff = 1000,
                         num_spec,
                         K,
                         A) {

  x <- hyper_pars$x
  ext_rate <- max(0, mu * (A ^ -x) * num_spec, na.rm = TRUE)
  ext_rate <- min(ext_rate, extcutoff, na.rm = TRUE)
  testit::assert(ext_rate >= 0)
  return(ext_rate)
}

#' Calculate anagenesis rate
#' @description Internal function.
#' Calculates the anagenesis rate given the current number of
#' immigrant species and the per capita rate.
#'
#' @inheritParams default_params_doc
#'
#' @export
#' @family rate calculations
#' @author Pedro Neves, Joshua Lambert
get_ana_rate <- function(laa,
                         num_immigrants) {

  ana_rate <- laa * num_immigrants

  testit::assert(is.numeric(ana_rate))
  testit::assert(ana_rate >= 0)
  return(ana_rate)
}

#' Calculate cladogenesis rate
#' @description Internal function.
#' Calculates the cladogenesis rate given the current number of
#' species in the system, the carrying capacity and the per capita cladogenesis
#' rate
#'
#' @inheritParams default_params_doc
#'
#' @export
#' @author Pedro Neves, Joshua Lambert
get_clado_rate <- function(lac,
                           hyper_pars,
                           num_spec,
                           K,
                           A) {
  testit::assert(are_hyper_pars(hyper_pars))

  d <- hyper_pars$d

  clado_rate <- max(
    0, lac * num_spec * (A ^ d) * (1 - num_spec / (K * A)), na.rm = TRUE
  )
  testit::assert(clado_rate >= 0)
  testit::assert(is.numeric(clado_rate))
  return(clado_rate)
}

#' Calculate immigration rate
#' @description Internal function.
#' Calculates the immigration rate given the current number of
#' species in the system, the carrying capacity
#'
#' @inheritParams default_params_doc
#'
#' @export
#' @family rate calculations
#' @author Pedro Neves, Joshua Lambert
get_immig_rate <- function(gam,
                           A,
                           num_spec,
                           K,
                           mainland_n) {

  immig_rate <- max(c(mainland_n * gam * (1 - (num_spec / (A * K))),
                      0), na.rm = TRUE)
  testit::assert(is.numeric(immig_rate))
  testit::assert(immig_rate >= 0)
  return(immig_rate)
}

#' Calculates when the next timestep will be.
#'
#' @param timeval current time of simulation
#' @param max_rates named list of max rates as returned by
#' \code{\link{update_rates}}.
#'
#' @return named list with numeric vector containing the time of the next
#' timestep and the change in time.
#' @author Pedro Neves
calc_next_timeval <- function(max_rates, timeval) {
  testit::assert(timeval >= 0)
  totalrate <- max_rates[[1]] + max_rates[[2]] + max_rates[[3]] + max_rates[[4]]
  dt <- stats::rexp(1, totalrate)
  timeval <- timeval + dt
  return(list(timeval = timeval, dt = dt))
}


#' Calculates when the next timestep will be, and if a shift has occured.
#'
#' @param timeval current time of simulation
#' @param max_rates named list of max rates as returned by
#' \code{\link{update_rates}}.
#' @param dynamic_shift_times numeric vector of times of rate shifts.
#'
#' @return named list with numeric vector containing the time of the next
#' timestep and the change in time.
#' @author Pedro Neves
calc_next_timeval_shift <- function(max_rates,
                                    timeval,
                                    dynamic_shift_times) {
  testit::assert(timeval >= 0)
  totalrate <- max_rates[[1]] + max_rates[[2]] + max_rates[[3]] + max_rates[[4]]
  dt <- stats::rexp(1, totalrate)
  timeval <- timeval + dt
  rate_shift <- FALSE

  if (timeval >= dynamic_shift_times[1]) {
    timeval <- dynamic_shift_times[1]
    dynamic_shift_times <- dynamic_shift_times[-1]
    rate_shift <- TRUE
  }

  out <- list(
    timeval = timeval,
    dt = dt,
    dynamic_shift_times = dynamic_shift_times,
    rate_shift = rate_shift
  )
  return(out)
}

