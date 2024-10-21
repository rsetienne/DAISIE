#' Calculates algorithm rates
#' @description Internal function that updates the all the rates and
#' max extinction horizon at time t.
#' @family rate calculations
#'
#' @inheritParams default_params_doc
#'
#' @seealso \code{\link{update_max_rates}()}
#' @keywords internal
#' @return a named list with the updated effective rates.
update_rates <- function(timeval,
                         total_time,
                         gam,
                         laa,
                         lac,
                         mu,
                         hyper_pars = hyper_pars,
                         area_pars = NULL,
                         peak = NULL,
                         island_ontogeny = NULL,
                         sea_level = NULL,
                         extcutoff,
                         K,
                         num_spec,
                         num_immigrants,
                         mainland_n,
                         trait_pars = NULL,
                         island_spec = NULL) {
  # Function to calculate rates at time = timeval. Returns list with each rate.
  # testit::assert(is.numeric(timeval))
  # testit::assert(is.numeric(total_time))
  # testit::assert(is.numeric(gam))
  # testit::assert(is.numeric(laa))
  # testit::assert(is.numeric(lac))
  # testit::assert(is.numeric(mu))
  # testit::assert(are_hyper_pars(hyper_pars))
  testit::assert(are_area_pars(area_pars))
  # testit::assert(is.numeric(island_ontogeny))
  # testit::assert(is.numeric(extcutoff) || is.null(extcutoff))
  # testit::assert(is.numeric(K))
  # testit::assert(is.numeric(num_spec) || is.null(num_spec))
  # testit::assert(is.numeric(num_immigrants) || is.null(num_immigrants))
  # testit::assert(is.numeric(mainland_n))
  # testit::assert(is.numeric(sea_level))

  if (!is.null(trait_pars)) {
    return(
      update_rates_trait(
        timeval = timeval,
        total_time = total_time,
        gam = gam,
        mu = mu,
        laa = laa,
        lac = lac,
        hyper_pars = hyper_pars,
        area_pars = area_pars,
        island_ontogeny = island_ontogeny,
        sea_level = sea_level,
        extcutoff = extcutoff,
        K = K,
        mainland_n = mainland_n,
        num_spec = num_spec,
        num_immigrants = num_immigrants,
        trait_pars = trait_pars,
        island_spec = island_spec
      )
    )
  }

  A <- island_area(
    timeval = timeval,
    total_time = total_time,
    area_pars = area_pars,
    peak = peak,
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

  # testit::assert(is.numeric(immig_rate))
  ext_rate <- get_ext_rate(
    mu = mu,
    hyper_pars = hyper_pars,
    extcutoff = extcutoff,
    num_spec = num_spec,
    A = A
  )

  # testit::assert(is.numeric(ext_rate))
  ana_rate <- get_ana_rate(
    laa = laa,
    num_immigrants = num_immigrants
  )
  # testit::assert(is.numeric(ana_rate))
  clado_rate <- get_clado_rate(
    lac = lac,
    hyper_pars = hyper_pars,
    num_spec = num_spec,
    K = K,
    A = A
  )

  # testit::assert(is.numeric(clado_rate))

  rates <- list(
    immig_rate = immig_rate,
    ext_rate = ext_rate,
    ana_rate = ana_rate,
    clado_rate = clado_rate
  )
  return(rates)
}

update_rates_trait <- function(timeval,
                               total_time,
                               gam,
                               laa,
                               lac,
                               mu,
                               hyper_pars = hyper_pars,
                               area_pars = NULL,
                               peak = NULL,
                               island_ontogeny = NULL,
                               sea_level = NULL,
                               extcutoff,
                               K,
                               num_spec,
                               num_immigrants,
                               mainland_n,
                               trait_pars = NULL,
                               island_spec = NULL) {
  # Function to calculate rates at time = timeval. Returns list with each rate.

  A <- island_area(
    timeval = timeval,
    total_time = total_time,
    area_pars = area_pars,
    peak = peak,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level
  )

  immig_rate <- get_immig_rate(
    gam = gam,
    A = A,
    num_spec = num_spec,
    K = K,
    mainland_n = mainland_n,
    trait_pars = trait_pars,
    island_spec = island_spec
  )

  ext_rate <- get_ext_rate(
    mu = mu,
    hyper_pars = hyper_pars,
    extcutoff = extcutoff,
    num_spec = num_spec,
    A = A,
    trait_pars = trait_pars,
    island_spec = island_spec
  )

  ana_rate <- get_ana_rate(
    laa = laa,
    num_immigrants = num_immigrants,
    trait_pars = trait_pars,
    island_spec = island_spec
  )
  clado_rate <- get_clado_rate(
    lac = lac,
    hyper_pars = hyper_pars,
    num_spec = num_spec,
    K = K,
    A = A,
    trait_pars = trait_pars,
    island_spec = island_spec
  )



  trans_rate <- get_trans_rate(trait_pars = trait_pars,
                               island_spec = island_spec)


  rates <- list(
    immig_rate = immig_rate$immig_rate1,
    ext_rate = ext_rate$ext_rate1,
    ana_rate = ana_rate$ana_rate1,
    clado_rate = clado_rate$clado_rate1,
    trans_rate = trans_rate$trans_rate1,
    immig_rate2 = immig_rate$immig_rate2,
    ext_rate2 = ext_rate$ext_rate2,
    ana_rate2 = ana_rate$ana_rate2,
    clado_rate2 = clado_rate$clado_rate2,
    trans_rate2 = trans_rate$trans_rate2,
    M2 = trait_pars$M2)

  return(rates)
}
#' Function to describe changes in area through time.
#'
#' @inheritParams default_params_doc
#'
#' @keywords internal
#' @family rate calculations
#' @author Pedro Neves, Joshua Lambert
#' @references
#' Valente, Luis M., Rampal S. Etienne, and Albert B. Phillimore.
#' "The effects of island ontogeny on species diversity and phylogeny."
#' Proceedings of the Royal Society of London B: Biological
#' Sciences 281.1784 (2014): 20133227.
island_area <- function(timeval,
                        total_time,
                        area_pars,
                        peak,
                        island_ontogeny,
                        sea_level) {
  # testit::assert(are_area_pars(area_pars))
  Tmax <- area_pars$total_island_age
  Amax <- area_pars$max_area
  Acurr <- area_pars$current_area
  proptime_max <- area_pars$proportional_peak_t
  ampl <- area_pars$sea_level_amplitude
  freq <- area_pars$sea_level_frequency
  theta <- area_pars$island_gradient_angle
  proptime <- timeval / Tmax
  proptime_curr <- total_time / Tmax
  theta <- theta * (pi / 180)
  # Constant ontogeny and sea-level
  if (island_ontogeny == 0 & sea_level == 0) {
    if (Amax != 1 || is.null(Amax)) {
      warning("Constant island area requires a maximum area of 1.")
    }
    return(1)
  }

  # Beta function ontogeny and constant sea-level
  if (island_ontogeny == 1 & sea_level == 0) {
    At <- calc_Abeta(proptime = proptime,
                     proptime_max = proptime_max,
                     peak = peak,
                     Amax = Amax)
    return(At)
  }

  if (island_ontogeny == 0 & sea_level == 1) {
    angular_freq <- 2 * pi * freq
    delta_sl <- ampl * cos((proptime_curr - proptime) * angular_freq)
    r_curr <- sqrt((Acurr) / pi)
    h_curr <- tan(theta) * r_curr
    h_delta <- max(0, h_curr - ampl + delta_sl)
    At <- pi * (h_delta / tan(theta)) ^ 2
    return(At)
  }
  if (island_ontogeny == 1 && sea_level == 1) {
    A_beta <- calc_Abeta(proptime,
                         proptime_max,
                         peak,
                         Amax)
    angular_freq <- 2 * pi * freq
    delta_sl <- ampl * cos((proptime_curr - proptime) * angular_freq)
    r_curr <- sqrt(A_beta / pi)
    h_curr <- tan(theta) * r_curr
    h_delta <- max(0, h_curr - ampl + delta_sl)
    At <- pi * (h_delta / tan(theta)) ^ 2
    return(At)
  }
}

#' Function to describe per-capita changes in extinction rate through time
#'
#' This function is only called directly inside the RHS of the ontogeny
#' likelihood functions. In all other cases \code{\link{get_ext_rate}()} is to
#' be called instead.
#'
#' @inheritParams default_params_doc
#'
#' @return Numeric with per capita extinction rate, given A(t), x, and mu0.
#' @noRd
get_ext_rate_per_capita <- function(mu,
                                    x,
                                    extcutoff = 1000,
                                    A = 1) {
  ext_rate_per_capita <- max(0, mu * (A ^ -x), na.rm = TRUE)
  ext_rate_per_capita <- min(ext_rate_per_capita, extcutoff, na.rm = TRUE)
  return(ext_rate_per_capita)
}


#' Calculate extinction rate
#'
#' @inheritParams default_params_doc
#'
#' @return A numeric, with the extinction rate given the base extinction rate,
#' if present, the hyperparemeter \code{x}, A(t) if time-dependent and traits
#' if running a traint model.
#'
#' @keywords internal
#' @family rate calculations
#' @references Valente, Luis M., Rampal S. Etienne, and Albert B. Phillimore.
#' "The effects of island ontogeny on species diversity and phylogeny."
#' Proceedings of the Royal Society of London B: Biological Sciences 281.1784
#' (2014): 20133227.
#' @author Pedro Neves, Joshua Lambert, Shu Xie
get_ext_rate <- function(mu,
                         hyper_pars,
                         extcutoff = 1000,
                         num_spec,
                         A = 1,
                         trait_pars = NULL,
                         island_spec = NULL) {

  x <- hyper_pars$x
  if (is.null(trait_pars)) {

    ext_rate <- num_spec * get_ext_rate_per_capita(
      mu = mu,
      x = x,
      extcutoff = extcutoff,
      A = A
    )
    ext_rate <- min(ext_rate, extcutoff, na.rm = TRUE)
    # testit::assert(ext_rate >= 0)
    return(ext_rate)
  } else {   ##species have two states
    if (is.matrix(island_spec) || is.null(island_spec)) {
      num_spec_trait1 <- length(which(island_spec[, 8] == "1"))
      num_spec_trait2 <- length(which(island_spec[, 8] == "2"))
    }
    ext_rate1 <- mu * num_spec_trait1
    ext_rate2 <- trait_pars$ext_rate2 * num_spec_trait2
    # testit::assert(is.numeric(ext_rate1))
    # testit::assert(is.numeric(ext_rate2))
    # testit::assert(ext_rate1 >= 0)
    # testit::assert(ext_rate2 >= 0)
    ext_list <- list(ext_rate1 = ext_rate1,
                     ext_rate2 = ext_rate2)
    return(ext_list)
  }
}



#' Calculate anagenesis rate
#' @description Internal function.
#' Calculates the anagenesis rate given the current number of
#' immigrant species and the per capita rate.
#'
#' @inheritParams default_params_doc
#'
#' @keywords internal
#' @family rate calculations
#' @author Pedro Neves, Joshua Lambert, Shu Xie
get_ana_rate <- function(laa,
                         num_immigrants,
                         island_spec = NULL,
                         trait_pars = NULL) {

  if(is.null(trait_pars)){
    ana_rate <- laa * num_immigrants

    # testit::assert(is.numeric(ana_rate))
    # testit::assert(ana_rate >= 0)
    return(ana_rate)
  } else {
    ana_rate1 = laa * length(intersect(which(island_spec[,4] == "I"),
                                       which(island_spec[,8] == "1")))
    ana_rate2 = trait_pars$ana_rate2 * length(
      intersect(which(island_spec[,4] == "I"),
                which(island_spec[,8] == "2"))
    )

    # testit::assert(is.numeric(ana_rate1))
    # testit::assert(ana_rate1 >= 0)
    # testit::assert(is.numeric(ana_rate2))
    # testit::assert(ana_rate2 >= 0)
    ana_list <- list(ana_rate1 = ana_rate1,
                     ana_rate2 = ana_rate2)
    return(ana_list)
  }
}

#' Calculate per-capita cladogenesis rate
#'
#' This function is only called directly inside the RHS of the ontogeny
#' likelihood functions. In all other cases \code{\link{get_clado_rate}()} is to
#' be called instead.
#'
#' @inheritParams default_params_doc
#'
#' @return Numeric with the per-capita cladogenesis rate given a base
#' cladogenesis rate, K, A and the d hyperparameter.
#' @noRd
get_clado_rate_per_capita <- function(lac,
                                      d,
                                      num_spec,
                                      K,
                                      A = 1) {
  clado_rate_per_capita <- lac * (A ^ d) * (1 - num_spec / (K * A))
  clado_rate_per_capita <- pmax(0, clado_rate_per_capita, na.rm = TRUE)

  return(clado_rate_per_capita)
}

#' Calculate cladogenesis rate
#' @description Internal function.
#' Calculates the cladogenesis rate given the current number of
#' species in the system, the carrying capacity and the per capita cladogenesis
#' rate
#'
#' @inheritParams default_params_doc
#'
#' @keywords internal
#' @author Pedro Neves, Joshua Lambert, Shu Xie
get_clado_rate <- function(lac,
                           hyper_pars,
                           num_spec,
                           K,
                           A,
                           trait_pars = NULL,
                           island_spec = NULL) {
  # testit::assert(are_hyper_pars(hyper_pars))

  d <- hyper_pars$d
  if (is.null(trait_pars)) {
    clado_rate_pc <- get_clado_rate_per_capita(
      lac = lac,
      d = d,
      num_spec = num_spec,
      K = K,
      A = A
    )
    clado_rate <- num_spec * clado_rate_pc
    # testit::assert(clado_rate >= 0)
    # testit::assert(is.numeric(clado_rate))
    return(clado_rate)
  }else{
    num_spec_trait1 <- length(which(island_spec[, 8] == "1"))
    num_spec_trait2 <- length(which(island_spec[, 8] == "2"))

    if ("K2" %in% names(trait_pars)) {
      clado_rate1 <- max(
        0, lac * num_spec_trait1 * (1 - num_spec_trait1 / K),
        na.rm = TRUE)
      clado_rate2 <- max(
        0, trait_pars$clado_rate2 * num_spec_trait2 * (1 - num_spec_trait2 / trait_pars$K2),
        na.rm = TRUE)
    } else {
      clado_rate1 <- max(
        0, lac * num_spec_trait1 * (1 - num_spec / K),
        na.rm = TRUE)
      clado_rate2 <- max(
        0, trait_pars$clado_rate2 * num_spec_trait2 * (1 - num_spec / K),
        na.rm = TRUE)
    }
    # testit::assert(clado_rate1 >= 0)
    # testit::assert(clado_rate2 >= 0)
    # testit::assert(is.numeric(clado_rate1))
    # testit::assert(is.numeric(clado_rate2))
    clado_list <- list(clado_rate1 = clado_rate1,
                       clado_rate2 = clado_rate2)
    return(clado_list)
  }
}

#' Calculate per-capita immigration rate
#'
#' This function is only called directly inside the RHS of the ontogeny
#' likelihood functions. In all other cases \code{\link{get_immig_rate}()} is to
#' be called instead.
#'
#' @inheritParams default_params_doc
#'
#' @return A numeric with the per-capita immigration rate given A(t) and K.
#' @noRd
get_immig_rate_per_capita <- function(gam,
                                      num_spec,
                                      K,
                                      A = 1) {
  immig_rate_per_capita <- pmax(
    0, gam * (1 - (num_spec / (K * A))), na.rm = TRUE
  )
  return(immig_rate_per_capita)
}

#' Calculate immigration rate
#' @description Internal function.
#' Calculates the immigration rate given the current number of
#' species in the system, the carrying capacity
#'
#' @inheritParams default_params_doc
#'
#' @keywords internal
#' @family rate calculations
#' @author Pedro Neves, Joshua Lambert
#' @references Valente, Luis M., Rampal S. Etienne, and Albert B. Phillimore.
#' "The effects of island ontogeny on species diversity and phylogeny."
#' Proceedings of the Royal Society of London B: Biological Sciences 281.1784 (2014): 20133227.
get_immig_rate <- function(gam,
                           A = 1,
                           num_spec,
                           K,
                           mainland_n,
                           trait_pars = NULL,
                           island_spec = NULL) {

  if (is.null(trait_pars)) {
    immig_rate <- mainland_n * get_immig_rate_per_capita(
      gam = gam,
      num_spec = num_spec,
      K = K,
      A = A
    )
    # testit::assert(is.numeric(immig_rate))
    # testit::assert(immig_rate >= 0)
    return(immig_rate)
  } else {
    mainland_n2 <- trait_pars$M2
    gam2 <- trait_pars$immig_rate2
    num_spec_trait1 <- length(which(island_spec[, 8] == "1"))
    num_spec_trait2 <- length(which(island_spec[, 8] == "2"))
    if ("K2" %in% names(trait_pars)) {
      immig_rate1 <- max(c(mainland_n * gam * (1 - (num_spec_trait1 / K)),
                           0), na.rm = TRUE)
      immig_rate2 <- max(c(mainland_n2 * gam2 * (1 - (num_spec_trait2 / trait_pars$K2)),
                           0), na.rm = TRUE)
    } else {
      immig_rate1 <- max(c(mainland_n * gam * (1 - (num_spec / K)),
                           0), na.rm = TRUE)
      immig_rate2 <- max(c(mainland_n2 * gam2 * (1 - (num_spec / K)),
                           0), na.rm = TRUE)
    }

    # testit::assert(is.numeric(immig_rate1))
    # testit::assert(immig_rate1 >= 0)
    # testit::assert(is.numeric(immig_rate2))
    # testit::assert(immig_rate2 >= 0)
    immig_list <- list(immig_rate1 = immig_rate1,
                       immig_rate2 = immig_rate2)
    return(immig_list)
  }
}

#' Calculate transition rate
#' @description Internal function.
#' Calculates the transition rate given the current number of
#' immigrant species and the per capita rate.
#' @param trait_pars A named list containing diversification rates considering
#' two trait states created by \code{\link{create_trait_pars}}:
#' \itemize{
#'   \item{[1]:A numeric with the per capita transition rate with state1}
#'   \item{[2]:A numeric with the per capita immigration rate with state2}
#'   \item{[3]:A numeric with the per capita extinction rate with state2}
#'   \item{[4]:A numeric with the per capita anagenesis rate with state2}
#'   \item{[5]:A numeric with the per capita cladogenesis rate with state2}
#'   \item{[6]:A numeric with the per capita transition rate with state2}
#'   \item{[7]:A numeric with the number of species with trait state 2 on mainland}
#' }
#' @param island_spec Matrix with current state of simulation containing number
#' of species.
#' @keywords internal
#' @family rates calculation
get_trans_rate <- function(trait_pars,
                           island_spec){

  # Make function accept island_spec matrix or numeric
  if (is.matrix(island_spec) || is.null(island_spec)) {
    num_spec_trait1 <- length(which(island_spec[, 8] == "1"))
    num_spec_trait2 <- length(which(island_spec[, 8] == "2"))
  }
  if ("K2" %in% names(trait_pars)) {
    trans_rate1 <- trait_pars$trans_rate * num_spec_trait1 *
      (1 - (num_spec_trait2 / trait_pars$K2))
    trans_rate2 <- trait_pars$trans_rate2 * num_spec_trait2 *
      (1 - (num_spec_trait1 / trait_pars$K))
  } else {
    trans_rate1 <- trait_pars$trans_rate * num_spec_trait1
    trans_rate2 <- trait_pars$trans_rate2 * num_spec_trait2
  }

  # testit::assert(is.numeric(trans_rate1))
  # testit::assert(trans_rate1 >= 0)
  # testit::assert(is.numeric(trans_rate2))
  # testit::assert(trans_rate2 >= 0)
  trans_list <- list(trans_rate1 = trans_rate1,
                     trans_rate2 = trans_rate2)
  return(trans_list)

}

#' Calculates when the next timestep will be.
#'
#' @param timeval current time of simulation
#' @param max_rates named list of max rates as returned by
#' \code{\link{update_rates}}.
#'
#' @return named list with numeric vector containing the time of the next
#' timestep and the change in time.
#'
#' @keywords internal
#'
#' @author Joshua Lambert, Pedro Neves, Shu Xie
calc_next_timeval <- function(max_rates, timeval, total_time) {
  # testit::assert(timeval >= 0)

  if (length(max_rates) == 4) {   ## no considering about two trait states
    totalrate <-
      max_rates[[1]] + max_rates[[2]] + max_rates[[3]] + max_rates[[4]]
  } else {
    totalrate <-
      max_rates[[1]] + max_rates[[2]] + max_rates[[3]] + max_rates[[4]] +
      max_rates[[5]] + max_rates[[6]] + max_rates[[7]] + max_rates[[8]] +
      max_rates[[9]] + max_rates[[10]]
  }
  if (totalrate != 0) {
      dt <- stats::rexp(1, totalrate)
      timeval <- timeval + dt
  } else {
      timeval <- total_time
  }
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
#' @keywords internal
#'
#' @author Joshua Lambert, Pedro Neves, Shu Xie
calc_next_timeval_shift <- function(max_rates,
                                    timeval,
                                    dynamic_shift_times,
                                    total_time) {
  # testit::assert(timeval >= 0)
  totalrate <- max_rates[[1]] + max_rates[[2]] + max_rates[[3]] + max_rates[[4]]
  if (totalrate != 0) {
    dt <- stats::rexp(1, totalrate)
    timeval <- timeval + dt
  } else {
    timeval <- total_time
  }
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

#' Calculates the area at a point in time from a beta function
#'
#' @inheritParams default_params_doc
#'
#' @keywords internal
#'
#' @author Joshua Lambert, Pedro Neves, Shu Xie
#'
#' @return Numeric
calc_Abeta <- function(proptime,
                       proptime_max,
                       peak,
                       Amax) {
  f <- proptime_max / (1 - proptime_max)
  a <- f * peak / (1 + f)
  b <- peak / (1 + f)
  At <- Amax * proptime ^ a *
    (1 - proptime) ^ b / ((a / (a + b)) ^ a * (b / (a + b)) ^ b)
  return(At)
}

#' Calculates the peak of ontogeny curve (beta function)
#'
#' @inheritParams default_params_doc
#'
#' @keywords internal
#'
#' @return numeric
calc_peak <- function(total_time,
                      area_pars) {
  Amax <- area_pars$max_area
  Acurr <- area_pars$current_area
  proptime_max <- area_pars$proportional_peak_t
  proptime_curr <- total_time / area_pars$total_island_age
  # testit::assert(Acurr <= Amax)
  # testit::assert(proptime_max <= 1 & proptime_max >= 0)
  # testit::assert(proptime_curr <= 1 & proptime_curr >= 0)

  Abeta2 <- function(x) {
    calc_Abeta(proptime_curr, proptime_max, x, Amax) - Acurr
  }
  peak <- stats::uniroot(Abeta2, c(0.01, 1000))$root
  # testit::assert(is.numeric(peak))
  # testit::assert(is.finite(peak))
  return(peak)
}

