#' Calculates algorithm rates
#' @description Internal function that updates the all the rates and
#' max extinction horizon at time t.
#' @family rates calculation
#'
#' @param timeval A numeric with the current time of simulation
#' @param totaltime A numeric with the total time of simulation
#' @param gam A numeric with the per capita immigration rate
#' @param mu A numeric with the per capita extinction rate in no ontogeny model
#' @param laa A numeric with the per capita anagenesis rate
#' @param lac A numeric with the per capita cladogenesis rate
#' @param ddmodel_sim A numeric determining which parameters are diversity-
#' dependent.
#' @param area_pars A named list containing area parameters as created
#' by create_area_pars:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#' }
#' @param ext_pars A numeric vector:
#' \itemize{
#'   \item{[1]: minimum extinction when area is at peak}
#'   \item{[2]: extinction rate when current area is 0.10 of maximum area}
#' }
#' @param island_ontogeny A string describing the type of island ontogeny.
#' Can be \code{NULL},
#' \code{"beta"} for a beta function describing area through time.
#' @param extcutoff A numeric with the cutoff for extinction rate
#' preventing it from being too
#' large and slowing down simulation. Should be big.
#' @param K A numeric with K (clade-specific carrying capacity)
#' @param island_spec A matrix containing state of system
#' @param mainland_n A numeirc with the total number of species present
#' in the mainland
#' @param t_hor A numeric with the time of horizon for max cladogenesis,
#' immigration and minimum extinction
#' @param hyper_pars A numeric vector for hyperparameters for the rate
#' calculations, \code{hyper_pars[1]} is d_0 the scaling parameter for
#' exponent for calculating cladogenesis rate, \code{hyper_pars[2]}
#' is x the exponent for calculating extinction rate,
#' \code{hyper_pars[3]} is alpha the exponent for calculating the
#' immigration rate, \code{hyper_pars[4]} is beta the exponent for
#' calculating the anagenesis rate.
#' @param sea_level a numeric describing the type of sea level.
#' @param dist_pars a numeric for the distance from the mainland.
update_rates <- function(timeval, totaltime,
                         gam, mu, laa, lac, ddmodel_sim = ddmodel_sim,
                         hyper_pars = hyper_pars,
                         area_pars = NULL,
                         dist_pars = NULL,
                         ext_pars = NULL,
                         island_ontogeny = NULL,
                         sea_level = NULL,
                         extcutoff,
                         K,
                         island_spec,
                         mainland_n,
                         t_hor = NULL) {
  # Function to calculate rates at time = timeval. Returns list with each rate.
  testit::assert(is.numeric(timeval))
  testit::assert(is.numeric(totaltime))
  testit::assert(is.numeric(gam))
  testit::assert(is.numeric(mu))
  testit::assert(is.numeric(laa))
  testit::assert(is.numeric(lac))
  testit::assert(is.numeric(ddmodel_sim))
  testit::assert(is.null(hyper_pars) || is.numeric(hyper_pars))
  testit::assert(is.null(area_pars) || are_area_pars(area_pars))
  testit::assert(is.null(dist_pars) || is.numeric(dist_pars))
  testit::assert(is.null(ext_pars) || is.numeric(ext_pars))
  testit::assert(is.numeric(island_ontogeny))
  testit::assert(is.numeric(extcutoff) || is.null(extcutoff))
  testit::assert(is.numeric(K))
  testit::assert(is.matrix(island_spec) || is.null(island_spec))
  testit::assert(is.numeric(mainland_n))
  testit::assert(is.numeric(t_hor) || is.null(t_hor))
  testit::assert(is.numeric(sea_level))
  immig_rate <- get_immig_rate(timeval = timeval,
                               totaltime = totaltime,
                               gam = gam,
                               ddmodel_sim = ddmodel_sim,
                               hyper_pars = hyper_pars,
                               area_pars = area_pars,
                               dist_pars = dist_pars,
                               island_ontogeny = island_ontogeny,
                               sea_level = sea_level,
                               island_spec = island_spec,
                               K = K,
                               mainland_n = mainland_n)
  testit::assert(is.numeric(immig_rate))
  ext_rate <- get_ext_rate(timeval = timeval,
                           mu = mu,
                           ddmodel_sim = ddmodel_sim,
                           hyper_pars = hyper_pars,
                           area_pars = area_pars,
                           ext_pars = ext_pars,
                           island_ontogeny = island_ontogeny,
                           sea_level = sea_level,
                           extcutoff = extcutoff,
                           island_spec = island_spec,
                           K = K)
  testit::assert(is.numeric(ext_rate))
  ana_rate <- get_ana_rate(laa = laa,
                           hyper_pars = hyper_pars,
                           dist_pars = dist_pars,
                           island_spec = island_spec)
  testit::assert(is.numeric(ana_rate))
  clado_rate <- get_clado_rate(timeval = timeval,
                               lac = lac,
                               ddmodel_sim = ddmodel_sim,
                               hyper_pars = hyper_pars,
                               area_pars = area_pars,
                               dist_pars = dist_pars,
                               island_ontogeny = island_ontogeny,
                               sea_level = sea_level,
                               island_spec = island_spec,
                               K = K)
  testit::assert(is.numeric(clado_rate))
  if (island_ontogeny == 0 && sea_level == 0) {
    immig_rate_max <- immig_rate
    testit::assert(is.numeric(immig_rate_max))
    ext_rate_max <- ext_rate
    testit::assert(is.numeric(ext_rate_max))
    clado_rate_max <- clado_rate
    testit::assert(is.numeric(clado_rate_max))
  } else if (t_hor > timeval) {
    ext_rate_max <- ext_rate
    testit::assert(is.numeric(ext_rate_max))
    immig_rate_max <- get_immig_rate(timeval = t_hor,
                                     totaltime = totaltime,
                                     gam = gam,
                                     ddmodel_sim = ddmodel_sim,
                                     mainland_n = mainland_n,
                                     hyper_pars = hyper_pars,
                                     area_pars = area_pars,
                                     dist_pars = dist_pars,
                                     island_ontogeny = island_ontogeny,
                                     sea_level = sea_level,
                                     island_spec = island_spec,
                                     K = K)
    testit::assert(is.numeric(immig_rate_max))
    clado_rate_max <- get_clado_rate(timeval = t_hor,
                                     lac = lac,
                                     ddmodel_sim = ddmodel_sim,
                                     hyper_pars = hyper_pars,
                                     area_pars = area_pars,
                                     dist_pars = dist_pars,
                                     island_ontogeny = island_ontogeny,
                                     sea_level = sea_level,
                                     island_spec = island_spec,
                                     K = K)
    testit::assert(is.numeric(clado_rate_max))
  } else {
    ext_rate_max <- get_ext_rate(timeval = t_hor,
                                 mu = mu,
                                 ddmodel_sim = ddmodel_sim,
                                 hyper_pars = hyper_pars,
                                 area_pars = area_pars,
                                 ext_pars = ext_pars,
                                 island_ontogeny = island_ontogeny,
                                 sea_level = sea_level,
                                 extcutoff = extcutoff,
                                 island_spec = island_spec,
                                 K = K)
    testit::assert(is.numeric(ext_rate_max) && ext_rate_max >= 0.0)
    immig_rate_max <- immig_rate
    testit::assert(is.numeric(immig_rate_max))
    clado_rate_max <- clado_rate
    testit::assert(is.numeric(clado_rate_max))
  }
  rates <- create_rates(
    immig_rate = immig_rate,
    ext_rate = ext_rate,
    ana_rate = ana_rate,
    clado_rate = clado_rate,
    ext_rate_max = ext_rate_max,
    immig_rate_max = immig_rate_max,
    clado_rate_max = clado_rate_max
  )
  return(rates)
}

#' Function to describe changes in area through time. Adapted from
#' Valente et al 2014 ProcB
#'
#' @param timeval current time of simulation
#' @param area_pars a named list containing area parameters as
#' created by create_area_pars:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#' }
#' @param island_ontogeny a string describing the type of island ontogeny.
#' Can be \code{NULL}, \code{"beta"} for a beta function describing area
#' through time.
#' @param sea_level a numeric describing the type of sea level.
#'
#' @export
#' @family rates calculation
#' @author Pedro Neves
#' @references Valente, Luis M., Rampal S. Etienne, and Albert B. Phillimore.
#' "The effects of island ontogeny on species diversity and phylogeny."
#' Proceedings of the Royal Society of London B: Biological
#' Sciences 281.1784 (2014): 20133227.
island_area <- function(timeval, area_pars, island_ontogeny, sea_level) {
  testit::assert(are_area_pars(area_pars))
  Tmax <- area_pars$total_island_age
  Amax <- area_pars$max_area
  Topt <- area_pars$proportional_peak_t
  peak <- area_pars$peak_sharpness
  ampl <- area_pars$sea_level_amplitude
  freq <- area_pars$sea_level_frequency
  proptime <- timeval / Tmax
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
    At_sine <- ampl * sin(proptime * angular_freq)
    At <- Amax + At_sine
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
    A_sine <- ampl * sin(proptime * angular_freq)
    At <- A_beta + A_sine
    return(At)
  }
}

#' Function to describe changes in extinction rate through time. From
#' Valente et al 2014 ProcB
#'
#' @param timeval current time of simulation
#' @param mu per capita extinction rate in no ontogeny model
#' @param ddmodel_sim A numeric determining which parameters are diversity-
#' dependent.
#' @param area_pars a named list containing area parameters as created by
#' create_area_pars:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#' }
#' @param ext_pars a numeric vector:
#' \itemize{
#'   \item{[1]: minimum extinction when area is at peak}
#'   \item{[2]: extinction rate when current area is 0.10 of maximum area}
#' }
#' @param island_ontogeny a string describing the type of island ontogeny.
#' Can be \code{NULL},
#' \code{"beta"} for a beta function describing area through time.
#' @param extcutoff cutoff for extinction rate preventing it from being too
#' large and slowing down simulation. Default is 1100
#' @param island_spec matrix containing state of system
#' @param sea_level a numeric describing sea level can be \code{NULL}
#' @param K carrying capacity
#' @param hyper_pars A numeric vector for hyperparameters for the rate
#' calculations, \code{hyper_pars[1]} is d_0 the scaling parameter for
#' exponent for calculating cladogenesis rate, \code{hyper_pars[2]}
#' is x the exponent for calculating extinction rate,
#' \code{hyper_pars[3]} is alpha the exponent for calculating the
#' immigration rate, \code{hyper_pars[4]} is beta the exponent for
#' calculating the anagenesis rate.
#'
#' @export
#' @seealso Does the same as \link{DAISIE_calc_clade_ext_rate}
#' @family rates calculation
#' @references Valente, Luis M., Rampal S. Etienne, and Albert B. Phillimore.
#' "The effects of island ontogeny on species diversity and phylogeny."
#' Proceedings of the Royal Society of London B: Biological Sciences 281.1784 (2014): 20133227.
#' @author Pedro Neves
get_ext_rate <- function(timeval,
                         mu,
                         ddmodel_sim,
                         hyper_pars,
                         area_pars,
                         ext_pars,
                         island_ontogeny,
                         sea_level,
                         extcutoff = 1100,
                         island_spec,
                         K) {
  testit::assert(is.numeric(island_ontogeny))
  testit::assert(is.numeric(sea_level))
  if (is.matrix(island_spec) || is.null(island_spec)) {
    N <- length(island_spec[, 1])
  } else if (is.numeric(island_spec)) {
    N <- island_spec
  }
  if (island_ontogeny == 0 && sea_level == 0) {
    if (ddmodel_sim == 0 || ddmodel_sim == 1 || ddmodel_sim == 11) {
      if (is.null(hyper_pars)) {
        ext_rate <- mu * N
      } else {
        X <- log(ext_pars[1] / ext_pars[2]) / log(0.1)
        A <- area_pars[1]
        ext_rate <- ext_pars[1] / ((A / area_pars$max_area) ^ X)
        ext_rate <- ext_rate * N
      }
    }
  }
  if (island_ontogeny != 0 || sea_level != 0) {
    X <- log(ext_pars[1] / ext_pars[2]) / log(0.1)
    ext_rate <-
      ext_pars[1] / ((island_area(timeval,
                               area_pars,
                               island_ontogeny,
                               sea_level) /
                     area_pars$max_area) ^ X)
    ext_rate[which(ext_rate > extcutoff)] <- extcutoff
    ext_rate <- ext_rate * N
  }
  testit::assert(is.numeric(ext_rate))
  testit::assert(ext_rate >= 0)
  return(ext_rate)
}

#' Calculate anagenesis rate
#' @description Internal function.
#' Calculates the anagenesis rate given the current number of
#' immigrant species and the per capita rate.
#'
#' @param laa per capita anagenesis rate
#' @param island_spec matrix with current state of system
#' @param hyper_pars A numeric vector for hyperparameters for the rate
#' calculations, \code{hyper_pars[1]} is d_0 the scaling parameter for
#' exponent for calculating cladogenesis rate, \code{hyper_pars[2]}
#' is x the exponent for calculating extinction rate,
#' \code{hyper_pars[3]} is alpha the exponent for calculating the
#' immigration rate, \code{hyper_pars[4]} is beta the exponent for
#' calculating the anagenesis rate.
#' @param dist_pars a numeric for the distance from the mainland.
#'
#' @seealso Does the same as \link{DAISIE_calc_clade_ana_rate}
#' @family rates calculation
#' @author Pedro Neves
get_ana_rate <- function(laa,
                         hyper_pars,
                         dist_pars,
                         island_spec) {
  if (is.null(hyper_pars)) {
  ana_rate <- laa * length(which(island_spec[, 4] == "I"))
  } else {
    dist <- dist_pars[1]
    beta <- hyper_pars[4]
    ana_rate <- laa * length(which(island_spec[, 4] == "I")) * dist ^ beta
  }
  return(ana_rate)
}

#' Calculate cladogenesis rate
#' @description Internal function.
#' Calculates the cladogenesis rate given the current number of
#' species in the system, the carrying capacity and the per capita cladogenesis
#' rate
#'
#' @param timeval current time of simulation
#' @param lac per capita cladogenesis rate
#' @param ddmodel_sim A numeric determining which parameters are diversity-
#' dependent.
#' @param area_pars a named list containing area parameters as created by create_area_pars:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#' }
#' @param island_ontogeny a numeric describing the type of island ontogeny.
#' Can be \code{NULL}, \code{1} for a beta function describing
#' area through time.
#' @param island_spec matrix with current state of system
#' @param sea_level a numeric describing sea level can be \code{NULL}
#' @param K carrying capacity
#' @param hyper_pars A numeric vector for hyperparameters for the rate
#' calculations, \code{hyper_pars[1]} is d_0 the scaling parameter for
#' exponent for calculating cladogenesis rate, \code{hyper_pars[2]}
#' is x the exponent for calculating extinction rate,
#' \code{hyper_pars[3]} is alpha the exponent for calculating the
#' immigration rate, \code{hyper_pars[4]} is beta the exponent for
#' calculating the anagenesis rate.
#' @param dist_pars a numeric for the distance from the mainland.
#'
#' @export
#' @seealso Does the same as \link{DAISIE_calc_clade_clado_rate}
#' @author Pedro Neves
#' @references Valente, Luis M., Rampal S. Etienne, and Albert B. Phillimore.
#' "The effects of island ontogeny on species diversity and phylogeny."
#' Proceedings of the Royal Society of London B: Biological Sciences 281.1784 (2014): 20133227.
get_clado_rate <- function(timeval,
                           lac,
                           ddmodel_sim = 11,
                           hyper_pars,
                           area_pars,
                           dist_pars,
                           island_ontogeny,
                           sea_level,
                           island_spec,
                           K) {
  # Make function accept island_spec matrix or numeric
  if (is.matrix(island_spec) || is.null(island_spec)) {
    N <- length(island_spec[, 1])
  } else if (is.numeric(island_spec)) {
    N <- island_spec
  }
  # No ontogeny scenario
    testit::assert(is.numeric(island_ontogeny))
    if (island_ontogeny == 0 && sea_level == 0) {
      if (ddmodel_sim == 0) {
        if (is.null(hyper_pars)) {
          clado_rate <- lac * N
        } else {
          d_0 <- hyper_pars[1]
          D <- dist_pars[1]
          clado_rate <- lac * N * A ^ d_0 * log(D)
        }
        testit::assert(is.numeric(clado_rate))
        testit::assert(clado_rate >= 0)
        return(clado_rate)
      }
      if (ddmodel_sim == 1 || ddmodel_sim == 11) {
        clado_rate <- max(c(N * lac * (1 - N / K), 0), na.rm = T)
        return(clado_rate)
      }
    # Ontogeny scenario
    }
    if (island_ontogeny != 0 || sea_level != 0) {
    clado_rate <- max(c(
      N * lac * island_area(timeval, area_pars, island_ontogeny, sea_level) *
        (1 - N / (island_area(
          timeval,
          area_pars,
          island_ontogeny,
          sea_level) * K)), 0), na.rm = T)
  }
  testit::assert(clado_rate >= 0)
  testit::assert(is.numeric(clado_rate))
  return(clado_rate)
}
#' Calculate immigration rate
#' @description Internal function.
#' Calculates the immigration rate given the current number of
#' species in the system, the carrying capacity
#'
#' @param timeval current time of simulation
#' @param totaltime total time of simulation
#' @param gam per capita immigration rate
#' @param ddmodel_sim A numeric determining which parameters are diversity-
#' dependent.
#' @param area_pars a named list containing area parameters as created by create_area_pars:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#' }
#' @param island_spec matrix with current state of system
#' @param K carrying capacity
#' @param sea_level a numeric describing sea level can be \code{NULL}
#' @param mainland_n total number of species present in the mainland
#' @param hyper_pars A numeric vector for hyperparameters for the rate
#' calculations, \code{hyper_pars[1]} is d_0 the scaling parameter for
#' exponent for calculating cladogenesis rate, \code{hyper_pars[2]}
#' is x the exponent for calculating extinction rate,
#' \code{hyper_pars[3]} is alpha the exponent for calculating the
#' immigration rate, \code{hyper_pars[4]} is beta the exponent for
#' calculating the anagenesis rate.
#' @param dist_pars a numeric for the distance from the mainland.
#' @param island_ontogeny a numeric describing the type of island ontogeny.
#' Can be \code{NULL}, \code{1} for a beta function describing
#' area through time.
#'
#' @seealso Does the same as \link{DAISIE_calc_clade_imm_rate}
#' @family rates calculation
#' @author Pedro Neves
#' @references Valente, Luis M., Rampal S. Etienne, and Albert B. Phillimore.
#' "The effects of island ontogeny on species diversity and phylogeny."
#' Proceedings of the Royal Society of London B: Biological Sciences 281.1784 (2014): 20133227.
get_immig_rate <- function(timeval,
                           totaltime,
                           gam,
                           ddmodel_sim = 11,
                           hyper_pars,
                           area_pars,
                           dist_pars,
                           island_ontogeny,
                           sea_level,
                           island_spec,
                           K,
                           mainland_n) {
  N <- length(island_spec[, 1])
  testit::assert(is.numeric(island_ontogeny))
  if (island_ontogeny == 0 && sea_level == 0) {
    if (ddmodel_sim == 0 || ddmodel_sim == 1) {
      if (is.null(hyper_pars)) {
        immig_rate <- gam * mainland_n
      } else {
        dist <- dist_pars[1]
        alpha <- hyper_pars[3]
        immig_rate <- (gam * dist ^ -alpha) / mainland_n
      }
    }
    if (ddmodel_sim == 11) {
      if (is.null(hyper_pars)) {
        immig_rate <- max(c(mainland_n * gam * (1 - N / K), 0), na.rm = T)
      } else {
        dist <- dist_pars[1]
        alpha <- hyper_pars[3]
        immig_rate <- ((gam * dist ^ -alpha) / mainland_n)
      }
    }
  }
  if (island_ontogeny != 0 || sea_level != 0) {
    immig_rate <- max(c(mainland_n * gam * (1 - N / (
      island_area(timeval,
                  area_pars,
                  island_ontogeny,
                  sea_level) * K)), 0), na.rm = T)
  }
  testit::assert(is.numeric(immig_rate))
  testit::assert(immig_rate >= 0)
  return(immig_rate)
}

#' Function to calculate and update horizon for maximum extinction rate
#' @description Internal function.
#' Calculates when the next horizon for maximum extinction will be in the
#' simulation
#'
#' @param timeval current time of simulation
#' @param totaltime total time of simulation
#' @param area_pars a named list containing area parameters as created by create_area_pars:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#' }
#' @param ext_multiplier reduces or increases distance of horizon to current
#' simulation time
#' @param island_ontogeny a string describing the type of island ontogeny.
#'  Can be \code{NULL}, \code{"beta"} for a beta function
#'   describing area through time.
#' @param ext effective extinction rate at timeval
#' @param t_hor time of horizon for max extinction
#' @param sea_level a numeric describing the type of sea level.
#'
#' @family rates calculation
#' @author Pedro Neves
get_t_hor <- function(timeval,
                      totaltime,
                      area_pars,
                      ext,
                      ext_multiplier,
                      island_ontogeny,
                      sea_level,
                      t_hor) {
  testit::assert(is.null(area_pars) || are_area_pars(area_pars))
  # Function calculates where the horizon for max(ext_rate) is.
  if (island_ontogeny == 0 & sea_level == 0) {
    testit::assert(totaltime > 0.0)
    t_hor <- totaltime
  }
  if (island_ontogeny == 1 & sea_level == 0) {
    if (is.null(t_hor)) {
      testit::assert(are_area_pars(area_pars))
      t_hor <- area_pars$proportional_peak_t * area_pars$total_island_age
    } else if (timeval >= t_hor) {
      # t_hor should dynamically be adjusted depending on parameter values.
      # Certain parameter combinations will always make it be > totaltime at
      # first calculation, slowing down the simulations
      t_hor <- t_hor + t_hor / 6 + ext_multiplier * (totaltime - timeval) * ext
      t_hor <- min(totaltime, t_hor)
    }
    testit::assert(t_hor > 0.0)
  }
  if (island_ontogeny == 0 & sea_level == 1) {
    if (is.null(t_hor)) {
      t_hor <- 1 / (area_pars$sea_level_frequency * 4)
    } else if (timeval >= t_hor &
               timeval < 1 / area_pars$sea_level_frequency + t_hor) {
      t_hor <- (1 / area_pars$sea_level_frequency) + t_hor
      t_hor <- min(totaltime, t_hor)
    }
  }
  if (island_ontogeny == 1 & sea_level == 1) {
    if (is.null(t_hor)) {
      max <- optimize(
      f = DAISIE::island_area,
      interval = c(0, totaltime),
      area_pars = area_pars,
      island_ontogeny = 1,
      sea_level = 0,
      maximum = TRUE
    )
      t_hor <- max$objective
    } else if (timeval >= t_hor) {
      t_hor <- t_hor + t_hor / 6 + ext_multiplier * (totaltime - timeval) * ext
    }
    }
  return(t_hor)
}

#' Calculates when the next timestep will be.
#'
#' @param rates list of numeric with probabilities of each event
#' @param timeval current time of simulation
#'
#' @return named list with numeric vector containing the time of the next
#' timestep and the change in time.
#' @author Pedro Neves
calc_next_timeval <- function(rates, timeval) {
  # Calculates when next event will happen
  testit::assert(are_rates(rates))
  testit::assert(timeval >= 0)
  totalrate <- rates$immig_rate_max + rates$ana_rate + rates$clado_rate_max +
    rates$ext_rate_max
  dt <- stats::rexp(1, totalrate)
  timeval <- timeval + dt
  return(list(timeval = timeval, dt = dt))
}


#' Calculate the clade-wide extinction rate
#'
#' @param ps_ext_rate per species extinction rate
#' @param n_species number of species in that clade
#'
#' @return the clade's extinction rate
#' @author Richel J.C. Bilderbeek
#' @examples
#'   testit::assert(
#'     DAISIE_calc_clade_ext_rate(
#'       ps_ext_rate = 0.2,
#'       n_species = 4
#'     ) == 0.8
#'   )
#' @export
DAISIE_calc_clade_ext_rate <- function(ps_ext_rate, n_species) {
  testit::assert(ps_ext_rate >= 0.0)
  testit::assert(n_species >= 0)
  ps_ext_rate * n_species
}

#' Calculate the clade-wide effective anagenesis rate.
#' With 'effective', this means that if an immigrant
#' undergoes anagenesis, it will become a new species.
#' Would such a species undergo anagenesis again, no net new
#' species is created; the species only gets renamed
#'
#' @param ps_ana_rate per species anagensis rate
#' @param n_immigrants number of immigrants in that clade
#'
#' @return the clade's effective anagenesis rate
#' @author Richel J.C. Bilderbeek
#' @examples
#'   testit::assert(
#'     DAISIE_calc_clade_ana_rate(
#'       ps_ana_rate = 0.3,
#'       n_immigrants = 5
#'     ) == 1.5
#'   )
#' @export
DAISIE_calc_clade_ana_rate <- function(ps_ana_rate, n_immigrants) {
  testit::assert(ps_ana_rate >= 0.0)
  testit::assert(n_immigrants >= 0)
  ps_ana_rate * n_immigrants
}

#' Calculate the clade-wide cladogenesis rate.
#'
#' @param ps_clado_rate per species cladogenesis rate
#' @param n_species number of species in that clade
#' @param carr_cap carrying capacity, number of species this clade will
#'   grow to
#' @return the clade's cladogenesis rate, which is at least zero. This
#'   rate will be zero if there are more species than the carrying capacity
#'   allows for
#' @note For clade-specific carrying capacity,
#'   each clade is simulated seperately in \code{\link{DAISIE_sim}}
#' @author Richel J.C. Bilderbeek
#' @examples
#'   testit::assert(
#'     DAISIE_calc_clade_clado_rate(
#'       ps_clado_rate = 0.2,
#'       n_species = 5,
#'       carr_cap = 10
#'     ) == 0.5
#'   )
#'   testit::assert(
#'     DAISIE_calc_clade_clado_rate(
#'       ps_clado_rate = 0.2,
#'       n_species = 2,
#'       carr_cap = 1
#'     ) == 0.0
#'   )
#' @export
DAISIE_calc_clade_clado_rate <- function(ps_clado_rate, n_species, carr_cap) {
  testit::assert(ps_clado_rate >= 0.0)
  testit::assert(n_species >= 0)
  testit::assert(carr_cap >= 0)
  return(max(
    0.0,
    n_species * ps_clado_rate * (1.0 - (n_species / carr_cap))
  ))
}

#' Calculate the clade-wide immigration rate.
#'
#' @param ps_imm_rate per species immigration rate
#' @param n_island_species number of species in that clade on the island
#' @param n_mainland_species number of species in that clade on the mainland
#' @param carr_cap carrying capacity, number of species this clade will
#'   grow to
#' @return the clade's immigration rate, which is at least zero. This
#'   rate will be zero if there are more species than the carrying capacity
#'   allows for
#' @author Richel J.C. Bilderbeek
#' @examples
#'   testit::assert(
#'     DAISIE_calc_clade_imm_rate(
#'       ps_imm_rate = 0.1,
#'       n_island_species = 5,
#'       n_mainland_species = 2,
#'       carr_cap = 10
#'     ) == 0.1
#'   )
#'   testit::assert(
#'     DAISIE_calc_clade_imm_rate(
#'       ps_imm_rate = 0.1,
#'       n_island_species = 5,
#'       n_mainland_species = 2,
#'       carr_cap = 1
#'     ) == 0.0
#'   )
#' @export
DAISIE_calc_clade_imm_rate <- function(
  ps_imm_rate,
  n_island_species,
  n_mainland_species,
  carr_cap
) {
  testit::assert(ps_imm_rate >= 0.0)
  testit::assert(n_island_species >= 0)
  testit::assert(n_mainland_species >= 0)
  testit::assert(carr_cap >= 0)
  return(max(
    0.0,
    n_mainland_species * ps_imm_rate * (1.0 - (n_island_species / carr_cap))
  ))
}
