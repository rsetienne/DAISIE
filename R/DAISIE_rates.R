#' Calculates algorithm rates
#' @description Internal function that updates the all the rates and
#' max extinction horizon at time t.
#' @family rate calculations
#'
#' @param timeval A numeric with the current time of simulation
#' @param totaltime A numeric with the total time of simulation
#' @param gam A numeric with the per capita immigration rate
#' @param laa A numeric with the per capita anagenesis rate
#' @param lac A numeric with the per capita cladogenesis rate
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
#' @param ext_pars A numeric vector:
#' \itemize{
#'   \item{[1]: minimum extinction when area is at peak}
#'   \item{[2]: extinction rate when current area is 0.10 of maximum area}
#' }
#' @param island_ontogeny A string describing the type of island ontogeny.
#' Can be \code{NULL},
#' \code{"beta"} for a beta function describing area through time.
#' @param extcutoff A numeric with the cutoff for extinction rate
#' preventing it from being too large and slowing down simulation.
#' Should be big.
#' @param K A numeric with K (clade-specific carrying capacity)
#' @param mainland_n A numeirc with the total number of species present
#' in the mainland
#' @param hyper_pars A numeric vector for hyperparameters for the rate
#' calculations:
#' \itemize{
#' \item{[1]: is d_0 the scaling parameter for exponent for calculating
#' cladogenesis rate}
#' \item{[2]: is x the exponent for calculating extinction rate}
#' \item{[3]: is alpha, the exponent for calculating the immigration rate}
#' \item{[4]: is beta the exponent for calculating the anagenesis rate.}
#' }
#' @param sea_level a numeric describing the type of sea level.
#' @param num_spec a numeric with the current number of species.
#' @param num_immigrants a numeric with the current number of non-endemic
#' species (a.k.a non-endemic species).
#' @param mu extinction rate
#' @param dist_pars a numeric for the distance from the mainland.
#'
#' @seealso \code{\link{update_max_rates}}
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
                         dist_pars = NULL,
                         ext_pars = NULL,
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
  testit::assert(are_dist_pars(dist_pars))
  testit::assert(is.null(ext_pars) || is.numeric(ext_pars))
  testit::assert(is.numeric(island_ontogeny))
  testit::assert(is.numeric(extcutoff) || is.null(extcutoff))
  testit::assert(is.numeric(K))
  testit::assert(is.numeric(num_spec) || is.null(num_spec))
  testit::assert(is.numeric(num_immigrants) || is.null(num_immigrants))
  testit::assert(is.numeric(mainland_n))
  testit::assert(is.numeric(sea_level))
  immig_rate <- get_immig_rate(
    timeval = timeval,
    totaltime = totaltime,
    gam = gam,
    hyper_pars = hyper_pars,
    area_pars = area_pars,
    dist_pars = dist_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    num_spec = num_spec,
    K = K,
    mainland_n = mainland_n
  )
  testit::assert(is.numeric(immig_rate))
  ext_rate <- get_ext_rate(
    timeval = timeval,
    mu = mu,
    hyper_pars = hyper_pars,
    area_pars = area_pars,
    ext_pars = ext_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    extcutoff = extcutoff,
    num_spec = num_spec,
    K = K
  )
  testit::assert(is.numeric(ext_rate))
  ana_rate <- get_ana_rate(
    laa = laa,
    hyper_pars = hyper_pars,
    dist_pars = dist_pars,
    num_immigrants = num_immigrants
  )
  testit::assert(is.numeric(ana_rate))
  clado_rate <- get_clado_rate(
    timeval = timeval,
    lac = lac,
    hyper_pars = hyper_pars,
    area_pars = area_pars,
    dist_pars = dist_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    num_spec = num_spec,
    K = K
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

#' Function to describe changes in area through time. Adapted from
#' Valente et al 2014 ProcB
#'
#' @param timeval current time of simulation
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
#' @param island_ontogeny a string describing the type of island ontogeny.
#' Can be \code{NULL}, \code{"beta"} for a beta function describing area
#' through time.
#' @param sea_level a numeric describing the type of sea level.
#'
#' @export
#' @family rate calculations
#' @author Pedro Neves, Joshua Lambert
#' @references Valente, Luis M., Rampal S. Etienne, and Albert B. Phillimore.
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
    r_zero <- sqrt((Amax * cos(theta)) / pi)
    h_zero <- tan(theta) * r_zero
    At <- pi * ((h_zero - delta_sl) ^ 2) * cos(theta) / (sin(theta)^2)
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
    r_zero <- sqrt((A_beta * cos(theta)) / pi)
    h_zero <- tan(theta) * r_zero
    At <- pi * ((h_zero - delta_sl) ^ 2) * cos(theta) / (sin(theta)^2)
    return(At)
  }
}

#' Function to describe changes in extinction rate through time. From
#' Valente et al 2014 ProcB
#'
#' @param timeval current time of simulation
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
#' @param sea_level a numeric describing sea level can be \code{NULL}
#' @param K carrying capacity
#' @param hyper_pars A numeric vector for hyperparameters for the rate
#' calculations:
#' \itemize{
#' \item{[1]: is d_0 the scaling parameter for exponent for calculating
#' cladogenesis rate}
#' \item{[2]: is x the exponent for calculating extinction rate}
#' \item{[3]: is alpha, the exponent for calculating the immigration rate}
#' \item{[4]: is beta the exponent for calculating the anagenesis rate.}
#' }
#' @param num_spec a numeric with the current number of species
#' @param mu extinction rate
#'
#' @export
#' @family rate calculations
#' @references Valente, Luis M., Rampal S. Etienne, and Albert B. Phillimore.
#' "The effects of island ontogeny on species diversity and phylogeny."
#' Proceedings of the Royal Society of London B: Biological Sciences 281.1784
#' (2014): 20133227.
#' @author Pedro Neves, Joshua Lambert
get_ext_rate <- function(timeval,
                         mu,
                         ext_pars,
                         hyper_pars,
                         area_pars,
                         island_ontogeny,
                         sea_level = 0,
                         extcutoff = 1000,
                         num_spec,
                         K) {
  testit::assert(is.numeric(island_ontogeny))
  testit::assert(is.numeric(sea_level))
  A <- island_area(
    timeval,
    area_pars,
    island_ontogeny,
    sea_level
  )
  if (island_ontogeny == 1 || sea_level == 1) {
    x <- log(ext_pars[1] / ext_pars[2]) / log(0.1)
  } else {
    x <- hyper_pars$x
    ext_pars[1] <- mu # Constant rate case
  }
  ext_rate <- ext_pars[1] / ((A / area_pars$max_area) ^ x)
  ext_rate <- ext_rate * num_spec
  ext_rate <- min(ext_rate, extcutoff, na.rm = TRUE)
  if (num_spec == 0) {
    ext_rate <- 0
  }
  testit::assert(ext_rate >= 0)
  return(ext_rate)
}

#' Calculate anagenesis rate
#' @description Internal function.
#' Calculates the anagenesis rate given the current number of
#' immigrant species and the per capita rate.
#'
#' @param laa per capita anagenesis rate
#' @param hyper_pars A numeric vector for hyperparameters for the rate
#' calculations:
#' \itemize{
#' \item{[1]: is d_0 the scaling parameter for exponent for calculating
#' cladogenesis rate}
#' \item{[2]: is x the exponent for calculating extinction rate}
#' \item{[3]: is alpha, the exponent for calculating the immigration rate}
#' \item{[4]: is beta the exponent for calculating the anagenesis rate.}
#' }
#' @param dist_pars a numeric for the distance from the mainland.
#' @param num_immigrants a numeric with the current number of non-endemic
#' species (a.k.a non-endemic species).
#'
#' @export
#' @family rate calculations
#' @author Pedro Neves, Joshua Lambert
get_ana_rate <- function(laa,
                         hyper_pars,
                         dist_pars,
                         num_immigrants) {
  D <- dist_pars$D
  beta <- hyper_pars$beta
  ana_rate <- laa * num_immigrants * D ^ beta

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
#' @param timeval current time of simulation
#' @param lac per capita cladogenesis rate
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
#' @param sea_level a numeric describing sea level can be \code{NULL}
#' @param K carrying capacity
#' @param hyper_pars A numeric vector for hyperparameters for the rate
#' calculations:
#' \itemize{
#' \item{[1]: is d_0 the scaling parameter for exponent for calculating
#' cladogenesis rate}
#' \item{[2]: is x the exponent for calculating extinction rate}
#' \item{[3]: is alpha, the exponent for calculating the immigration rate}
#' \item{[4]: is beta the exponent for calculating the anagenesis rate.}
#' }
#' @param dist_pars a numeric for the distance from the mainland.
#' @param num_spec a numeric with the ccurrent number of species.
#'
#' @export
#' @author Pedro Neves
#' @references Valente, Luis M., Rampal S. Etienne, and Albert B. Phillimore.
#' "The effects of island ontogeny on species diversity and phylogeny."
#' Proceedings of the Royal Society of London B: Biological Sciences 281.1784 (2014): 20133227.
get_clado_rate <- function(timeval,
                           lac,
                           hyper_pars = NULL,
                           area_pars,
                           dist_pars,
                           island_ontogeny,
                           sea_level = 0,
                           num_spec,
                           K) {
  testit::assert(is.numeric(island_ontogeny))
  testit::assert(is.numeric(sea_level))
  testit::assert(are_hyper_pars(hyper_pars))
  testit::assert(are_dist_pars(dist_pars))

  A <- DAISIE::island_area(
    timeval = timeval,
    area_pars = area_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level
  )
  d_0 <- hyper_pars$d_0
  D <- dist_pars$D

  clado_rate <- max(
    0, lac * num_spec * A ^ (d_0 * log(D)) * (1 - num_spec / (K * A)),
    na.rm = TRUE
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
#' @param timeval current time of simulation
#' @param totaltime total time of simulation
#' @param gam per capita immigration rate
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
#' @param K carrying capacity
#' @param sea_level a numeric describing sea level can be \code{NULL}
#' @param mainland_n total number of species present in the mainland
#' @param hyper_pars A numeric vector for hyperparameters for the rate
#' calculations:
#' \itemize{
#' \item{[1]: is d_0 the scaling parameter for exponent for calculating
#' cladogenesis rate}
#' \item{[2]: is x the exponent for calculating extinction rate}
#' \item{[3]: is alpha, the exponent for calculating the immigration rate}
#' \item{[4]: is beta the exponent for calculating the anagenesis rate.}
#' }
#' @param dist_pars a numeric for the distance from the mainland.
#' @param island_ontogeny a numeric describing the type of island ontogeny.
#' Can be \code{NULL}, \code{1} for a beta function describing
#' area through time.
#' @param num_spec a numeric with the current number of species.
#'
#' @export
#' @family rate calculations
#' @author Pedro Neves, Joshua Lambert
#' @references Valente, Luis M., Rampal S. Etienne, and Albert B. Phillimore.
#' "The effects of island ontogeny on species diversity and phylogeny."
#' Proceedings of the Royal Society of London B: Biological Sciences 281.1784 (2014): 20133227.
get_immig_rate <- function(timeval,
                           totaltime,
                           gam,
                           hyper_pars,
                           area_pars,
                           dist_pars,
                           island_ontogeny,
                           sea_level,
                           num_spec,
                           K,
                           mainland_n) {
  testit::assert(is.numeric(island_ontogeny))
  testit::assert(is.numeric(sea_level))

  D <- dist_pars$D
  alpha <- hyper_pars$alpha
  A <- DAISIE::island_area(
    timeval = timeval,
    area_pars = area_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level
  )
  immig_rate <- max(c(mainland_n * gam * D^-alpha  * (1 - (num_spec / (A * K))),
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


