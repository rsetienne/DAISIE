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
#' @param trait_pars A named list containing diversification rates considering two trait states:
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
                         mainland_n,
                         trait_pars = NULL,
                         island_spec = NULL) {
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
  testit::assert(are_trait_pars(trait_pars))
  testit::assert(is.null(ext_pars) || is.numeric(ext_pars))
  testit::assert(is.numeric(island_ontogeny))
  testit::assert(is.numeric(extcutoff) || is.null(extcutoff))
  testit::assert(is.numeric(K))
  testit::assert(is.numeric(num_spec) || is.null(num_spec))
  testit::assert(is.numeric(num_immigrants) || is.null(num_immigrants))
  testit::assert(is.numeric(mainland_n))
  testit::assert(is.numeric(sea_level))
  
  if (!is.null(trait_pars)) {
    return(
      update_rates_trait(
        timeval = timeval,
        totaltime = totaltime,
        gam = gam, 
        mu = mu, 
        laa = laa, 
        lac = lac,
        hyper_pars = hyper_pars,
        area_pars = area_pars,
        dist_pars = dist_pars,
        ext_pars = ext_pars,
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


update_rates_trait <- function(timeval,
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
                               mainland_n,
                               trait_pars = NULL,
                               island_spec) {
  # Function to calculate rates at time = timeval. Returns list with each rate.
  
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
    mainland_n = mainland_n,
    trait_pars = trait_pars,
    island_spec = island_spec)
  testit::assert(is.list(immig_rate))
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
    K = K,
    trait_pars = trait_pars,
    island_spec = island_spec
  )
  testit::assert(is.list(ext_rate))
  ana_rate <- get_ana_rate(
    laa = laa,
    hyper_pars = hyper_pars,
    dist_pars = dist_pars,
    num_immigrants = num_immigrants,
    trait_pars = trait_pars,
    island_spec = island_spec
  )
  testit::assert(is.list(ana_rate))
  clado_rate <- get_clado_rate(
    timeval = timeval,
    lac = lac,
    hyper_pars = hyper_pars,
    area_pars = area_pars,
    dist_pars = dist_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    num_spec = num_spec,
    K = K,
    trait_pars = trait_pars,
    island_spec = island_spec
  )
  testit::assert(is.list(clado_rate))
  
  testit::assert(!is.null(trait_pars))
  trans_rate <- get_trans_rate(trait_pars = trait_pars,
                               island_spec = island_spec)
  testit::assert(is.list(trans_rate))
  # trait_pars <- create_trait_pars(trans_rate = trans_rate$trans_rate1,
  #                                    immig_rate2 = immig_rate$immig_rate2,
  #                                    ext_rate2 = ext_rate$ext_rate2,
  #                                    ana_rate2 = ana_rate$ana_rate2,
  #                                    clado_rate2 = clado_rate$clado_rate2,
  #                                    trans_rate2 = trans_rate$trans_rate2,
  #                                    M2 = trait_pars$M2)
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
                         K,
                         trait_pars = NULL,
                         island_spec = NULL) {
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
  if(is.null(trait_pars)){
    ext_rate <- ext_pars[1] / ((A / area_pars$max_area) ^ x)
    ext_rate <- ext_rate * num_spec
    ext_rate <- min(ext_rate, extcutoff, na.rm = TRUE)
    if (num_spec == 0) {
      ext_rate <- 0
    }
    testit::assert(ext_rate >= 0)
    return(ext_rate)
  }else{   ##species have two states
    if (is.matrix(island_spec) || is.null(island_spec)) {
      num_spec_trait1 <- length(which(island_spec[, 8] == "1"))
      num_spec_trait2 <- length(which(island_spec[, 8] == "2"))
    } else if (is.numeric(island_spec)) {
      stop("Different trait states cannot be separated,please transform to matrix form.")
    }
    ext_rate1 <- mu * num_spec_trait1
    ext_rate2 <- trait_pars$ext_rate2 * num_spec_trait2
    testit::assert(is.numeric(ext_rate1))
    testit::assert(is.numeric(ext_rate2))
    testit::assert(ext_rate1 >= 0)
    testit::assert(ext_rate2 >= 0)
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
#' @export
#' @family rate calculations
#' @author Pedro Neves, Joshua Lambert
get_ana_rate <- function(laa,
                         hyper_pars,
                         dist_pars,
                         num_immigrants,
                         island_spec = NULL,
                         trait_pars = NULL) {
  D <- dist_pars$D
  beta <- hyper_pars$beta
  if(is.null(trait_pars)){
    ana_rate <- laa * num_immigrants * D ^ beta
    testit::assert(is.numeric(ana_rate))
    testit::assert(ana_rate >= 0)
    return(ana_rate) 
  }else{
    ana_rate1 = laa * length(intersect(which(island_spec[,4] == "I"), 
                                       which(island_spec[,8] == "1"))) * D ^ beta
    ana_rate2 = trait_pars$ana_rate2 * length(intersect(which(island_spec[,4] == "I"), 
                                                   which(island_spec[,8] == "2"))) * D ^ beta
    testit::assert(is.numeric(ana_rate1))
    testit::assert(ana_rate1 >= 0)
    testit::assert(is.numeric(ana_rate2))
    testit::assert(ana_rate2 >= 0)
    ana_list <- list(ana_rate1 = ana_rate1,
                     ana_rate2 = ana_rate2)
    return(ana_list)
  }
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
                           K,
                           trait_pars = NULL,
                           island_spec = NULL) {
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
  if(is.null(trait_pars)){
  clado_rate <- max(
    0, lac * num_spec * A ^ (d_0 * log(D)) * (1 - num_spec / (K * A)),
    na.rm = TRUE
  )
  testit::assert(clado_rate >= 0)
  testit::assert(is.numeric(clado_rate))
  return(clado_rate)
  }else{
    num_spec_trait1 <- length(which(island_spec[, 8] == "1"))
    num_spec_trait2 <- length(which(island_spec[, 8] == "2"))
    clado_rate1 <- max(
      0, lac * num_spec_trait1 * A ^ (d_0 * log(D)) * (1 - num_spec / (K * A)),
      na.rm = TRUE)
    clado_rate2 <- max(
      0, trait_pars$clado_rate2 * num_spec_trait2 * A ^ (d_0 * log(D)) * (1 - num_spec / (K * A)),
      na.rm = TRUE
    )
    testit::assert(clado_rate1 >= 0)
    testit::assert(clado_rate2 >= 0)
    testit::assert(is.numeric(clado_rate1))
    testit::assert(is.numeric(clado_rate2))
    clado_list <- list(clado_rate1 = clado_rate1,
                       clado_rate2 = clado_rate2)
    return(clado_list)
  }
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
                           mainland_n,
                           trait_pars = NULL,
                           island_spec = NULL) {
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
  if(is.null(trait_pars)){
    immig_rate <- max(c(mainland_n * gam * D^-alpha  * (1 - (num_spec / (A * K))),
                        0), na.rm = TRUE)
    testit::assert(is.numeric(immig_rate))
    testit::assert(immig_rate >= 0)
    return(immig_rate) 
  }else{
    mainland_n2 <- trait_pars$M2
    gam2 <- trait_pars$immig_rate2
    immig_rate1 <- max(c(mainland_n * gam * D^-alpha  * (1 - (num_spec / (A * K))),
                        0), na.rm = TRUE)
    immig_rate2 <- max(c(mainland_n2 * gam2 * D^-alpha  * (1 - (num_spec / (A * K))),
                        0), na.rm = TRUE)
    testit::assert(is.numeric(immig_rate1))
    testit::assert(immig_rate1 >= 0)
    testit::assert(is.numeric(immig_rate2))
    testit::assert(immig_rate2 >= 0)
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
#' @family rates calculation
get_trans_rate <- function(trait_pars,
                           island_spec){
  if(is.null(trait_pars)){
    stop("Transition rate only calculate when exists more than one trait state.") #or trans_rate = NULL
  }else{
    # Make function accept island_spec matrix or numeric
    if (is.matrix(island_spec) || is.null(island_spec)) {
      num_spec_trait1 <- length(which(island_spec[, 8] == "1"))
      num_spec_trait2 <- length(which(island_spec[, 8] == "2"))
    } else if (is.numeric(island_spec)) {
      stop("Different trait states cannot be separated,please transform to matrix form.")
    }
    trans_rate1 <- trait_pars$trans_rate * num_spec_trait1
    trans_rate2 <- trait_pars$trans_rate2 * num_spec_trait2
    testit::assert(is.numeric(trans_rate1))
    testit::assert(trans_rate1 >= 0)
    testit::assert(is.numeric(trans_rate2))
    testit::assert(trans_rate2 >= 0)
    trans_list <- list(trans_rate1 = trans_rate1,
                       trans_rate2 = trans_rate2)
    return(trans_list)
  }
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
  
  if(length(max_rates) <= 4){   ## no considering about two trait states
    totalrate <- max_rates[[1]] + max_rates[[2]] + max_rates[[3]] + max_rates[[4]] 
  }else{
    totalrate <- max_rates[[1]] + max_rates[[2]] + max_rates[[3]] + max_rates[[4]] + 
      max_rates[[5]] + max_rates[[6]] + max_rates[[7]] + max_rates[[8]] + 
      max_rates[[9]] + max_rates[[10]] 
  }
  dt <- stats::rexp(1, totalrate)
  timeval <- timeval + dt
  return(list(timeval = timeval, dt = dt))
}

