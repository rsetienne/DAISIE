#' Calculates algorithm rates
#' @description Internal function that updates the all the rates and
#' max extinction horizon at time t.
#' @family rates calculation
#' @param timeval A numeric with the current time of simulation
#' @param totaltime A numeric with the total time of simulation
#' @param gam A numeric with the per capita immigration rate
#' @param mu A numeric with the per capita extinction rate in no ontogeny model
#' @param laa A numeric with the per capita anagenesis rate
#' @param lac A numeric with the per capita cladogenesis rate
#' @param Apars A named list containing area parameters as created by create_area_params:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#' }
#' @param Epars A numeric vector:
#' \itemize{
#'   \item{[1]: minimum extinction when area is at peak}
#'   \item{[2]: extinction rate when current area is 0.10 of maximum area}
#' }
#' @param Tpars A named list containing diversification rates considering two trait states:
#' \itemize{
#'   \item{[1]:A numeric with the per capita transition rate with state1}
#'   \item{[2]:A numeric with the per capita immigration rate with state2}
#'   \item{[3]:A numeric with the per capita extinction rate with state2}
#'   \item{[4]:A numeric with the per capita anagenesis rate with state2}
#'   \item{[5]:A numeric with the per capita cladogenesis rate with state2}
#'   \item{[6]:A numeric with the per capita transition rate with state2}
#'   \item{[7]:A numeric with the number of species with trait state 2 on mainland}
#' }
#' @param island_ontogeny A string describing the type of island ontogeny.
#' Can be \code{NULL},
#' \code{"beta"} for a beta function describing area through time,
#'  or \code{"linear"} for a linear function
#' @param extcutoff A numeric with the cutoff for extinction rate preventing it from being too
#' large and slowing down simulation. Should be big.
#' @param K A numeric with K (clade-specific carrying capacity)
#' @param island_spec A matrix containing state of system
#' @param mainland_n A numeirc with the total number of species present in the mainland
#' @param t_hor A numeric with the time of horizon for max cladogenesis, immigration and minimum extinction
update_rates <- function(timeval, totaltime,
                         gam, mu, laa, lac, 
                         Apars = NULL, Epars = NULL,Tpars = NULL,
                         island_ontogeny = 0,
                         extcutoff,
                         K,
                         island_spec, mainland_n, t_hor = NULL) {
  # Function to calculate rates at time = timeval. Returns list with each rate.
  testit::assert(is.numeric(timeval))
  testit::assert(is.numeric(totaltime))
  testit::assert(is.numeric(gam))
  testit::assert(is.numeric(mu))
  testit::assert(is.numeric(laa))
  testit::assert(is.numeric(lac))
  testit::assert(is.null(Tpars) || are_trait_state_params(Tpars))
  testit::assert(is.null(Apars) || are_area_params(Apars))
  testit::assert(is.null(Epars) || is.numeric(Epars))
  testit::assert(is.numeric(island_ontogeny))
  testit::assert(is.numeric(extcutoff) || is.null(extcutoff))
  testit::assert(is.numeric(K))
  testit::assert(is.matrix(island_spec) || is.null(island_spec))
  testit::assert(is.numeric(mainland_n))
  testit::assert(is.numeric(t_hor) || is.null(t_hor))

  if (!is.null(Tpars)) {
    return(
      update_rates_trait(
        timeval = timeval, 
        totaltime = totaltime,
        gam = gam, mu = mu, laa = laa, lac = lac, 
        Apars = Apars, 
        Epars = Epars,
        Tpars = Tpars,
        island_ontogeny = island_ontogeny,
        extcutoff = extcutoff,
        K = K,
        island_spec = island_spec, 
        mainland_n = mainland_n, 
        t_hor = t_hor
      )
    )
  }

  immig_rate <- get_immig_rate(timeval = timeval,
                               totaltime = totaltime,
                               gam = gam,
                               Apars = Apars,
                               Tpars = Tpars,
                               island_ontogeny = island_ontogeny,
                               island_spec = island_spec,
                               K = K,
                               mainland_n = mainland_n)

  ext_rate <- get_ext_rate(timeval = timeval,
                           mu = mu,
                           Tpars = Tpars,
                           Apars = Apars,
                           Epars = Epars,
                           island_ontogeny = island_ontogeny,
                           extcutoff = extcutoff,
                           island_spec = island_spec,
                           K = K)

  ana_rate <- get_ana_rate(laa = laa,
                           Tpars = Tpars,
                           island_spec = island_spec)


  clado_rate <- get_clado_rate(timeval = timeval,
                               lac = lac,
                               Apars = Apars,
                               Tpars = Tpars,
                               island_ontogeny = island_ontogeny,
                               island_spec = island_spec,
                               K = K)

  if ((island_ontogeny) == 0) {
    immig_rate_max <- immig_rate
    ext_rate_max <- ext_rate
    clado_rate_max <- clado_rate
    testit::assert(is.numeric(immig_rate_max))
    testit::assert(is.numeric(ext_rate_max))
    testit::assert(is.numeric(clado_rate_max))
  } else if ((Apars$proportional_peak_t * Apars$total_island_age) > timeval) {
    ext_rate_max <- ext_rate
    testit::assert(is.numeric(ext_rate_max))
    testit::assert(is.null(Tpars))
    immig_rate_max <- get_immig_rate(timeval = Apars$proportional_peak_t * Apars$total_island_age,
                                     totaltime = totaltime,
                                     gam = gam,
                                     mainland_n = mainland_n,
                                     Apars = Apars,
                                     Tpars = Tpars,
                                     island_ontogeny = island_ontogeny,
                                     island_spec = island_spec,
                                     K = K)
    testit::assert(is.numeric(immig_rate_max))
    
    clado_rate_max <- get_clado_rate(timeval = Apars$proportional_peak_t * Apars$total_island_age, # SHOULD BE GENERALIZED
                                     lac = lac,
                                     Apars = Apars,
                                     Tpars = Tpars,
                                     island_ontogeny = island_ontogeny,
                                     island_spec = island_spec,
                                     K = K)
    testit::assert(is.numeric(clado_rate_max))
    
  } else {
    # Ontogeny, max rate is t_hor, which in this case is totaltime (from hor)
    ext_rate_max <- get_ext_rate(timeval = t_hor,
                                 mu = mu,
                                 Apars = Apars,
                                 Epars = Epars,
                                 Tpars = Tpars,
                                 island_ontogeny = island_ontogeny,
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
    clado_rate_max = clado_rate_max,
    Tpars = Tpars)

  return(rates)
}


update_rates_trait <- function(timeval, totaltime,
                         gam, mu, laa, lac, 
                         Apars = NULL, 
                         Epars = NULL,
                         Tpars = NULL,
                         island_ontogeny = 0,
                         extcutoff,
                         K,
                         island_spec, mainland_n, t_hor = NULL) {
  # Function to calculate rates at time = timeval. Returns list with each rate.
  testit::assert(is.numeric(timeval))
  testit::assert(is.numeric(totaltime))
  testit::assert(is.numeric(gam))
  testit::assert(is.numeric(mu))
  testit::assert(is.numeric(laa))
  testit::assert(is.numeric(lac))
  testit::assert(is.null(Tpars) || are_trait_state_params(Tpars))
  testit::assert(is.null(Apars) || are_area_params(Apars))
  testit::assert(is.null(Epars) || is.numeric(Epars))
  testit::assert(is.numeric(island_ontogeny))
  testit::assert(is.numeric(extcutoff) || is.null(extcutoff))
  testit::assert(is.numeric(K))
  testit::assert(is.matrix(island_spec) || is.null(island_spec))
  testit::assert(is.numeric(mainland_n))
  testit::assert(is.numeric(t_hor) || is.null(t_hor))
  
  immig_rate <- get_immig_rate(timeval = timeval,
                               totaltime = totaltime,
                               gam = gam,
                               Apars = Apars,
                               Tpars = Tpars,
                               island_ontogeny = island_ontogeny,
                               island_spec = island_spec,
                               K = K,
                               mainland_n = mainland_n)
  
  ext_rate <- get_ext_rate(timeval = timeval,
                           mu = mu,
                           Tpars = Tpars,
                           Apars = Apars,
                           Epars = Epars,
                           island_ontogeny = island_ontogeny,
                           extcutoff = extcutoff,
                           island_spec = island_spec,
                           K = K)
  
  ana_rate <- get_ana_rate(laa = laa,
                           Tpars = Tpars,
                           island_spec = island_spec)
  
  
  clado_rate <- get_clado_rate(timeval = timeval,
                               lac = lac,
                               Apars = Apars,
                               Tpars = Tpars,
                               island_ontogeny = island_ontogeny,
                               island_spec = island_spec,
                               K = K)
  if ((island_ontogeny) == 0) {
    immig_rate_max <- max(unlist(immig_rate))
    ext_rate_max <- max(unlist(ext_rate))
    clado_rate_max <- max(unlist(clado_rate))
    testit::assert(is.numeric(immig_rate_max))
    testit::assert(is.numeric(ext_rate_max))
    testit::assert(is.numeric(clado_rate_max))
  }
  ###  Currently don't consider ontogeny and two trait states simutaneously!
  
  # } else if ((Apars$proportional_peak_t * Apars$total_island_age) > timeval) {
  #   ext_rate_max <- max(unlist(ext_rate))
  #   testit::assert(is.numeric(ext_rate_max))
  #   testit::assert(is.null(Tpars))
  #   immig_rate_max <- get_immig_rate(timeval = Apars$proportional_peak_t * Apars$total_island_age,
  #                                    totaltime = totaltime,
  #                                    gam = gam,
  #                                    mainland_n = mainland_n,
  #                                    Apars = Apars,
  #                                    Tpars = Tpars,
  #                                    island_ontogeny = island_ontogeny,
  #                                    island_spec = island_spec,
  #                                    K = K)
  #   testit::assert(is.numeric(immig_rate_max))
  #   
  #   clado_rate_max <- get_clado_rate(timeval = Apars$proportional_peak_t * Apars$total_island_age, # SHOULD BE GENERALIZED
  #                                    lac = lac,
  #                                    Apars = Apars,
  #                                    Tpars = Tpars,
  #                                    island_ontogeny = island_ontogeny,
  #                                    island_spec = island_spec,
  #                                    K = K)
  #   testit::assert(is.numeric(clado_rate_max))
  #   
  # } else {
  #   testit::assert(is.null(Tpars))
  #   # Ontogeny, max rate is t_hor, which in this case is totaltime (from hor)
  #   ext_rate_max <- get_ext_rate(timeval = t_hor,
  #                                mu = mu,
  #                                Apars = Apars,
  #                                Epars = Epars,
  #                                Tpars = Tpars,
  #                                island_ontogeny = island_ontogeny,
  #                                extcutoff = extcutoff,
  #                                island_spec = island_spec,
  #                                K = K)
  #   
  #   testit::assert(is.numeric(ext_rate_max) && ext_rate_max >= 0.0)
  #   immig_rate_max <- immig_rate
  #   testit::assert(is.numeric(immig_rate_max))
  #   clado_rate_max <- clado_rate
  #   testit::assert(is.numeric(clado_rate_max))
  # }
  testit::assert(!is.null(Tpars))
  trans_rate <- get_trans_rate(Tpars = Tpars,
                               island_spec = island_spec)
  Tpars <- create_trait_state_params(trans_rate = trans_rate$trans_rate1,
                                     immig_rate2 = immig_rate$immig_rate2,
                                     ext_rate2 = ext_rate$ext_rate2,
                                     ana_rate2 = ana_rate$ana_rate2,
                                     clado_rate2 = clado_rate$clado_rate2,
                                     trans_rate2 = trans_rate$trans_rate2,
                                     M2 = Tpars$M2)
  immig_rate <- immig_rate$immig_rate1
  ana_rate <- ana_rate$ana_rate1
  ext_rate <- ext_rate$ext_rate1
  clado_rate <- clado_rate$clado_rate1
  rates <- create_rates(
    immig_rate = immig_rate,
    ext_rate = ext_rate,
    ana_rate = ana_rate,
    clado_rate = clado_rate,
    ext_rate_max = ext_rate_max,
    immig_rate_max = immig_rate_max,
    clado_rate_max = clado_rate_max,
    Tpars = Tpars)
  
  return(rates)
}



#' Function to describe changes in area through time. Adapted from
#' Valente et al 2014 ProcB
#' @param timeval current time of simulation
#' @param Apars a named list containing area parameters as created by create_area_params:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#' }
#' @param island_ontogeny a string describing the type of island ontogeny. Can be \code{NULL},
#' \code{"beta"} for a beta function describing area through time,
#'  or \code{"linear"} for a linear function
#' @export
#' @family rates calculation
#' @author Pedro Neves
#' @references Valente, Luis M., Rampal S. Etienne, and Albert B. Phillimore.
#' "The effects of island ontogeny on species diversity and phylogeny."
#' Proceedings of the Royal Society of London B: Biological Sciences 281.1784 (2014): 20133227.
island_area <- function(timeval, Apars, island_ontogeny) {
  testit::assert(are_area_params(Apars))

  Tmax <- Apars$total_island_age
  Amax <- Apars$max_area
  Topt <- Apars$proportional_peak_t
  peak <- Apars$peak_sharpness
  proptime <- timeval/Tmax

  # Constant
  if(island_ontogeny == 0)
  {
    if(Amax != 1 || is.null(Amax))
    {
      warning('Constant ontogeny requires a maximum area of 1.')
    }
    return(1)
  }

  # Linear decline
  if (island_ontogeny == 1) {
    b <- Amax # intercept (peak area)
    m <- -(b / Topt) # slope
    At <- m * timeval + b
    return(At)
  }

  # Beta function
  if (island_ontogeny == 2) {
    f <- Topt / (1 - Topt)
    a <- f * peak / (1 + f)
    b <- peak / (1 + f)
    At <-
      Amax * proptime ^ a * (1 - proptime) ^ b / ((a / (a + b)) ^ a * (b / (a + b)) ^ b)
    return(At)
  }
}


#' Function to describe changes in extinction rate through time. From
#' Valente et al 2014 ProcB
#' @param timeval current time of simulation
#' @param mu per capita extinction rate in no ontogeny model
#' @param Apars a named list containing area parameters as created by create_area_params:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#' }
#' @param Epars a numeric vector:
#' \itemize{
#'   \item{[1]: minimum extinction when area is at peak}
#'   \item{[2]: extinction rate when current area is 0.10 of maximum area}
#' }
#' @param Tpars A named list containing diversification rates considering two trait states:
#' \itemize{
#'   \item{[1]:A numeric with the per capita transition rate with state1}
#'   \item{[2]:A numeric with the per capita immigration rate with state2}
#'   \item{[3]:A numeric with the per capita extinction rate with state2}
#'   \item{[4]:A numeric with the per capita anagenesis rate with state2}
#'   \item{[5]:A numeric with the per capita cladogenesis rate with state2}
#'   \item{[6]:A numeric with the per capita transition rate with state2}
#'   \item{[7]:A numeric with the number of species with trait state 2 on mainland}
#' }
#' @param island_ontogeny a string describing the type of island ontogeny.
#' Can be \code{NULL},
#' \code{"beta"} for a beta function describing area through time,
#'  or \code{"linear"} for a linear function
#' @param extcutoff cutoff for extinction rate preventing it from being too
#' large and slowing down simulation. Default is 1100
#' @param island_spec matrix containing state of system
#' @param K carrying capacity
#' @export
#' @seealso Does the same as \link{DAISIE_calc_clade_ext_rate}
#' @family rates calculation
#' @references Valente, Luis M., Rampal S. Etienne, and Albert B. Phillimore.
#' "The effects of island ontogeny on species diversity and phylogeny."
#' Proceedings of the Royal Society of London B: Biological Sciences 281.1784 (2014): 20133227.
#' @author Pedro Neves
get_ext_rate <- function(timeval,
                         mu,
                         Tpars,
                         Apars,
                         Epars,
                         island_ontogeny,
                         extcutoff = 1100,
                         island_spec,
                         K){
  testit::assert(is.numeric(island_ontogeny))
  if (is.null(Tpars)){    ##without considering trait states
    # Make function accept island_spec matrix or numeric
    if (is.matrix(island_spec) || is.null(island_spec)) {
      N <- length(island_spec[, 1])
    } else if (is.numeric(island_spec)) {
      N <- island_spec
    }
    if (island_ontogeny == 0) {
      ext_rate <- mu * N
      testit::assert(is.numeric(ext_rate))
      testit::assert(ext_rate >= 0)
      return(ext_rate)
    } else {
      X <- log(Epars[1] / Epars[2]) / log(0.1)
      ext_rate <-
        Epars[1] / ((island_area(timeval, Apars, island_ontogeny) /
                       Apars$max_area)^X)
      ext_rate[which(ext_rate > extcutoff)] <- extcutoff
      ext_rate <- ext_rate * N
      testit::assert(is.numeric(ext_rate))
      testit::assert(ext_rate >= 0)
      ext_rate
    }
  }else{   ##species have two states
    if (is.matrix(island_spec) || is.null(island_spec)) {
      N1 <- length(which(island_spec[, 8] == "1"))
      N2 <- length(which(island_spec[, 8] == "2"))
    } else if (is.numeric(island_spec)) {
      stop("Different trait states cannot be separated,please transform to matrix form.")
    }
    ext_rate1 <- mu * N1
    ext_rate2 <- Tpars$ext_rate2 * N2
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
#' @param laa per capita anagenesis rate
#' @param Tpars A named list containing diversification rates considering two trait states:
#' \itemize{
#'   \item{[1]:A numeric with the per capita transition rate with state1}
#'   \item{[2]:A numeric with the per capita immigration rate with state2}
#'   \item{[3]:A numeric with the per capita extinction rate with state2}
#'   \item{[4]:A numeric with the per capita anagenesis rate with state2}
#'   \item{[5]:A numeric with the per capita cladogenesis rate with state2}
#'   \item{[6]:A numeric with the per capita transition rate with state2}
#'   \item{[7]:A numeric with the number of species with trait state 2 on mainland}
#' }
#' @param island_spec matrix with current state of system
#' @seealso Does the same as \link{DAISIE_calc_clade_ana_rate}
#' @family rates calculation
#' @author Pedro Neves
get_ana_rate <- function(laa,
                         Tpars = NULL,
                         island_spec) {
  if(is.null(Tpars)){
    ana_rate = laa * length(which(island_spec[,4] == "I"))
    testit::assert(is.numeric(ana_rate))
    testit::assert(ana_rate >= 0)
    ana_rate
  }else{
    ana_rate1 = laa * length(intersect(which(island_spec[,4] == "I"), which(island_spec[,8] == "1")))
    ana_rate2 = Tpars$ana_rate2 * length(intersect(which(island_spec[,4] == "I"), which(island_spec[,8] == "2")))
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
#' @param timeval current time of simulation
#' @param lac per capita cladogenesis rate
#' @param Apars a named list containing area parameters as created by create_area_params:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#' }
#' @param Tpars A named list containing diversification rates considering two trait states:
#' \itemize{
#'   \item{[1]:A numeric with the per capita transition rate with state1}
#'   \item{[2]:A numeric with the per capita immigration rate with state2}
#'   \item{[3]:A numeric with the per capita extinction rate with state2}
#'   \item{[4]:A numeric with the per capita anagenesis rate with state2}
#'   \item{[5]:A numeric with the per capita cladogenesis rate with state2}
#'   \item{[6]:A numeric with the per capita transition rate with state2}
#'   \item{[7]:A numeric with the number of species with trait state 2 on mainland}
#' }
#' @param island_ontogeny a string describing the type of island ontogeny.
#' Can be \code{NULL},
#' \code{"beta"} for a beta function describing area through time,
#'  or \code{"linear"} for a linear function
#' @param island_spec matrix with current state of system
#' @param K carrying capacity
#' @export
#' @seealso Does the same as \link{DAISIE_calc_clade_clado_rate}
#' @author Pedro Neves
#' @references Valente, Luis M., Rampal S. Etienne, and Albert B. Phillimore.
#' "The effects of island ontogeny on species diversity and phylogeny."
#' Proceedings of the Royal Society of London B: Biological Sciences 281.1784 (2014): 20133227.
get_clado_rate <- function(timeval,
                           lac,
                           Tpars = NULL,
                           Apars,
                           island_ontogeny,
                           island_spec,
                           K) {
  if(is.null(Tpars)) {
    # Make function accept island_spec matrix or numeric
    if (is.matrix(island_spec) || is.null(island_spec)) {
      N <- length(island_spec[,1])
    } else if (is.numeric(island_spec)) {
      N <- island_spec
    }
    # No ontogeny scenario
    testit::assert(is.numeric(island_ontogeny))
    if (island_ontogeny == 0) {
      clado_rate <- max(c(N * lac * (1 - N / K), 0), na.rm = T)
      testit::assert(clado_rate >= 0)
      testit::assert(is.numeric(clado_rate))
      return(clado_rate)
      # Ontogeny scenario
    } else {
      clado_rate <-  max(c(
        N * lac * island_area(timeval, Apars, island_ontogeny) *
          (1 - N / (island_area(
            timeval,
            Apars,
            island_ontogeny) * K)), 0), na.rm = T)
      testit::assert(clado_rate >= 0)
      testit::assert(is.numeric(clado_rate))
      return(clado_rate)
    }
  }else{
    # Shu's work
    if (is.matrix(island_spec) || is.null(island_spec)) {
      N1 <- length(which(island_spec[, 8] == "1"))
      N2 <- length(which(island_spec[, 8] == "2"))
    } else if (is.numeric(island_spec)) {
      stop("Different trait states cannot be separated,please transform to matrix form.")
    }
    clado_rate1 <- max(c(N1 * lac * (1 - (N1 + N2) / K), 0), na.rm = T)
    clado_rate2 <- max(c(N2 * Tpars$clado_rate2 * (1 - (N1 + N2) / K), 0), na.rm = T)
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
#' @param timeval current time of simulation
#' @param totaltime total time of simulation
#' @param gam per capita immigration rate
#' @param Apars a named list containing area parameters as created by create_area_params:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#' }
#' @param Tpars A named list containing diversification rates considering two trait states:
#' \itemize{
#'   \item{[1]:A numeric with the per capita transition rate with state1}
#'   \item{[2]:A numeric with the per capita immigration rate with state2}
#'   \item{[3]:A numeric with the per capita extinction rate with state2}
#'   \item{[4]:A numeric with the per capita anagenesis rate with state2}
#'   \item{[5]:A numeric with the per capita cladogenesis rate with state2}
#'   \item{[6]:A numeric with the per capita transition rate with state2}
#'   \item{[7]:A numeric with the number of species with trait state 2 on mainland}
#' }
#' @param island_ontogeny a string describing the type of island ontogeny.
#' Can be \code{NULL},
#' \code{"beta"} for a beta function describing area through time,
#'  or \code{"linear"} for a linear function
#' @param island_spec matrix with current state of system
#' @param K carrying capacity
#' @param mainland_n total number of species present in the mainland
#' @seealso Does the same as \link{DAISIE_calc_clade_imm_rate}
#' @family rates calculation
#' @author Pedro Neves
#' @references Valente, Luis M., Rampal S. Etienne, and Albert B. Phillimore.
#' "The effects of island ontogeny on species diversity and phylogeny."
#' Proceedings of the Royal Society of London B: Biological Sciences 281.1784 (2014): 20133227.
get_immig_rate <- function(timeval,
                           totaltime,
                           gam,
                           Apars,
                           Tpars = NULL,
                           island_ontogeny,
                           island_spec,
                           K,
                           mainland_n) {
  N <- length(island_spec[, 1])
  testit::assert(is.numeric(island_ontogeny))
  if(is.null(Tpars)){
    if (island_ontogeny == 0) {
      immig_rate <- max(
        c(mainland_n * gam * (1 - length(island_spec[, 1]) / K), 0),
        na.rm = T
      )
      testit::assert(is.numeric(immig_rate))
      testit::assert(immig_rate >= 0)
      return(immig_rate)
    } else {
      immig_rate <- max(c(mainland_n * gam * (1 - length(island_spec[, 1]) / (    ##Trai-DAISIE need to calculate ntotal
        island_area(timeval,
                    Apars,
                    island_ontogeny) * K)), 0), na.rm = T)
    }
    testit::assert(is.numeric(immig_rate))
    testit::assert(immig_rate >= 0)
    immig_rate
  }else{
    if (is.matrix(island_spec) || is.null(island_spec)) {
      N1 <- length(which(island_spec[, 8] == "1"))
      N2 <- length(which(island_spec[, 8] == "2"))
    } else if (is.numeric(island_spec)) {
      stop("Different trait states cannot be separated,please transform to matrix form.")
    }
    mainland_n2 <- Tpars$M2
    gam2 <- Tpars$immig_rate2
    immig_rate1 <- max(
      c(mainland_n * gam * (1 - (N1 + N2) / K), 0),
      na.rm = T
    )
    immig_rate2 <- max(
      c(mainland_n2 * gam2 * (1 - (N1 + N2) / K), 0),
      na.rm = T
    )
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
#' @param Tpars A named list containing diversification rates considering two trait states:
#' \itemize{
#'   \item{[1]:A numeric with the per capita transition rate with state1}
#'   \item{[2]:A numeric with the per capita immigration rate with state2}
#'   \item{[3]:A numeric with the per capita extinction rate with state2}
#'   \item{[4]:A numeric with the per capita anagenesis rate with state2}
#'   \item{[5]:A numeric with the per capita cladogenesis rate with state2}
#'   \item{[6]:A numeric with the per capita transition rate with state2}
#'   \item{[7]:A numeric with the number of species with trait state 2 on mainland}
#' }
#' @param island_spec matrix with current state of system
#' @family rates calculation
get_trans_rate <- function(Tpars,
                           island_spec){
  if(is.null(Tpars)){
    stop("Transition rate only calculate when exists more than one trait state.") #or trans_rate = NULL
  }else{
    # Make function accept island_spec matrix or numeric
    if (is.matrix(island_spec) || is.null(island_spec)) {
      N1 <- length(which(island_spec[, 8] == "1"))
      N2 <- length(which(island_spec[, 8] == "2"))
    } else if (is.numeric(island_spec)) {
      stop("Different trait states cannot be separated,please transform to matrix form.")
    }
    trans_rate1 <- Tpars$trans_rate * N1
    trans_rate2 <- Tpars$trans_rate2 * N2
    testit::assert(is.numeric(trans_rate1))
    testit::assert(trans_rate1 >= 0)
    testit::assert(is.numeric(trans_rate2))
    testit::assert(trans_rate2 >= 0)
    trans_list <- list(trans_rate1 = trans_rate1,
                       trans_rate2 = trans_rate2)
    return(trans_list)
  }
}



#' Function to calculate and update horizon for maximum extinction rate
#' @description Internal function.
#' Calculates when the next horizon for maximum extinction will be in the
#' simulation
#' @param timeval current time of simulation
#' @param totaltime total time of simulation
#' @param Apars a named list containing area parameters as created by create_area_params:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#' }
#' @param Tpars A named list containing diversification rates considering two trait states:
#' \itemize{
#'   \item{[1]:A numeric with the per capita transition rate with state1}
#'   \item{[2]:A numeric with the per capita immigration rate with state2}
#'   \item{[3]:A numeric with the per capita extinction rate with state2}
#'   \item{[4]:A numeric with the per capita anagenesis rate with state2}
#'   \item{[5]:A numeric with the per capita cladogenesis rate with state2}
#'   \item{[6]:A numeric with the per capita transition rate with state2}
#'   \item{[7]:A numeric with the number of species with trait state 2 on mainland}
#' }
#' @param ext_multiplier reduces or increases distance of horizon to current
#' simulation time
#' @param island_ontogeny a string describing the type of island ontogeny.
#'  Can be \code{NULL}, \code{"beta"} for a beta function
#'   describing area through time, or \code{"linear"} for a linear function
#' @param ext effective extinction rate at timeval
#' @param t_hor time of horizon for max extinction
#'
#' @family rates calculation
#' @author Pedro Neves
get_t_hor <- function(timeval,
                      totaltime,
                      Tpars,
                      Apars,
                      ext,
                      ext_multiplier,
                      island_ontogeny,
                      t_hor) {

  ################~~~TODO~~~#####################
  # Use optimize (optimize(island_area, interval = c(0, 10), maximum = TRUE, Apars = create_area_params(1000, 0.1, 1, 17), island_ontogeny = "beta"))
  # to select maximum point to identify maximum of function
  ###############################################
  testit::assert(is.null(Apars) || are_area_params(Apars))
  # Function calculates where the horizon for max(ext_rate) is.
  if (island_ontogeny == 0) {
    testit::assert(totaltime > 0.0)
    t_hor <- totaltime
  } else if(!is.null(Tpars)){
    testit::assert(totaltime > 0.0)
    t_hor <- totaltime
  }else{

    if (is.null(t_hor)) {
      testit::assert(are_area_params(Apars))
      # This is the time at which Amax is reached
      t_hor <- Apars$proportional_peak_t * Apars$total_island_age
    } else if (timeval >= t_hor) {
      # t_hor should dynamically be adjusted depending on parameter values.
      # Certain parameter combinations will always make it be > totaltime at
      # first calculation, slowing down the simulations
      t_hor <- t_hor + t_hor / 6 + ext_multiplier * (totaltime - timeval) * ext
      t_hor <- min(totaltime, t_hor)
    }
    testit::assert(t_hor > 0.0)
  }
  return(t_hor)
}

#' Calculates when the next timestep will be.
#'
#' @param rates list of numeric with probabilities of each event
#' @param timeval current time of simulation
#' @param Tpars A named list containing diversification rates considering two trait states:
#' \itemize{
#'   \item{[1]:A numeric with the per capita transition rate with state1}
#'   \item{[2]:A numeric with the per capita immigration rate with state2}
#'   \item{[3]:A numeric with the per capita extinction rate with state2}
#'   \item{[4]:A numeric with the per capita anagenesis rate with state2}
#'   \item{[5]:A numeric with the per capita cladogenesis rate with state2}
#'   \item{[6]:A numeric with the per capita transition rate with state2}
#'   \item{[7]:A numeric with the number of species with trait state 2 on mainland}
#' }
#' @return named list with numeric vector containing the time of the next
#' timestep and the change in time.
#' @author Pedro Neves
calc_next_timeval <- function(rates = rates,
                              timeval = timeval,
                              Tpars = Tpars) {
  # Calculates when next event will happen
  testit::assert(are_rates(rates))
  testit::assert(timeval >= 0)

  if(is.null(Tpars)){
    totalrate <- rates$immig_rate_max + 
                 rates$ana_rate + 
                 rates$clado_rate_max + 
                 rates$ext_rate_max
  }else{
    totalrate <- rates$immig_rate +
                 rates$ana_rate +
                 rates$clado_rate +
                 rates$ext_rate +
                 rates$trans_rate +
                 rates$immig_rate2 +
                 rates$ana_rate2 +
                 rates$clado_rate2 +
                 rates$ext_rate2 +
                 rates$trans_rate2
  }
    dt <- stats::rexp(1, totalrate)
    timeval <- timeval + dt

  return(list(timeval = timeval, dt = dt))
}


#' Calculate the clade-wide extinction rate
#' @param ps_ext_rate per species extinction rate
#' @param n_species number of species in that clade
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
#' @param ps_ana_rate per species anagensis rate
#' @param n_immigrants number of immigrants in that clade
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

