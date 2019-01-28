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
#' @param island_ontogeny A string describing the type of island ontogeny. 
#' Can be \code{NULL},
#' \code{"quadratic"} for a beta function describing area through time,
#'  or \code{"linear"} for a linear function
#' @param extcutoff A numeric with the cutoff for extinction rate preventing it from being too 
#' large and slowing down simulation. Should be big.
#' @param K A numeric with K (clade-specific carrying capacity)
#' @param island_spec A matrix containing state of system
#' @param mainland_n A numeirc with the total number of species present in the mainland
#' @param t_hor A numeric with the time of horizon for max cladogenesis, immigration and minimum extinction
update_rates <- function(timeval, totaltime,
                         gam, mu, laa, lac, Apars = NULL, Epars = NULL,
                         island_ontogeny = NULL, 
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
  testit::assert(is.null(Apars) || are_area_params(Apars))
  testit::assert(is.null(Epars) || is.numeric(Epars))
  testit::assert(is.character(island_ontogeny) || is.null(island_ontogeny))
  testit::assert(is.numeric(extcutoff) || is.null(extcutoff))
  testit::assert(is.numeric(K))
  testit::assert(is.matrix(island_spec) || is.null(island_spec))
  testit::assert(is.numeric(mainland_n))
  testit::assert(is.numeric(t_hor) || is.null(t_hor))
  
  immig_rate <- get_immig_rate(timeval = timeval,
                               totaltime = totaltime,
                               gam = gam,
                               Apars = Apars,
                               island_ontogeny = island_ontogeny,
                               island_spec = island_spec,
                               K = K,
                               mainland_n = mainland_n)
  testit::assert(is.numeric(immig_rate))
  
  ext_rate <- get_ext_rate(timeval = timeval,
                           mu = mu,
                           Apars = Apars,
                           Epars = Epars, 
                           island_ontogeny = island_ontogeny, 
                           extcutoff = extcutoff,
                           island_spec = island_spec,
                           K = K)
  testit::assert(is.numeric(ext_rate))
  
  ana_rate <- get_ana_rate(laa = laa,
                           island_spec = island_spec)
  testit::assert(is.numeric(ana_rate))
  
  clado_rate <- get_clado_rate(timeval = timeval, 
                               totaltime = totaltime,
                               lac = lac,
                               Apars = Apars,
                               island_ontogeny = island_ontogeny,
                               island_spec = island_spec,
                               K = K)
  testit::assert(is.numeric(clado_rate))
  
  if (is.null(island_ontogeny)) {
    
    immig_rate_max <- immig_rate
    testit::assert(is.numeric(immig_rate_max))
    ext_rate_max <- ext_rate
    testit::assert(is.numeric(ext_rate_max))
    clado_rate_max <- clado_rate
    testit::assert(is.numeric(clado_rate_max))
    
  } else if ((Apars$proportional_peak_t * Apars$total_island_age) > timeval) {
    ext_rate_max <- ext_rate
    testit::assert(is.numeric(ext_rate_max))
    
    immig_rate_max <- get_immig_rate(timeval = Apars$proportional_peak_t * Apars$total_island_age,
                                     totaltime = totaltime,
                                     gam = gam,
                                     mainland_n = mainland_n,
                                     Apars = Apars,
                                     island_ontogeny = island_ontogeny, 
                                     island_spec = island_spec,
                                     K = K)
    testit::assert(is.numeric(immig_rate_max))
    
    clado_rate_max <- get_clado_rate(timeval = Apars$proportional_peak_t * Apars$total_island_age, # SHOULD BE GENERALIZED
                                     totaltime = totaltime,
                                     lac = lac,
                                     Apars = Apars,
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
    clado_rate_max = clado_rate_max
  )
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
#' \code{"quadratic"} for a beta function describing area through time,
#'  or \code{"linear"} for a linear function
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
  if (is.null(island_ontogeny)) {
    return(Apars$max_area)
  }	
  if (island_ontogeny == "quadratic") {

    f <- Topt / (1 - Topt)
    a <- f * peak / (1 + f)
    b <- peak / (1 + f) 
    At <- Amax * proptime ^ a * (1 - proptime) ^ b / ((a / (a + b)) ^ a * (b / (a + b)) ^ b)
    return(At)}
  
  #Linear decline
  if (island_ontogeny == "linear") {
    b <- Amax # intercept (peak area)
    m <- -(b / Topt) # slope
    At <- m * timeval + b
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
#' @param island_ontogeny a string describing the type of island ontogeny. 
#' Can be \code{NULL},
#' \code{"quadratic"} for a beta function describing area through time,
#'  or \code{"linear"} for a linear function
#' @param extcutoff cutoff for extinction rate preventing it from being too 
#' large and slowing down simulation
#' @param island_spec matrix containing state of system
#' @param K carrying capacity
#' @seealso Does the same as \link{DAISIE_calc_clade_ext_rate}
#' @family rates calculation
#' @references Valente, Luis M., Rampal S. Etienne, and Albert B. Phillimore. 
#' "The effects of island ontogeny on species diversity and phylogeny." 
#' Proceedings of the Royal Society of London B: Biological Sciences 281.1784 (2014): 20133227.
#' @author Pedro Neves
get_ext_rate <- function(timeval, 
                         mu,
                         Apars,
                         Epars, 
                         island_ontogeny, 
                         extcutoff,
                         island_spec,
                         K){
  # Epars[1] and Epars[2] (mu_min, mu_p) must be user specified
  if (is.null(island_ontogeny)) {
    extrate <- mu * length(island_spec[,1])
    testit::assert(is.numeric(extrate))
    return(extrate)
    
    } else {
      
    X <- log(Epars[1] / Epars[2]) / log(0.1)
    extrate <- Epars[1] / ((island_area(timeval, Apars, island_ontogeny) / Apars$max_area)^X)
    extrate[which(extrate > extcutoff)] <- extcutoff
    extrate[which(extrate > extcutoff)] <- extcutoff
    extrate <- extrate * length(island_spec[,1])

    testit::assert(is.numeric(extrate))
    extrate
  }
}

#' Calculate anagenesis rate
#' @description Internal function. 
#' Calculates the anagenesis rate given the current number of
#' immigrant species and the per capita rate.
#' @param laa per capita anagenesis rate
#' @param island_spec matrix with current state of system
#' @seealso Does the same as \link{DAISIE_calc_clade_ana_rate}
#' @family rates calculation
#' @author Pedro Neves
get_ana_rate <- function(laa, island_spec) {
  ana_rate = laa * length(which(island_spec[,4] == "I"))
  ana_rate
} 

#' Calculate cladogenesis rate
#' @description Internal function. 
#' Calculates the cladogenesis rate given the current number of
#' species in the system, the carrying capacity and the per capita cladogenesis
#' rate
#' @param timeval current time of simulation
#' @param totaltime total time of simulation
#' @param lac per capita cladogenesis rate
#' @param Apars a named list containing area parameters as created by create_area_params:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the 
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#' }
#' @param island_ontogeny a string describing the type of island ontogeny. 
#' Can be \code{NULL},
#' \code{"quadratic"} for a beta function describing area through time,
#'  or \code{"linear"} for a linear function
#' @param island_spec matrix with current state of system
#' @param K carrying capacity
#'
#' @seealso Does the same as \link{DAISIE_calc_clade_clado_rate}
#' @author Pedro Neves
#' @references Valente, Luis M., Rampal S. Etienne, and Albert B. Phillimore. 
#' "The effects of island ontogeny on species diversity and phylogeny." 
#' Proceedings of the Royal Society of London B: Biological Sciences 281.1784 (2014): 20133227.
get_clado_rate <- function(timeval, 
                           totaltime,
                           lac,
                           Apars,
                           island_ontogeny,
                           island_spec,
                           K) {
  # No ontogeny scenario
  if (is.null(island_ontogeny)) {
    clado_rate <- max(c(length(island_spec[,1])
                        * (lac * (1 - length(island_spec[, 1]) / K)),
                        0),
                      na.rm = T)

    return(clado_rate)
    
    # Ontogeny scenario
  } else {
    
    clado_rate <-  max(c(
      length(island_spec[, 1]) * lac *
        island_area(timeval,
                    Apars,
                    island_ontogeny) *
        (1 - length(island_spec[, 1]) /
           (island_area(timeval,
                        Apars,
                        island_ontogeny) * K)), 0), na.rm = T)
    clado_rate
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
#' @param island_ontogeny a string describing the type of island ontogeny. 
#' Can be \code{NULL},
#' \code{"quadratic"} for a beta function describing area through time,
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
get_immig_rate <- function(
  timeval,
  totaltime,
  gam,
  Apars,
  island_ontogeny,
  island_spec,
  K, 
  mainland_n
) {
  if (is.null(island_ontogeny)) {
    immig_rate <- max(c(mainland_n 
                       * gam * (1 - length(island_spec[, 1]) / K), 0), na.rm = T)
    return(immig_rate)
  } else {
    
    immig_rate <- max(c(mainland_n * gam * (1 - length(island_spec[, 1]) / (
        island_area(timeval,
                    Apars,
                    island_ontogeny) * K)), 0), na.rm = T)
  }
  immig_rate
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
#' @param ext_multiplier reduces or increases distance of horizon to current
#' simulation time
#' @param island_ontogeny a string describing the type of island ontogeny.
#'  Can be \code{NULL}, \code{"quadratic"} for a beta function
#'   describing area through time, or \code{"linear"} for a linear function
#' @param ext effective extinction rate at timeval
#' @param dt change in timeval
#' @param old_timeval timeval before update at current loop
#' @param t_hor time of horizon for max extinction
#'
#' @family rates calculation
#' @author Pedro Neves
get_t_hor <- function(timeval,
                     totaltime,
                     Apars,
                     ext,
                     ext_multiplier,
                     island_ontogeny,
                     t_hor,
                     dt,
                     old_timeval) {

  ################~~~TODO~~~#####################
  # Use optimize (optimize(island_area, interval = c(0, 10), maximum = TRUE, Apars = create_area_params(1000, 0.1, 1, 17), island_ontogeny = "quadratic"))
  # to select maximum point to identify maximum of function
  ###############################################
  testit::assert(is.null(Apars) || are_area_params(Apars))
  # Function calculates where the horizon for max(ext_rate) is.
  if (is.null(island_ontogeny)) {
    testit::assert(totaltime > 0.0)
    return(totaltime)
  } else {
    
    if (is.null(t_hor)) {
      testit::assert(are_area_params(Apars))
      t_hor <- Apars$proportional_peak_t * Apars$total_island_age
      testit::assert(t_hor > 0.0)
      return(t_hor)
      
    } else if (timeval >= t_hor) {
      # t_hor should dynamically be adjusted depending on parameter values.
      # Certain parameter combinations will always make it be > totaltime at 
      # first calculation, slowing down the simulations
      t_hor <- t_hor + t_hor / 6 + ext_multiplier * (totaltime - timeval) * ext
      # t_hor <- old_timeval + ext_multiplier * (totaltime - timeval + dt)
      t_hor <- min(totaltime, t_hor)
      testit::assert(t_hor > 0.0)
      return(t_hor)
    }
  }
}

#' Calculates when the next timestep will be.
#'
#' @param rates list of numeric with probabilities of each event
#' @param timeval current time of simulation
#' @return named list with numeric vector containing the time of the next
#' timestep and the change in time.
#' @author Pedro Neves
calc_next_timeval <- function(rates, timeval) {
  # Calculates when next event will happen
  testit::assert(are_rates(rates))
  testit::assert(timeval >= 0)
  totalrate <- rates$immig_rate_max + rates$ana_rate + rates$clado_rate_max + rates$ext_rate_max
  dt <- rexp(1, totalrate)
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
  max(
    0.0,
    n_species * ps_clado_rate * (1.0 - (n_species / carr_cap))
  )
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
  max(
    0.0,
     n_mainland_species * ps_imm_rate * (1.0 - (n_island_species / carr_cap))
  )
}


