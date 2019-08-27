#' Internal function of the DAISIE simulation
#'
#' @param time Simulated amount of time
#' @param mainland_n A numeric stating the number of mainland species, that
#'   is, the number of species that can potentially colonize the island.
#'   If \code{\link{DAISIE_sim}} uses a clade-specific diversity dependence,
#'   this value is set to 1. 
#'   If \code{\link{DAISIE_sim}} uses an island-specific diversity dependence,
#'   this value is set to the number of mainland species.
#' @param pars A numeric vector:
#' \itemize{
#'   \item{[1]: cladogenesis rate}
#'   \item{[2]: extinction rate}
#'   \item{[3]: carrying capacity}
#'   \item{[4]: immigration rate}
#'   \item{[5]: anagenesis rate}
#' }
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
#' @param single_trait_state Boolean describing if trait states considered in the model
#' Default is \code{FALSE}
#' @param island_ontogeny A string describing the type of island ontogeny. Can be \code{NULL},
#' \code{beta} for a beta function describing area through time,
#' @param keep_final_state logical indicating if final state of simulation 
#' should be returned. Default is \code{FALSE}
#' @param island_spec A matrix with species on island (state of system at each time point)
DAISIE_sim_core <- function(
  time,
  mainland_n,
  pars,
  ddmodel = c(1,0,1),
  island_type = "oceanic",
  nonoceanic = NULL,
  Apars = NULL,
  Epars = NULL,
  Tpars = NULL,
  single_trait_state = TRUE,
  island_ontogeny = 0,
  keep_final_state = FALSE,
  island_spec = NULL
) {
  testit::assert(is.logical(keep_final_state))
  testit::assert(length(pars) == 5)
  testit::assert(is.null(Apars) || are_area_params(Apars))
  
  # testit::assert(is.null(island_spec) || is.matrix(island_spec))
  
  if (pars[4] == 0) {
    stop('Rate of colonisation is zero. Island cannot be colonised.')
  }  
  
  if (!is.null(Apars) && island_ontogeny == "const") {
    stop("Apars specified for constant island_ontogeny. Set Apars to NULL.")
  }
  
  if ((is.null(Epars) || is.null(Apars)) && (island_ontogeny != 0 && island_ontogeny != "const")) {
    stop(
      "Island ontogeny specified but Area parameters and/or extinction 
      parameters not available. Please either set island_ontogeny to NULL, or 
      specify Apars and Epars."
    )
  }
  if (!is.null(Tpars)) {
    return(
      DAISIE_sim_core_shu(
        time = time,
      mainland_n = mainland_n,
      pars = pars,
      ddmodel = ddmodel,
      island_type = island_type,
      nonoceanic = nonoceanic,
      Apars = Apars,
      Epars = Epars,
      Tpars = Tpars,
      single_trait_state = single_trait_state,
      island_ontogeny = island_ontogeny,
      keep_final_state = keep_final_state,
      island_spec = island_spec
      )
    )
  }
  
  
  timeval <- 0
  totaltime <- time
  lac <- pars[1]
  mu <- pars[2]
  K <- pars[3]
  gam <- pars[4]
  laa <- pars[5]
  extcutoff <- max(1000, 1000 * (laa + lac + gam))
  testit::assert(is.numeric(extcutoff))
  ext_multiplier <- 0.5
  testit::assert((totaltime <= Apars$total_island_age) || is.null(Apars))
  # Make island_ontogeny be numeric
  island_ontogeny <- translate_island_ontogeny(island_ontogeny)
  
  #### Start Gillespie ####
  
  # Start output and tracking objects
  if (is.null(island_spec)) {
    island_spec = c()
    stt_table <- matrix(ncol = 4)
    colnames(stt_table) <- c("Time","nI","nA","nC")
    stt_table[1,] <- c(totaltime,0,0,0)
  } else {
    # stt_table <- matrix(stt_table[nrow(stt_table), ], nrow = 1, ncol = 4)
    stt_table <- matrix(ncol = 4)
    colnames(stt_table) <- c("Time","nI","nA","nC")
    stt_table[1, 1] <- totaltime
    stt_table[1, 2] <- length(which(island_spec[, 4] == "I"))
    stt_table[1, 3] <- length(which(island_spec[, 4] == "A"))
    stt_table[1, 4] <- length(which(island_spec[, 4] == "C"))
  }
  
  mainland_spec <- seq(1, mainland_n, 1)
  maxspecID <- mainland_n
  
  testit::assert(is.null(Apars) || are_area_params(Apars))
  # Pick t_hor (before timeval, to set Amax t_hor)
  t_hor <- get_t_hor(
    timeval = 0,
    totaltime = totaltime,
    Apars = Apars,
    ext = 0,
    ext_multiplier = ext_multiplier,
    island_ontogeny = island_ontogeny, 
    t_hor = NULL
  )
  
  while (timeval < totaltime) {
    # Calculate rates
    rates <- update_rates(
      timeval = timeval,
      totaltime = totaltime,
      gam = gam,
      mu = mu,
      laa = laa,
      lac = lac,
      Apars = Apars,
      Epars = Epars,
      island_ontogeny = island_ontogeny,
      extcutoff = extcutoff,
      K = K,
      island_spec = island_spec,
      mainland_n = mainland_n,
      t_hor = t_hor
    )
    
    testit::assert(timeval >= 0)
    timeval_and_dt <- calc_next_timeval(rates = rates, timeval = timeval)
    timeval <- timeval_and_dt$timeval
    dt <- timeval_and_dt$dt
    
    if (timeval <= t_hor) {
      testit::assert(are_rates(rates))
      
      # Determine event
      possible_event <- DAISIE_sample_event(
        rates = rates,
        island_ontogeny = island_ontogeny
      )
      
      updated_state <- DAISIE_sim_update_state(
        timeval = timeval, 
        totaltime = totaltime,
        possible_event = possible_event,
        maxspecID = maxspecID,
        mainland_spec = mainland_spec,
        island_spec = island_spec,
        stt_table = stt_table
      )
      
      island_spec <- updated_state$island_spec
      maxspecID <- updated_state$maxspecID
      stt_table <- updated_state$stt_table
    } else {
      #### After t_hor is reached ####
      
      timeval <- t_hor
      t_hor <- get_t_hor(
        timeval = timeval,
        totaltime = totaltime,
        Apars = Apars,
        ext = rates$ext_rate,
        ext_multiplier = ext_multiplier,
        island_ontogeny = island_ontogeny, 
        t_hor = t_hor
      )
    }
    # TODO Check if this is redundant, or a good idea
    if (rates$ext_rate_max >= extcutoff && length(island_spec[,1]) == 0) {
      timeval <- totaltime
    }
  }
  
  # Finalize stt_table 
  stt_table <- rbind(
    stt_table, 
    c(
      0, 
      stt_table[nrow(stt_table), 2],
      stt_table[nrow(stt_table), 3],
      stt_table[nrow(stt_table), 4]
    )
  )
  
  island <- DAISIE_create_island(
    stt_table = stt_table,
    totaltime = totaltime,
    island_spec = island_spec,
    mainland_n = mainland_n,
    keep_final_state = keep_final_state
  )
  return(island)  
}  
  
  
  
  
  
  
  
  
DAISIE_sim_core_shu <- function(
  time,
  mainland_n,
  pars,
  ddmodel = c(1,0,1),
  island_type = "oceanic",
  nonoceanic = NULL,
  Apars = NULL,
  Epars = NULL,
  Tpars = NULL,
  single_trait_state = TRUE,
  island_ontogeny = 0,
  keep_final_state = FALSE,
  island_spec = NULL
) {  
  testit::assert(is.logical(keep_final_state))
  testit::assert(length(pars) == 5)
  testit::assert(is.null(Apars) || are_area_params(Apars))
  testit::assert(is.null(Tpars) || is.numeric(Tpars))
  # testit::assert(is.null(island_spec) || is.matrix(island_spec))
  
  if (pars[4] == 0 && island_type == "oceanic") {
    stop('Rate of colonisation is zero. Island cannot be colonised.')
  }  
  
  if (!is.null(Apars) && island_ontogeny == "const") {
    stop("Apars specified for constant island_ontogeny. Set Apars to NULL.")
  }
  
  if ((is.null(Epars) || is.null(Apars)) && (island_ontogeny != 0 && island_ontogeny != "const")) {
    stop(
      "Island ontogeny specified but Area parameters and/or extinction 
      parameters not available. Please either set island_ontogeny to NULL, or 
      specify Apars and Epars."
    )
    if(!is.null(Tpars) && single_trait_state == TRUE){
      stop("Considering single trait state. Set Tpars to NULL")
    }
    if(is.null(Tpars) && single_trait_state == FALSE){
      stop("Considering more than one trait state. Please either set single_trait_state to TRUE, or 
         specify Tpars.")
    }
  }
  
  timeval <- 0
  totaltime <- time
  lac <- pars[1]
  mu <- pars[2]
  K <- pars[3]
  gam <- pars[4]
  laa <- pars[5]
  extcutoff <- max(1000, 1000 * (laa + lac + gam))
  testit::assert(is.numeric(extcutoff))
  ext_multiplier <- 0.5
  testit::assert((totaltime <= Apars$total_island_age) || is.null(Apars))
  # Make island_ontogeny be numeric
  island_ontogeny <- translate_island_ontogeny(island_ontogeny)
  
  #### Start Gillespie ####
  #Considering two trait states
  if(single_trait_state == FALSE){
    testit::assert(!"I never get here")
    # Start output and tracking objects
    if (is.null(island_spec)) {
      island_spec = c()
      stt_table <- matrix(ncol = 7)
      colnames(stt_table) <- c("Time","nI","nA","nC","nI2","nA2","nC2")
      stt_table[1,] <- c(totaltime,0,0,0,0,0,0)
    } else {
      # stt_table <- matrix(stt_table[nrow(stt_table), ], nrow = 1, ncol = 7)
      stt_table <- matrix(ncol = 7)
      colnames(stt_table) <- c("Time","nI","nA","nC","nI2","nA2","nC2")
      stt_table[1, 1] <- totaltime
      stt_table[1, 2] <- length(intersect(which(island_spec[, 4] == "I"),which(island_spec[, 8] == "1")))
      stt_table[1, 3] <- length(intersect(which(island_spec[, 4] == "A"),which(island_spec[, 8] == "1")))
      stt_table[1, 4] <- length(intersect(which(island_spec[, 4] == "C"),which(island_spec[, 8] == "1")))
      stt_table[1, 5] <- length(intersect(which(island_spec[, 4] == "I"),which(island_spec[, 8] == "2")))
      stt_table[1, 6] <- length(intersect(which(island_spec[, 4] == "A"),which(island_spec[, 8] == "2")))
      stt_table[1, 7] <- length(intersect(which(island_spec[, 4] == "C"),which(island_spec[, 8] == "2")))
    }
    mainland_n2 <- Tpars$M2
    mainland_ntotal <- mainland_n + mainland_n2
    mainland_spec1 <- seq(1,mainland_n,1)     
    mainland_spec2 <- seq(mainland_n +1 , mainland_ntotal,1)
    mainland_spec <- seq(1,mainland_ntotal,1)
    maxspecID <- mainland_ntotal
    
    testit::assert(is.null(Apars) || are_area_params(Apars))
    # Pick t_hor (before timeval, to set Amax t_hor)
    t_hor <- get_t_hor(
      timeval = 0,
      totaltime = totaltime,
      Apars = Apars,
      ext = 0,
      ext_multiplier = ext_multiplier,
      island_ontogeny = island_ontogeny, 
      t_hor = NULL
    )
    
    while (timeval < totaltime) {
      # Calculate rates
      rates <- update_rates(
        timeval = timeval,
        totaltime = totaltime,
        gam = gam,
        mu = mu,
        laa = laa,
        lac = lac,
        Tpars = Tpars,
        Apars = Apars,
        Epars = Epars,
        single_trait_state = single_trait_state,
        island_ontogeny = island_ontogeny,
        extcutoff = extcutoff,
        K = K,
        island_spec = island_spec,
        mainland_n = mainland_n,
        t_hor = t_hor
      )
      
      testit::assert(timeval >= 0)
      timeval_and_dt <- calc_next_timeval(rates = rates,
                                          single_trait_state = single_trait_state,
                                          timeval = timeval)
      timeval <- timeval_and_dt$timeval
      dt <- timeval_and_dt$dt
      
      if (timeval <= t_hor) {
        testit::assert(are_rates(rates))
        
        # Determine event
        possible_event <- DAISIE_sample_event(
          rates = rates,
          island_ontogeny = island_ontogeny
        )
        
        updated_state <- DAISIE_sim_update_state_trait(
          timeval = timeval, 
          totaltime = totaltime,
          possible_event = possible_event,
          maxspecID = maxspecID,
          mainland_spec0 = mainland_spec0,
          mainland_spec1 = mainland_spec1,
          island_spec = island_spec,
          stt_table = stt_table
        )
        
        island_spec <- updated_state$island_spec
        maxspecID <- updated_state$maxspecID
        stt_table <- updated_state$stt_table
      } else {
        #### After t_hor is reached ####
        
        timeval <- t_hor
        t_hor <- get_t_hor(
          timeval = timeval,
          totaltime = totaltime,
          Apars = Apars,
          ext = rates$ext_rate,
          ext_multiplier = ext_multiplier,
          island_ontogeny = island_ontogeny, 
          t_hor = t_hor
        )
      }
    }
    # TODO Check if this is redundant, or a good idea
    if (rates$ext_rate_max >= extcutoff && length(island_spec[,1]) == 0) {
      timeval <- totaltime
    }
    
    
    # Finalize stt_table 
    stt_table <- rbind(
      stt_table, 
      c(
        0, 
        stt_table[nrow(stt_table), 2],
        stt_table[nrow(stt_table), 3],
        stt_table[nrow(stt_table), 4],
        stt_table[nrow(stt_table), 5],
        stt_table[nrow(stt_table), 6],
        stt_table[nrow(stt_table), 7]
      )
    )
  }
  island <- DAISIE_create_island(
    stt_table = stt_table,
    totaltime = totaltime,
    island_spec = island_spec,
    mainland_n = mainland_n,
    keep_final_state = keep_final_state
  )
  return(island)
  }
