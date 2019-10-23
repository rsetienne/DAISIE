DAISIE_sim_core_new <- function(
  time,
  mainland_n,
  pars,
  ddmodel = 11,
  island_type = "oceanic",
  nonoceanic_params = NULL,
  k_dist_params = NULL,
  island_ontogeny = NULL,
  Apars = NULL,
  Epars = NULL,
  keep_final_state = FALSE,
  island_spec = NULL
) {
  testit::assert(is.logical(keep_final_state))
  testit::assert(is.null(Apars) || are_area_params(Apars))
  if (pars[4] == 0 && island_type == "oceanic") {
    stop("Island has no species and the rate of
	colonisation is zero. Island cannot be colonised.")
  }
  if (!is.null(Apars) && island_ontogeny == "const") {
    stop("Apars specified for constant island_ontogeny. Set Apars to NULL.")
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
  testit::assert((totaltime <= Apars$total_island_age) || is.null(Apars))
  island_ontogeny <- translate_island_ontogeny(island_ontogeny)
  if ((is.null(Epars) || is.null(Apars)) && (island_ontogeny != 0)) {
    stop ("Island ontogeny specified but area parameters and/or extinction
	parameters not available. Please either set island_ontogeny to NULL, or
	specify Apars and Epars.")
  }
  if (island_type == "nonoceanic") {
    nonoceanic_sample <- DAISIE_nonoceanic_spec(prob_samp = nonoceanic_params[1],
                                                prob_nonend = nonoceanic_params[2],
                                                mainland_n = mainland_n)
    init_nonend_spec <- nonoceanic_sample[[1]]
    innit_end_spec <- nonoceanic_sample[[2]]
    mainland_spec <- nonoceanic_sample[[3]]
  }
  if (island_type == "oceanic") {
    mainland_spec <- seq(1,mainland_n,1)
    init_nonend_spec <- 0
    init_end_spec <- 0
  }
  maxspecID <- mainland_n

  if(!is.null(k_dist_params)) {
    K <- rgamma(1, shape = k_dist_params[[1]], rate = k_dist_params[[2]])
  }


  #### Start Gillespie ####
  # Start output  and tracking objects
  if (is.null(island_spec)) {
    island_spec <- c()
    stt_table <- matrix(ncol = 4)
    colnames(stt_table) <- c("Time","nI","nA","nC")
    if (island_type == "oceanic") {
      stt_table[1,] <- c(totaltime,0,0,0)
    } else {
      stt_table[1, ] <- c(totaltime,
                          length(init_nonend_spec),
                          length(init_end_spec),
                          0)
      if (length(init_nonend_spec) == 0) {
        init_nonend_spec <- 0
      }
      if (length(init_end_spec) == 0) {
        init_end_spec <- 0
      }
      if (length(mainland_spec) == 0) {
        mainland_spec <- 0
      }
      if (length(init_nonend_spec) == 1 &&
          init_nonend_spec != 0 || length(init_nonend_spec) > 1) {
        for (i in 1:length(init_nonend_spec)) {
          island_spec <- rbind(island_spec,
                               c(init_nonend_spec[i],
                                 init_nonend_spec[i],
                                 timeval,
                                 "I",
                                 NA,
                                 NA,
                                 NA))
        }
      }
      if (length(init_end_spec) == 1 &&
          init_end_spec != 0 || length(init_end_spec) > 1) {
        for(j in 1:length(init_end_spec)) {
          island_spec <- rbind(island_spec,
                               c(init_end_spec[j],
                                 init_end_spec[j],
                                 timeval,
                                 "A",
                                 NA,
                                 NA,
                                 NA))
        }
      }
    }
  }
  #if starting using keep_final_state
  if (keep_final_state == TRUE) {
    #stt_table <- matrix(stt_table[nrow[stt_table), ], nrow = 1, ncol = 4)
    stt_table <- matrix(ncol = 4)
    colnames(stt_table) <- c("Time", "nI", "nA", "nC")
    stt_table[1, 1] <- totaltime
    stt_table[1, 2] <- length(which(island_spec[, 4] == "I"))
    stt_table[1, 3] <- length(which(island_spec[, 4] == "A"))
    stt_table[1, 4] <- length(which(island_spec[, 4] == "C"))
  }
  testit::assert(is.null(Apars) || are_area_params(Apars))
  # Pick t_hor (before timeval, to set Amax t_hor)
  t_hor <- get_t_hor(
    timeval = 0,
    totaltime = totaltime,
    Apars = Apars,
    ext = 0,
    ext_multiplier = ext_multiplier,
    island_ontogeny = island_ontogeny,
    t_hor = NULL)
  while(timeval < totaltime) {
    island_area <- island_area(timeval, Apars, island_ontogeny)
    if(island_area > land_bridge_threshold) {
      lac <- pars[1]
      mu <- pars[2]
      K <- pars[3]
      gam <- pars[4]
      laa <- pars[5]
    } else {
      lac <- pars[6]
      mu <- pars[7]
      K <- pars[8]
      gam <- pars[9]
      laa <- pars[10]
    }
    # Calculate rates
    rates <- update_rates(
      timeval = timeval,
      totaltime = totaltime,
      gam = gam,
      mu = mu,
      laa = laa,
      lac = lac,
      ddmodel_sim = ddmodel_sim,
      Apars = Apars,
      Epars = Epars,
      island_ontogeny = island_ontogeny,
      extcutoff = extcutoff,
      K = K,
      island_spec = island_spec,
      mainland_n = mainland_n,
      t_hor = t_hor)
    timeval_and_dt <- calc_next_timeval(rates, timeval)
    timeval <- timeval_and_dt$timeval
    dt <- timeval_and_dt$dt ####UP TO HERE
    if ( timeval < pars[11] & ((timeval + dt) >= pars[11])  )
    {
      lac <- pars[6]
      mu <- pars[7]
      K <- pars[8]
      gam <- pars[9]
      laa <- pars[10]
      ext_rate <- mu * length(island_spec[,1])
      ana_rate <- laa * length(which(island_spec[,4] == "I"))
      clado_rate <- max(c(length(island_spec[,1]) * (lac * (1 -length(island_spec[,1])/K)),0),na.rm = T)
      immig_rate <- max( c(mainland_n * gam * (1 - length(island_spec[,1])/K), 0), na.rm = TRUE )
      totalrate <- ext_rate + clado_rate + ana_rate + immig_rate
      dt <- stats::rexp(1, totalrate)
      timeval <- pars[11] + dt
    } else
    {
      timeval <- timeval + dt
    }

    possible_event <- sample(1:4,1,replace = FALSE,c(immig_rate,ext_rate,ana_rate,clado_rate))

    ##############
    if(timeval <= totaltime)
    {
      new_state <- DAISIE_sim_update_state(timeval = timeval,
                                           totaltime = totaltime,
                                           possible_event = possible_event,
                                           maxspecID = maxspecID,
                                           mainland_spec = mainland_spec,
                                           island_spec = island_spec,
                                           stt_table = stt_table)
      island_spec <- new_state$island_spec
      maxspecID <- new_state$maxspecID
    }
    stt_table <- rbind(stt_table,
                       c(totaltime - timeval,
                         length(which(island_spec[,4] == "I")),
                         length(which(island_spec[,4] == "A")),
                         length(which(island_spec[,4] == "C"))
                       )
    )
  }

  stt_table[nrow(stt_table),1] <- 0

  #############
  ### if there are no species on the island branching_times = island_age, stac = 0, missing_species = 0
  if(length(island_spec[,1]) == 0)
  {
    island <- list(stt_table = stt_table, branching_times = totaltime, stac = 0, missing_species = 0)
  } else
  {
    cnames <- c("Species","Mainland Ancestor","Colonisation time (BP)",
                "Species type","branch_code","branching time (BP)","Anagenetic_origin")
    colnames(island_spec) <- cnames

    ### set ages as counting backwards from present
    island_spec[,"branching time (BP)"] <- totaltime - as.numeric(island_spec[,"branching time (BP)"])
    island_spec[,"Colonisation time (BP)"] <- totaltime - as.numeric(island_spec[,"Colonisation time (BP)"])

    if(mainland_n == 1)
    {
      island <- DAISIE_ONEcolonist(totaltime,island_spec,stt_table)
    } else if(mainland_n > 1)
    {
      ### number of colonists present
      colonists_present <- sort(as.numeric(unique(island_spec[,'Mainland Ancestor'])))
      number_colonists_present <- length(colonists_present)

      island_clades_info <- list()
      for(i in 1:number_colonists_present)
      {
        subset_island <- island_spec[which(island_spec[,'Mainland Ancestor']==colonists_present[i]),]
        if(class(subset_island) != 'matrix')
        {
          subset_island <- rbind(subset_island[1:7])
          colnames(subset_island) <- cnames
        }
        island_clades_info[[i]] <- DAISIE_ONEcolonist(totaltime,island_spec=subset_island,stt_table=NULL)
        island_clades_info[[i]]$stt_table <- NULL
      }
      island <- list(stt_table = stt_table, taxon_list = island_clades_info)
    }
  }
  return(island)
}
