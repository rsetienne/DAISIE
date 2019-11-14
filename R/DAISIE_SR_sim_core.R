DAISIE_SR_sim_core <- function(time,
                               mainland_n,
                               pars,
                               ddmodel_sim = 11,
                               island_type = "oceanic",
                               nonoceanic_pars = NULL,
                               k_dist_pars = NULL,
                               island_ontogeny = 0,
                               sea_level = 0,
                               area_pars = NULL,
                               ext_pars = NULL,
                               shift_times = NULL) {
  timeval <- 0
  totaltime <- time
  lac <- pars[1]
  mu <- pars[2]
  K <- pars[3]
  gam <- pars[4]
  laa <- pars[5]

  if(pars[4] == 0) {
    stop('Rate of colonisation is zero. Island cannot be colonised.')
  }
  if (island_type == "nonoceanic") {
    nonoceanic_sample <- DAISIE_nonoceanic_spec(prob_samp = nonoceanic_pars[1],
                                                prob_nonend = nonoceanic_pars[2],
                                                mainland_n = mainland_n)
    init_nonend_spec_vec <- nonoceanic_sample[[1]]
    init_end_spec_vec <- nonoceanic_sample[[2]]
    mainland_spec <- nonoceanic_sample[[3]]
  }
  if (island_type == "oceanic") {
    mainland_spec <- seq(1, mainland_n, 1)
    init_nonend_spec <- 0
    init_end_spec <- 0
  }
  maxspecID <- mainland_n
  island_spec <- c()
  stt_table <- matrix(ncol = 4)
  colnames(stt_table) <- c("Time", "nI", "nA", "nC")
  if (island_type == "oceanic") {
    stt_table[1, ] <- c(totaltime, 0, 0, 0)
  } else {
    nonoceanic_tables <- DAISIE_nonoceanic_stt_table(stt_table,
                                                     totaltime,
                                                     timeval,
                                                     init_nonend_spec_vec,
                                                     init_end_spec_vec,
                                                     mainland_spec,
                                                     island_spec
    )
    stt_table <- nonoceanic_tables$stt_table
    init_nonend_spec <- nonoceanic_tables$init_nonend_spec
    init_end_spec <- nonoceanic_tables$init_end_spec
    mainland_spec <- nonoceanic_tables$mainland_spec
    island_spec <- nonoceanic_tables$island_spec
  }

  t_hor <- get_t_hor(
    timeval = 0,
    totaltime = totaltime,
    area_pars = area_pars,
    ext = 0,
    ext_multiplier = ext_multiplier,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    t_hor = NULL)
  while (timeval < totaltime) {
    land_bridge_period <- land_bridge_periods(timeval, shift_times)
    if (land_bridge_period[[1]] == FALSE) {
      lac <- pars[1]
      mu <- pars[2]
      K <- pars[3]
      gam <- pars[4]
      laa <- pars[5]
    }
    if (land_bridge_period[[1]] == TRUE) {
      lac <- pars[6]
      mu <- pars[7]
      K <- pars[8]
      gam <- pars[9]
      laa <- pars[10]
    }
    ext_rate <- mu * length(island_spec[,1])
    ana_rate <- laa * length(which(island_spec[,4] == "I"))
    clado_rate <- max(c(length(island_spec[,1]) * (lac * (1 -length(island_spec[,1])/K)),0),na.rm = T)
    immig_rate <- max(c(mainland_n * gam * (1 - length(island_spec[,1])/K),0),na.rm = T)
    totalrate <- ext_rate + clado_rate + ana_rate + immig_rate
    dt <- stats::rexp(1,totalrate)
    next_time_step <- timeval + dt
    land_bridge_period_plus_dt <- land_bridge_periods(next_time_step, shift_times)
    if (land_bridge_period[[1]] == FALSE & land_bridge_periods_plus_dt[[1]] == TRUE) {
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
      timeval <- shift_times[land_bridge_period[[2]]] + dt
    }
    if (land_bridge_period[[1]] == TRUE & land_bridge_period_plus_dt[[1]] == FALSE) {
      lac <- pars[1]
      mu <- pars[2]
      K <- pars[3]
      gam <- pars[4]
      laa <- pars[5]
      ext_rate <- mu * length(island_spec[,1])
      ana_rate <- laa * length(which(island_spec[,4] == "I"))
      clado_rate <- max(c(length(island_spec[,1]) * (lac * (1 -length(island_spec[,1])/K)),0),na.rm = T)
      immig_rate <- max(c(mainland_n * gam * (1 - length(island_spec[,1])/K),0),na.rm = T)
      totalrate <- ext_rate + clado_rate + ana_rate + immig_rate
      dt <- stats::rexp(1,totalrate)
      timeval <- shift_times[land_bridge_period[[2]]] + dt
    } else {
      timeval <- timeval + dt
    }

    possible_event <- sample(1:4,1,replace = FALSE,c(immig_rate,ext_rate,ana_rate,clado_rate))

    ##############
    if (timeval <= totaltime) {
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
    island <- list(stt_table = stt_table,
                   branching_times = totaltime,
                   stac = 0,
                   missing_species = 0,
                   init_nonend_spec = 0, #needs changing only here
                   init_end_spec = 0,    #to make the code pass the
                   carrying_capacity = 0) #test and work with
                                             #format_CS)
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
      island <- list(stt_table = stt_table,
                     taxon_list = island_clades_info,
                     init_nonend_spec = 0, #needs changing only here
                     init_end_spec = 0,    #to make the code pass the
                     carrying_capacity = 0) #test and work with
    }                                           #format_CS
  }
  return(island)
}
