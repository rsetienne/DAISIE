#' Internal function of the DAISIE simulation
#'
#' @inheritParams default_params_doc
DAISIE_sim_core_trait_dependent <- function(
  time,
  mainland_n,
  pars,
  nonoceanic_pars = c(0, 0),
  island_ontogeny = 0,
  sea_level = 0,
  hyper_pars = NULL,
  area_pars = NULL,
  dist_pars = NULL,
  ext_pars = NULL,
  extcutoff = 1000,
  trait_pars = NULL
) {
  
  #### Initialization ####
  timeval <- 0
  totaltime <- time
  island_ontogeny <- translate_island_ontogeny(island_ontogeny)
  sea_level <- translate_sea_level(sea_level)
  
  testit::assert(length(pars) == 5)
  if (!is.null(area_pars) &&
      (island_ontogeny == 0 && sea_level == 0)) {
    stop("area_pars specified for constant island_ontogeny and sea_level.
         Run DAISIE_sim_constant_rate instead.")
  }
  testit::assert(are_area_pars(area_pars))
  if (pars[4] == 0 && trait_pars$immig_rate2 == 0) {
    stop("Island has no species and the rate of
    colonisation is zero. Island cannot be colonised.")
  }
  if ((is.null(ext_pars) || is.null(area_pars)) &&
      (island_ontogeny != 0 || sea_level != 0)) {
    stop("Island ontogeny and/or sea level specified but area parameters
    and/or extinction parameters not available. Please either set
    island_ontogeny and sea_level to NULL, or specify area_pars and ext_pars.")
  }
  testit::assert(is.numeric(extcutoff))
  default_metapars <- create_default_pars(
    area_pars = area_pars,
    hyper_pars = hyper_pars,
    dist_pars = dist_pars,
    totaltime = totaltime,
    pars = pars)
  hyper_pars <- default_metapars$hyper_pars
  dist_pars <- default_metapars$dist_pars
  area_pars <- default_metapars$area_pars
  
  testit::assert(are_hyper_pars(hyper_pars = hyper_pars))
  testit::assert(are_area_pars(area_pars = area_pars))
  testit::assert(are_dist_pars(dist_pars = dist_pars))
  testit::assert((totaltime <= area_pars$total_island_age) ||
                   is.null(area_pars))
  nonoceanic_sample <- DAISIE_nonoceanic_spec(
    prob_samp = nonoceanic_pars[1],
    prob_nonend = nonoceanic_pars[2],
    mainland_n = mainland_n)
  
  ####  what is the useage of maxspecID and how to set M1 and M2??####
  mainland_n2 <- trait_pars$M2
  mainland_ntotal <- mainland_n + mainland_n2
  testit::assert(mainland_ntotal > 0)
  if(mainland_n != 0){
    mainland_spec <- seq(1, mainland_n, 1)
  }else{
    mainland_spec = c()
  }
  maxspecID <- mainland_ntotal
  
  island_spec <- c()
  stt_table <- matrix(ncol = 7)
  colnames(stt_table) <- c("Time","nI","nA","nC","nI2","nA2","nC2")
  init_nonend_spec <- nonoceanic_sample$init_nonend_spec
  init_end_spec <- nonoceanic_sample$init_end_spec
  stt_table[1,] <- c(totaltime,0,0,0,0,0,0)
  # spec_tables <- list(stt_table = stt_table,
  #                     init_nonend_spec = init_nonend_spec,
  #                     init_end_spec = init_end_spec,
  #                     mainland_spec = mainland_spec,
  #                     island_spec = island_spec)
  # stt_table <- spec_tables$stt_table
  # mainland_spec <- spec_tables$mainland_spec
  # island_spec <- spec_tables$island_spec
  lac <- pars[1]
  mu <- pars[2]
  K <- pars[3]
  gam <- pars[4]
  laa <- pars[5]
  
  num_spec <- length(island_spec[, 1])
  num_spec_trait1 <- length(which(island_spec[,8] == "1"))
  num_spec_trait2 <- length(which(island_spec[,8] == "2"))
  num_immigrants <- length(which(island_spec[, 4] == "I"))
  num_immigrants_trait1 <- length(intersect(which(island_spec[, 4] == "I"),
                                            which(island_spec[, 8] == "1")))
  num_immigrants_trait2 <- length(intersect(which(island_spec[, 4] == "I"),
                                            which(island_spec[, 8] == "2")))
  
  #### Start Monte Carlo iterations ####
  while (timeval < totaltime) {
    rates <- update_rates(
      timeval = timeval,
      totaltime = totaltime,
      gam = gam,
      laa = laa,
      lac = lac,
      mu = mu,
      hyper_pars = hyper_pars,
      area_pars = area_pars,
      dist_pars = dist_pars,
      ext_pars = ext_pars,
      K = K,
      num_spec = num_spec,
      num_immigrants = num_immigrants,
      mainland_n = mainland_n,
      extcutoff = extcutoff,
      island_ontogeny = 0,
      sea_level = 0,
      island_spec = island_spec,
      trait_pars = trait_pars
    )
    testit::assert(are_rates(rates))
    timeval_and_dt <- calc_next_timeval(
      max_rates = rates,
      timeval = timeval
    )
    timeval <- timeval_and_dt$timeval
    
    if (timeval < totaltime) {
      rates <- update_rates(
        timeval = timeval,
        totaltime = totaltime,
        gam = gam,
        laa = laa,
        lac = lac,
        mu = mu,
        hyper_pars = hyper_pars,
        area_pars = area_pars,
        dist_pars = dist_pars,
        ext_pars = NULL,
        K = K,
        num_spec = num_spec,
        num_immigrants = num_immigrants,
        mainland_n = mainland_n,
        extcutoff = NULL,
        island_ontogeny = 0,
        sea_level = 0,
        island_spec = island_spec,
        trait_pars = trait_pars
      )
      testit::assert(are_rates(rates))
      possible_event <- DAISIE_sample_event_trait_dependent(
        rates = rates
      )
      
      updated_state <- DAISIE_sim_update_state_trait_dependent(
        timeval = timeval,
        totaltime = totaltime,
        possible_event = possible_event,
        maxspecID = maxspecID,
        mainland_spec = mainland_spec,
        island_spec = island_spec,
        stt_table = stt_table,
        trait_pars = trait_pars
      )
      
      island_spec <- updated_state$island_spec
      maxspecID <- updated_state$maxspecID
      stt_table <- updated_state$stt_table
      num_spec <- length(island_spec[, 1])
      num_immigrants <- length(which(island_spec[, 4] == "I"))
    }
  }
  #### Finalize STT ####
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
  island <- DAISIE_create_island(
    stt_table = stt_table,
    totaltime = totaltime,
    island_spec = island_spec,
    mainland_n = mainland_n,
    trait_pars = trait_pars)
  return(island)
  }
