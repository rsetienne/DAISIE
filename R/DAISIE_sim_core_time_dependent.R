#' Internal function of the DAISIE simulation
#'
#' @inheritParams default_params_doc
DAISIE_sim_core_time_dependent <- function(
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
  extcutoff = 1000
) {
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
  if (pars[4] == 0 && nonoceanic_pars[1] == 0) {
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
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    area_pars = area_pars,
    hyper_pars = hyper_pars,
    dist_pars = dist_pars,
    ext_pars = ext_pars,
    totaltime = totaltime,
    pars = pars)
  hyper_pars <- default_metapars$hyper_pars
  dist_pars <- default_metapars$dist_pars
  ext_pars <- default_metapars$ext_pars
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
  maxspecID <- mainland_n
  island_spec <- c()
  stt_table <- matrix(ncol = 4)
  colnames(stt_table) <- c("Time", "nI", "nA", "nC")
  spec_tables <- DAISIE_spec_tables(stt_table,
                                    totaltime,
                                    timeval,
                                    nonoceanic_sample,
                                    island_spec)
  stt_table <- spec_tables$stt_table
  mainland_spec <- spec_tables$mainland_spec
  island_spec <- spec_tables$island_spec

  lac <- pars[1]
  mu <- pars[2]
  K <- pars[3]
  gam <- pars[4]
  laa <- pars[5]

  num_spec <- length(island_spec[, 1])
  num_immigrants <- length(which(island_spec[, 4] == "I"))

  global_max_area_time <- get_global_max_area_time(
    totaltime = totaltime,
    area_pars = area_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level
  )
  global_min_area_time <- get_global_min_area_time(
    totaltime = totaltime,
    area_pars = area_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level
  )
  testit::assert(is.numeric(global_max_area_time))
  testit::assert(is.finite(global_max_area_time))
  testit::assert(is.numeric(global_min_area_time))
  testit::assert(is.finite(global_min_area_time))

  #### Start Monte Carlo ####
  while (timeval < totaltime) {
    max_rates <- update_max_rates(
      timeval = timeval,
      totaltime = totaltime,
      gam = gam,
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
      num_spec = num_spec,
      num_immigrants = num_immigrants,
      mainland_n = mainland_n,
      global_min_area_time = global_min_area_time,
      global_max_area_time = global_max_area_time
    )

    timeval_and_dt <- calc_next_timeval(
      max_rates = max_rates,
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
        mu = numeric(),
        hyper_pars = hyper_pars,
        area_pars = area_pars,
        dist_pars = dist_pars,
        ext_pars = ext_pars,
        island_ontogeny = island_ontogeny,
        sea_level = sea_level,
        extcutoff = extcutoff,
        K = K,
        num_spec = num_spec,
        num_immigrants = num_immigrants,
        mainland_n = mainland_n
      )
      testit::assert(are_rates(rates))
      possible_event <- DAISIE_sample_event_time_dependent(
        max_rates = max_rates
      )
      updated_state <- DAISIE_sim_update_state_time_dependent(
        timeval = timeval,
        totaltime = totaltime,
        possible_event = possible_event,
        maxspecID = maxspecID,
        mainland_spec = mainland_spec,
        island_spec = island_spec,
        stt_table = stt_table,
        rates = rates,
        max_rates = max_rates
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
      stt_table[nrow(stt_table), 4]
    )
  )
  island <- DAISIE_create_island(
    stt_table = stt_table,
    totaltime = totaltime,
    island_spec = island_spec,
    mainland_n = mainland_n)
  return(island)
}
