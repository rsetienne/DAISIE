#' Internal function of the DAISIE simulation
#'
#' @inheritParams default_params_doc
DAISIE_sim_core_constant_rate_shift <- function(
  time,
  mainland_n,
  pars,
  nonoceanic_pars = c(0, 0),
  hyper_pars = NULL,
  area_pars = NULL,
  dist_pars = NULL,
  shift_times
) {
  timeval <- 0
  totaltime <- time
  testit::assert(length(pars) == 10 && !is.null(shift_times))
  shift_times <- totaltime - shift_times
  shift_times <- sort(shift_times)
  shift_times <- c(shift_times, Inf)
  dynamic_shift_times <- shift_times

  if (pars[4] == 0 && nonoceanic_pars[1] == 0) {
    stop("Island has no species and the rate of
    colonisation is zero. Island cannot be colonised.")
  }
  default_metapars <- create_default_pars(
    area_pars = area_pars,
    hyper_pars = hyper_pars,
    dist_pars = dist_pars,
    totaltime = totaltime,
    pars = pars
  )
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

  if (shift_times[1] != 0) {
    lac <- pars[1]
    mu <- pars[2]
    K <- pars[3]
    gam <- pars[4]
    laa <- pars[5]
    rate_set <- 1
  } else {
    # Land bridge
    lac <- pars[6]
    mu <- pars[7]
    K <- pars[8]
    gam <- pars[9]
    laa <- pars[10]
    rate_set <- 2
  }
  num_spec <- length(island_spec[, 1])
  num_immigrants <- length(which(island_spec[, 4] == "I"))
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
      K = K,
      num_spec = num_spec,
      num_immigrants = num_immigrants,
      mainland_n = mainland_n,
      island_ontogeny = 0,
      sea_level = 0,
      extcutoff = NULL
    )
    testit::assert(are_rates(rates))

    timeval_and_dt <- calc_next_timeval(
      max_rates = rates,
      timeval = timeval
    )

    timeval <- timeval_and_dt$timeval

    if (timeval >= dynamic_shift_times[1]) {

      # First set of rates for island
      if (rate_set == 2) {
        lac <- pars[1]
        mu <- pars[2]
        K <- pars[3]
        gam <- pars[4]
        laa <- pars[5]
        rate_set <- 1
      } else { # Second set of rates for land bridge
        lac <- pars[6]
        mu <- pars[7]
        K <- pars[8]
        gam <- pars[9]
        laa <- pars[10]
        rate_set <- 2
      }
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
        K = K,
        num_spec = num_spec,
        num_immigrants = num_immigrants,
        mainland_n = mainland_n,
        island_ontogeny = 0,
        sea_level = 0,
        extcutoff = NULL
      )
      timeval_and_dt <- calc_next_timeval(
        max_rates = rates,
        timeval = timeval
      )

      dt <- timeval_and_dt$dt
      timeval <- dynamic_shift_times[1] + dt
      dynamic_shift_times <- dynamic_shift_times[-1]
    }

    possible_event <- DAISIE_sample_event_constant_rate(
      rates = rates
    )

      if (timeval <= totaltime) {

        # Update system

        updated_state <- DAISIE_sim_update_state_constant_rate(
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
