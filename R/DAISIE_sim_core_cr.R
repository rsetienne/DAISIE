#' Internal function of the DAISIE simulation
#'
#' @inheritParams default_params_doc
#' @keywords internal
DAISIE_sim_core_cr <- function(
  time,
  mainland_n,
  pars,
  nonoceanic_pars,
  hyper_pars,
  area_pars,
  use_rcpp = FALSE
) {

  #### Initialization ####

  total_time <- time

  testit::assert(length(pars) == 5 || length(pars) == 10)

  testit::assert(is.null(area_pars) || are_area_pars(area_pars))
  if (pars[4] == 0 && nonoceanic_pars[1] == 0) {
    stop("Island has no species and the rate of
    colonisation is zero. Island cannot be colonised.")
  }

  nonoceanic_sample <- DAISIE_nonoceanic_spec(
    prob_samp = nonoceanic_pars[1],
    prob_nonend = nonoceanic_pars[2],
    mainland_n = mainland_n)

  maxspecID <- mainland_n
  island_spec <- c()
  stt_table <- matrix(ncol = 4)
  colnames(stt_table) <- c("Time", "nI", "nA", "nC")
  spec_tables <- DAISIE_spec_tables(stt_table,
                                    total_time,
                                    0, # timeval
                                    nonoceanic_sample,
                                    island_spec,
                                    maxspecID)

  #### Start Monte Carlo iterations ####

  if (use_rcpp) {
    results <- .Call("daisie_sim_rcpp", total_time,
                               pars,
                               hyper_pars,
                               area_pars,
                               42,
                               mainland_n,
                               0, # island ontogeny
                               0) # sea level
    stt_table = results$stt_table
    colnames(stt_table) <- c("Time","nI","nA","nC")
    island_spec = results$island_spec
    for (i in 5:7) { # Rcpp writes "NA"
      island_spec[island_spec[,i] == "NA", i] <- NA
    }
  } else {
     results <- timeloop(total_time,
                         pars,
                         hyper_pars,
                         area_pars,
                         mainland_n,
                         spec_tables)
     stt_table <- results$stt_table
     island_spec <- results$island_spec
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
    total_time = total_time,
    island_spec = island_spec,
    mainland_n = mainland_n)
  ordered_stt_times <- sort(island$stt_table[, 1], decreasing = TRUE)
  testit::assert(all(ordered_stt_times == island$stt_table[, 1]))
  return(island)
}

timeloop <- function(total_time,
                     pars,
                     hyper_pars,
                     area_pars,
                     mainland_n,
                     spec_tables) {

  lac <- pars[1]
  mu  <- pars[2]
  K   <- pars[3]
  gam <- pars[4]
  laa <- pars[5]

  maxspecID <- spec_tables$maxspecID
  stt_table <- spec_tables$stt_table
  mainland_spec <- spec_tables$mainland_spec
  island_spec <- spec_tables$island_spec

  num_spec <- length(island_spec[, 1])
  num_immigrants <- length(which(island_spec[, 4] == "I"))

  timeval <- 0

  while (timeval < total_time) {
    rates <- update_rates(
      timeval = timeval,
      total_time = total_time,
      gam = gam,
      laa = laa,
      lac = lac,
      mu = mu,
      hyper_pars = hyper_pars,
      area_pars = area_pars,
      K = K,
      num_spec = num_spec,
      num_immigrants = num_immigrants,
      mainland_n = mainland_n,
      extcutoff = NULL,
      island_ontogeny = 0,
      sea_level = 0
    )
    testit::assert(are_rates(rates))

    timeval_and_dt <- calc_next_timeval(
      max_rates = rates,
      timeval = timeval,
      total_time = total_time
    )

    timeval <- timeval_and_dt$timeval

    if (timeval < total_time) {
      rates <- update_rates(
        timeval = timeval,
        total_time = total_time,
        gam = gam,
        laa = laa,
        lac = lac,
        mu = mu,
        hyper_pars = hyper_pars,
        area_pars = area_pars,
        K = K,
        num_spec = num_spec,
        num_immigrants = num_immigrants,
        mainland_n = mainland_n,
        extcutoff = NULL,
        island_ontogeny = 0,
        sea_level = 0
      )
      testit::assert(are_rates(rates))
      possible_event <- DAISIE_sample_event_cr(
        rates = rates
      )

      updated_state <- DAISIE_sim_update_state_cr(
        timeval = timeval,
        total_time = total_time,
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

  return(list(stt_table = stt_table, island_spec = island_spec))
}
