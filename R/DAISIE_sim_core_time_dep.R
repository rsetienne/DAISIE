#' Internal function of the DAISIE simulation
#'
#' @inheritParams default_params_doc
#' @param area_windows Optional list as returned by `get_area_windows()`.
#'   If `NULL` it is computed. Giving it avoids recomputation when the core
#'   is called repeatedly.
#' @keywords internal
DAISIE_sim_core_time_dep <- function(
  time,
  mainland_n,
  pars,
  nonoceanic_pars,
  island_ontogeny = 0,
  sea_level = 0,
  hyper_pars,
  area_pars,
  peak,
  Amax,
  Amin,
  extcutoff = 1000,
  area_windows = NULL,
  use_rcpp = FALSE
) {
  timeval <- 0
  total_time <- time
  island_ontogeny <- translate_island_ontogeny(island_ontogeny)
  sea_level <- translate_sea_level(sea_level)

  testit::assert(length(pars) == 5)
  if (!is.null(area_pars) &&
      (island_ontogeny == 0 && sea_level == 0)) {
    stop("area_pars specified for constant island_ontogeny and sea_level.
         Run DAISIE_sim_cr instead.")
  }
  testit::assert(are_area_pars(area_pars))
  if (pars[4] == 0 && nonoceanic_pars[1] == 0) {
    stop("Island has no species and the rate of
    colonisation is zero. Island cannot be colonised.")
  }
  if (is.null(area_pars) && (island_ontogeny != 0 || sea_level != 0)) {
    stop("Island ontogeny and/or sea level specified but area parameters not
    available. Please either set island_ontogeny and sea_level to NULL, or
    specify area_pars or sea_level.")
  }
  testit::assert(is.numeric(extcutoff))

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
                                    timeval,
                                    nonoceanic_sample,
                                    island_spec,
                                    maxspecID)
  maxspecID <- spec_tables$maxspecID
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

  #### Start Monte Carlo ####
  # Thinning with piecewise-constant bounds: within each time window the
  # candidate events are drawn with the maximum rates for the minimum and
  # maximum area of that window (see get_area_windows()).
  if (is.null(area_windows)) {
    area_windows <- get_area_windows(
      total_time = total_time,
      area_pars = area_pars,
      peak = peak,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      Amin = Amin,
      Amax = Amax
    )
  }
  if (use_rcpp) {
    results <- run_sim_engine(list(
      mode = "time_dep",
      total_time = total_time,
      pars = pars,
      d = hyper_pars$d,
      x = hyper_pars$x,
      mainland_n = mainland_n,
      mainland_spec = as.integer(mainland_spec),
      area_pars = area_pars,
      island_ontogeny = as.integer(island_ontogeny),
      sea_level = as.integer(sea_level),
      peak = if (is.null(peak)) 0 else peak,
      Amax = Amax,
      Amin = Amin,
      extcutoff = extcutoff,
      area_windows = area_windows,
      island_spec = daisie_engine_island_spec(island_spec),
      maxspecID = maxspecID,
      stt_start = stt_table[1, ]
    ))
    island_spec <- results$island_spec
    stt_table <- results$stt_table
  } else {
    n_win <- length(area_windows$end)
    win <- 1L
    cache_key <- c(-1, -1, -1)
    max_rates <- NULL
    while (timeval < total_time) {
      while (win < n_win && timeval >= area_windows$end[win]) {
        win <- win + 1L
      }
      if (num_spec == 0) {
        # Without species only immigration is possible, at a rate that does not
        # depend on the area, so no area window is needed.
        win_end <- total_time
        win_Amin <- Amin
        win_Amax <- Amax
        key <- c(0, 0, 0)
      } else {
        win_end <- area_windows$end[win]
        win_Amin <- area_windows$Amin[win]
        win_Amax <- area_windows$Amax[win]
        key <- c(num_spec, num_immigrants, win)
      }
      if (!all(key == cache_key)) {
        max_rates <- update_max_rates(
          gam = gam,
          laa = laa,
          lac = lac,
          mu = mu,
          hyper_pars = hyper_pars,
          extcutoff = extcutoff,
          K = K,
          num_spec = num_spec,
          num_immigrants = num_immigrants,
          mainland_n = mainland_n,
          Amin = win_Amin,
          Amax = win_Amax
        )
        cache_key <- key
      }

      timeval_and_dt <- calc_next_timeval(
        max_rates = max_rates,
        timeval = timeval,
        total_time = total_time
      )
      if (timeval_and_dt$timeval >= win_end) {
        # no candidate event in this window: move to the start of the next
        timeval <- win_end
        next
      }
      timeval <- timeval_and_dt$timeval

      rates <- update_rates(
        timeval = timeval,
        total_time = total_time,
        gam = gam,
        laa = laa,
        lac = lac,
        mu = mu,
        hyper_pars = hyper_pars,
        area_pars = area_pars,
        peak = peak,
        island_ontogeny = island_ontogeny,
        sea_level = sea_level,
        extcutoff = extcutoff,
        K = K,
        num_spec = num_spec,
        num_immigrants = num_immigrants,
        mainland_n = mainland_n
      )
      if (rates$ext_rate > max_rates$ext_max_rate * (1 + 1e-8) ||
          rates$immig_rate > max_rates$immig_max_rate * (1 + 1e-8) ||
          rates$clado_rate > max_rates$clado_max_rate * (1 + 1e-8)) {
        stop("Rate exceeds its upper bound in the time-dependent simulation.")
      }
      possible_event <- DAISIE_sample_event_time_dep(
        max_rates = max_rates
      )
      updated_state <- DAISIE_sim_update_state_time_dep(
        timeval = timeval,
        total_time = total_time,
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
    total_time = total_time,
    island_spec = island_spec,
    mainland_n = mainland_n)
  # ordered_stt_times <- sort(island$stt_table[, 1], decreasing = TRUE)
  # testit::assert(all(ordered_stt_times == island$stt_table[, 1]))
  return(island)
}
