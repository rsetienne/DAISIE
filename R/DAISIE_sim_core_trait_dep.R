#' Internal function of the DAISIE simulation
#'
#' @inheritParams default_params_doc
#' @keywords internal
DAISIE_sim_core_trait_dep <- function(
    time,
    mainland_n,
    pars,
    island_ontogeny = 0,
    sea_level = 0,
    hyper_pars,
    area_pars,
    extcutoff = 1000,
    trait_pars = NULL
) {

  #### Initialization ####
  timeval <- 0
  total_time <- time
  island_ontogeny <- translate_island_ontogeny(island_ontogeny)
  sea_level <- translate_sea_level(sea_level)

  if(is.null(trait_pars)){
    stop("A second set of rates should be contain considering two trait states.
         If only one state,run DAISIE_sim_cr instead.")
  }
  testit::assert(length(pars) == 5)

  if (pars[4] == 0 && trait_pars$immig_rate2 == 0) {
    stop("Island has no species and the rate of
    colonisation is zero. Island cannot be colonised.")
  }

  mainland_n2 <- trait_pars$M2
  mainland_ntotal <- mainland_n + mainland_n2
  testit::assert(mainland_ntotal > 0)
  if(mainland_n != 0){
    mainland_spec <- seq(1, mainland_n, 1)
  }else{
    mainland_spec <- c()
  }
  maxspecID <- mainland_ntotal

  island_spec <- c()
  stt_table <- matrix(ncol = 7)
  colnames(stt_table) <- c("Time","nI","nA","nC","nI2","nA2","nC2")
  stt_table[1,] <- c(total_time,0,0,0,0,0,0)
  lac <- pars[1]
  mu <- pars[2]
  K <- pars[3]
  gam <- pars[4]
  laa <- pars[5]

  num_spec <- length(island_spec[, 1])
  num_immigrants <- length(which(island_spec[, 4] == "I"))

  #### Start Monte Carlo iterations ####
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
      extcutoff = extcutoff,
      island_ontogeny = 0,
      sea_level = 0,
      island_spec = island_spec,
      trait_pars = trait_pars
    )
    timeval_and_dt <- calc_next_timeval(
      max_rates = rates,
      timeval = timeval,
      total_time = total_time
    )
    timeval <- timeval_and_dt$timeval

    if (timeval < total_time) {
      possible_event <- DAISIE_sample_event_trait_dep(
        rates = rates
      )

      updated_state <- DAISIE_sim_update_state_trait_dep(
        timeval = timeval,
        total_time = total_time,
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
    total_time = total_time,
    island_spec = island_spec,
    mainland_n = mainland_n,
    trait_pars = trait_pars)
  # ordered_stt_times <- sort(island$stt_table[, 1], decreasing = TRUE)
  # testit::assert(all(ordered_stt_times == island$stt_table[, 1]))
  return(island)
}
