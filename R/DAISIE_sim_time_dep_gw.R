#' Simulates island replicates with an guild-wide (GW) diversity-dependent
#' time-dependent process
#'
#' @inheritParams default_params_doc
#'
#' @return A list. The highest level of the least corresponds to each individual
#' replicate. See return for `DAISIE_sim_time_dep()` for details.
DAISIE_sim_time_dep_gw <- function(total_time,
                                   M,
                                   pars,
                                   replicates,
                                   area_pars,
                                   hyper_pars,
                                   nonoceanic_pars,
                                   num_guilds,
                                   sample_freq,
                                   island_ontogeny,
                                   sea_level,
                                   peak,
                                   Amax,
                                   Amin,
                                   extcutoff,
                                   verbose) {
  island_replicates <- list()
  if (!is.numeric(num_guilds)) {
    stop("num_guilds must be numeric")
  }
  guild_size <- M / num_guilds
  testit::assert(num_guilds < M)
  testit::assert(M %% num_guilds == 0)
  for (rep in 1:replicates) {
    island_replicates[[rep]] <- list()
    full_list <- list()
    for (m_spec in 1:num_guilds) {
      full_list[[m_spec]]  <- DAISIE_sim_core_time_dep(
        time = total_time,
        mainland_n = guild_size,
        pars = pars,
        nonoceanic_pars = nonoceanic_pars,
        island_ontogeny = island_ontogeny,
        sea_level = sea_level,
        hyper_pars = hyper_pars,
        area_pars = area_pars,
        peak = peak,
        Amax = Amax,
        Amin = Amin,
        extcutoff = extcutoff
      )
    }
    island_replicates[[rep]] <- full_list
    if (verbose == TRUE) {
      message("Island replicate ", rep)
    }
  }
  island_replicates <- DAISIE_format_GW(
    island_replicates = island_replicates,
    time = total_time,
    M = M,
    sample_freq = sample_freq,
    num_guilds = num_guilds,
    verbose = verbose)

  return(island_replicates)
}
