#' Simulates island replicates with an island-wide (IW) diversity-dependent
#' time-dependent process
#'
#' @inheritParams default_params_doc
#'
#' @return A list. The highest level of the least corresponds to each individual
#' replicate. See return for `DAISIE_sim_time_dep()` for details.
DAISIE_sim_time_dep_iw <- function(total_time,
                                   M,
                                   pars,
                                   replicates,
                                   area_pars,
                                   hyper_pars,
                                   nonoceanic_pars,
                                   sample_freq,
                                   island_ontogeny,
                                   sea_level,
                                   peak,
                                   Amax,
                                   Amin,
                                   extcutoff,
                                   cond,
                                   verbose) {
  island_replicates <- list()
  for (rep in 1:replicates) {
    if (cond == 0) {
      number_present <- -1
    } else {
      number_present <- 0
    }
    while (number_present < cond) {
      island_replicates[[rep]] <- DAISIE_sim_core_time_dep(
        time = total_time,
        mainland_n = M,
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
      stac_vec <- unlist(island_replicates)[which(names(unlist(island_replicates)) == "taxon_list.stac")]
      present <- which(stac_vec != 0)
      number_present <- length(present)
    }
    if (verbose == TRUE) {
      print(paste("Island replicate ", rep, sep = ""))
    }
  }
  island_replicates <- DAISIE_format_IW(
    island_replicates = island_replicates,
    time = total_time,
    M = M,
    sample_freq = sample_freq,
    verbose = verbose)

  return(island_replicates)
}
