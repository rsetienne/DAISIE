#' Simulates island replicates with an island-wide (IW) diversity-dependent
#' constant-rate process
#'
#' @inheritParams default_params_doc
#' @return A list. The highest level of the least corresponds to each individual
#' replicate. See return for `DAISIE_sim_cr()` for details.
DAISIE_sim_cr_iw <- function(total_time,
                             M,
                             pars,
                             replicates,
                             nonoceanic_pars,
                             sample_freq,
                             hyper_pars,
                             area_pars,
                             cond,
                             verbose) {
  island_replicates <- list()
  for (rep in seq_len(replicates)) {
    if (cond == 0) {
      number_present <- -1
    } else {
      number_present <- 0
    }
    while (number_present < cond) {
      island_replicates[[rep]] <- DAISIE_sim_core_cr(
        time = total_time,
        mainland_n = M,
        pars = pars,
        nonoceanic_pars = nonoceanic_pars,
        hyper_pars = hyper_pars,
        area_pars = area_pars
      )
      temp_island_replicates <- island_replicates[[rep]]
      stac_vec <- unlist(temp_island_replicates)[which(names(unlist(temp_island_replicates)) == "taxon_list.stac")]
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
