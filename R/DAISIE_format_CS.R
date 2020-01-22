#' Wrapper function around for the correct DAISIE_format_CS
#'
#' @param island_replicates DAISIE_sim_core simulation output
#' @param time Numeric double with total time of simulation.
#' @param M Int stating number of mainland species.
#' @param sample_freq Int stating how often results are
#' sampled for plotting.
#' @param verbose Logical controling if progress is printed to console.
#'
#' @return List with CS DAISIE simulation output
DAISIE_format_CS <- function(island_replicates,
                             time,
                             M,
                             sample_freq = 25,
                             verbose = TRUE) {
  totaltime <- time
  if (is.infinite(sample_freq)) {
    several_islands <- DAISIE_format_CS_full_stt(
      island_replicates = island_replicates,
      time = totaltime,
      M = M,
      verbose = verbose
    )
  } else {
    several_islands <- DAISIE_format_CS_sampled_stt(
      island_replicates = island_replicates,
      sample_freq = sample_freq,
      time = totaltime,
      M = M,
      verbose = verbose
    )
  }
  return(several_islands)
}
