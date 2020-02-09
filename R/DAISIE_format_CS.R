#' Wrapper function around for the DAISIE_format_CS_full_stt and
#' DAISIE_format_CS_sampled_stt
#'
#' @inheritParams default_params_doc
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
