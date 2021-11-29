#' Wrapper function around for the DAISIE_format_CS_full_stt and
#' DAISIE_format_CS_sampled_stt
#'
#' @inheritParams default_params_doc
#'
#' @keywords internal
#'
#' @return List with CS DAISIE simulation output
DAISIE_format_CS <- function(island_replicates,
                             time,
                             M,
                             sample_freq = 25,
                             verbose = TRUE,
                             trait_pars = NULL) {
  total_time <- time
  testit::assert(
    !is.na(sample_freq) && !is.null(sample_freq) && sample_freq >= 1
  )
  if (is.infinite(sample_freq)) {
    several_islands <- DAISIE_format_CS_full_stt(
      island_replicates = island_replicates,
      time = total_time,
      M = M,
      verbose = verbose,
      trait_pars = trait_pars
    )
  } else {
    several_islands <- DAISIE_format_CS_sampled_stt(
      island_replicates = island_replicates,
      time = total_time,
      sample_freq = sample_freq,
      M = M,
      verbose = verbose,
      trait_pars = trait_pars
    )
  }
  return(several_islands)
}
