#' Simulate speciation, immigration and extinction on both
#' mainland and island. If there are no processes,
#' this results in equivalent results as Valente et al., 2015.
#' @inheritParams DAISIE_sim
#' @note This function is still a stub: \code{mainland_params} is
#'   unused and a warning is emitted if this value is non-NULL
#' @author Richel J.C. Bilderbeek
DAISIE_sim_with_mainland <- function(
  time,
  M,
  pars,
  replicates,
  mainland_params = NULL,
  divdepmodel = 'CS',
  prop_type2_pool = NA,
  replicates_apply_type2 = TRUE,
  sample_freq = 25
) {
  if (!is.null(mainland_params)) {
    stop("Mainland speciation not implemented yet")
  }
  return(NULL)
}