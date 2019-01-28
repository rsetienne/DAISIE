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
    warning("Mainland speciation not implemented yet")
  }
  DAISIE_sim(
    time = time,
    M = M,
    pars = pars,
    replicates = replicates,
    mainland_params = NULL,
    divdepmodel = divdepmodel,
    prop_type2_pool = prop_type2_pool,
    replicates_apply_type2 = replicates_apply_type2,
    sample_freq = sample_freq,
    plot_sims = FALSE,
    island_ontogeny = NULL,
    Apars = NULL,
    Epars = NULL,
    verbose = FALSE,
    keep_final_state = FALSE,
    stored_data = NULL
  )
}