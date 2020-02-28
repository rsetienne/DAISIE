#' Depracated function for single-shift simulation model. Please use
#' DAISIE_sim_constant_rate_shift.
#'
#' @inheritParams default_params_doc
#' @param pars A numeric vector containing the model parameters:
#' \itemize{
#'   \item{pars[1]: lambda^c (cladogenesis rate)}
#'   \item{pars[2]: mu (extinction rate)}
#'   \item{pars[3]: K (carrying capacity), set K=Inf for diversity
#'   independence.}
#'   \item{pars[4]: gamma (immigration rate)}
#'   \item{pars[5]: lambda^a (anagenesis rate)}
#'   \item{pars[6]: lambda^c (cladogenesis rate) for either type 2 species
#'   or rate set 2 in rate shift model}
#'   \item{pars[7]: mu (extinction rate) for either type 2 species or rate
#'   set 2 in rate shift model}
#'   \item{pars[8]: K (carrying capacity) for either type 2 species or rate
#'   set 2 in rate shift model, set K=Inf for diversity independence.}
#'   \item{pars[9]: gamma (immigration rate) for either type 2 species or
#'   rate set 2 in rate shift model}
#'   \item{pars[10]: lambda^a (anagenesis rate) for either type 2 species
#'   or rate set 2 in rate shift model}
#'   \item{pars[11]: time of rate shift before the present}
#' }
#'
#' @return simulation output
#' @export
DAISIE_SR_sim <- function(time,
                          M,
                          pars,
                          replicates,
                          sample_freq = 25,
                          plot_sims = TRUE,
                          verbose = TRUE) {
  shift_times <- pars[11]
  pars <- pars[1:10]
  island_replicates <- DAISIE_sim_constant_rate_shift(time = time,
                                                      M = M,
                                                      pars = pars,
                                                      replicates = replicates,
                                                      shift_times = shift_times,
                                                      sample_freq = sample_freq,
                                                      plot_sims = plot_sims,
                                                      verbose = verbose)
  return(island_replicates)
}
