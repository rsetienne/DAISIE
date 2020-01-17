#' Create parameters for DAISIE_sim
#'
#' @inheritParams default_params_doc
#'
#' @return list of parameters to run the model
#' @export
#'
#' @examples
#' DAISIE_create_sim_pars(
#' time = 10,
#' M = 1000,
#' pars = c(2, 2, 40, 0.1, 1),
#' replicates = 100,
#' divdepmodel = "CS",
#' island_type = "oceanic",
#' nonoceanic_pars = NULL,
#' prop_type2_pool = NA,
#' replicates_apply_type2 = TRUE,
#' sample_freq = 25,
#' plot_sims = TRUE,
#' island_ontogeny = "const",
#' area_pars = NULL,
#' ext_pars = NULL,
#' verbose = TRUE)

DAISIE_create_sim_pars <- function(time = 10,
                              M = 1000,
                              pars = c(2, 2, 40, 0.1, 1),
                              replicates = 100,
                              divdepmodel = "CS",
                              island_type = "oceanic",
                              nonoceanic_pars = NULL,
                              prop_type2_pool = NA,
                              replicates_apply_type2 = TRUE,
                              sample_freq = 25,
                              plot_sims = TRUE,
                              island_ontogeny = "const",
                              area_pars = NULL,
                              ext_pars = NULL,
                              verbose = TRUE) {
  if (pars[4] == 0 && island_type == "oceanic") {
    stop("Immigration rate is zero with no initial species.")
  }
  if (island_type == "oceanic" && !is.null(nonoceanic_pars)) {
    warning("Nonoceanic parameters have been specified with an oceanic
    island. Set nonoceanic_pars to NULL")
  }
  if (island_type == "nonoceanic" && is.null(nonoceanic_pars)) {
    stop("Nonoceanic island has no parameters.")
  }
  pars <- list(time = time,
                 M = M,
                 pars = pars,
                 replicates = replicates,
                 divdepmodel = divdepmodel,
                 island_type = island_type,
                 nonoceanic_pars = nonoceanic_pars,
                 prop_type2_pool = prop_type2_pool,
                 replicates_apply_type2 = replicates_apply_type2,
                 sample_freq = sample_freq,
                 plot_sims = plot_sims,
                 island_ontogeny = island_ontogeny,
                 area_pars = area_pars,
                 ext_pars = ext_pars,
                 verbose = verbose)
    return(pars)
}
