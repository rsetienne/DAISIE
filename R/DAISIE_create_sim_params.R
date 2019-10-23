#' Create parameters for DAISIE_sim
#'
#' @inheritParams default_params_doc
#'
#' @return list of parameters to run the model
#' @export
#'
#' @examples
#' DAISIE_create_sim_params(
#' time = 10,
#' M = 1000,
#' pars = c(2, 2, 40, 0.1, 1),
#' replicates = 100,
#' mainland_params = NULL,
#' divdepmodel = "CS",
#' ddmodel_sim = 11,
#' island_type = "oceanic",
#' nonoceanic_params = NULL,
#' prop_type2_pool = NA,
#' replicates_apply_type2 = TRUE,
#' sample_freq = 25,
#' plot_sims = TRUE,
#' island_ontogeny = "const",
#' Apars = NULL,
#' Epars = NULL,
#' keep_final_state = FALSE,
#' stored_data = NULL,
#' verbose = TRUE)

DAISIE_create_sim_params <- function(time = 10,
                              M = 1000,
                              pars = c(2, 2, 40, 0.1, 1),
                              replicates = 100,
                              mainland_params = NULL,
                              divdepmodel = "CS",
                              ddmodel_sim = 11,
                              island_type = "oceanic",
                              nonoceanic_params = NULL,
                              prop_type2_pool = NA,
                              replicates_apply_type2 = TRUE,
                              sample_freq = 25,
                              plot_sims = TRUE,
                              island_ontogeny = "const",
                              Apars = NULL,
                              Epars = NULL,
                              keep_final_state = FALSE,
                              stored_data = NULL,
                              verbose = TRUE) {
  if (pars[4] == 0 && island_type == "oceanic") {
    stop("Immigration rate is zero with no initial species.")
  }
  if (island_type == "oceanic" && !is.null(nonoceanic_params)) {
    warning("Nonoceanic parameters have been specified with an oceanic
    island. Set nonoceanic_params to NULL")
  }
  if (island_type == "nonoceanic" && is.null(nonoceanic_params)) {
    stop("Nonoceanic island has no parameters.")
  }
  params <- list(time = time,
                 M = M,
                 pars = pars,
                 replicates = replicates,
                 mainland_params = mainland_params,
                 divdepmodel = divdepmodel,
                 ddmodel_sim = ddmodel_sim,
                 island_type = island_type,
                 nonoceanic_params = nonoceanic_params,
                 prop_type2_pool = prop_type2_pool,
                 replicates_apply_type2 = replicates_apply_type2,
                 sample_freq = sample_freq,
                 plot_sims = plot_sims,
                 island_ontogeny = island_ontogeny,
                 Apars = Apars,
                 Epars = Epars,
                 keep_final_state = keep_final_state,
                 stored_data = stored_data,
                 verbose = verbose)
    return(params)
}
