#' Create parameters for DAISIE_sim
#'
#' @param time Length of the simulation in time units.
#' @param M The size of mainland pool.
#' @param pars model parameters.
#' @param replicates number of simulation replicates.
#' @param mainland_params parameters for simulation mainland processes.
#' @param divdepmodel governing diversity-dependent model
#' @param ddmodel parameters to be diversity-dependent
#' @param island_type oceanic or non-oceanic island
#' @param nonoceanic parameters for non-oceanic island model
#' @param prop_type2_pool Fraction of mainland species that belongs to the
#' second subset of species.
#' @param replicates_apply_type2 Default replicates_apply_type2 = TRUE runs 
#' simulations until the number of islands where a type 2 species has 
#' colonised is equal to the specified number of replicates.
#' @param sample_freq number of units times should be divided by.
#' @param plot_sims Default = TRUE plots species-through-time (STT) plots.
#' @param island_ontogeny a string describing the type of island ontogeny. 
#' Can be \code{"const"}, \code{"beta"}.
#' @param Apars A numeric vector:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: vale from 0 to 1 indicating where in the island's history the
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#' }
#' @param Epars A numeric vector:
#' \itemize{
#'   \item{[1]: minimum extinction when area is at peak}
#'   \item{[2]: extinction rate when current area is 0.10 of maximum area}
#' } 
#' @param keep_final_state logical indicating if final state of simulation 
#' should be returned. Default is \code{FALSE}.  
#' @param stored_data output of DAISIE_sim function when run with 
#' keep_final_state. If not \code{NULL}.
#' @param verbose \code{Default = TRUE} Give intermediate output.
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
#' ddmodel = c(1, 0, 1),
#' island_type = "oceanic",
#' nonoceanic = NULL,
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
                              ddmodel = c(1, 0, 1),
                              island_type = "oceanic",
                              nonoceanic = NULL,
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
  if (island_type == "oceanic" && !is.null(nonoceanic)) {
    warning("Nonoceanic parameters have been specified with an oceanic island.
         Set nonoceanic to NULL")
  }
  if (island_type == "nonoceanic" && is.null(nonoceanic)) {
    stop("Nonoceanic island has no parameters.")
  }
  if (!all(is.numeric(ddmodel)) || length(ddmodel) != 3) {
    stop("ddmodel has to be a numeric vector of length 3")
  }
  params <- list(time = time,
                 M = M,
                 pars = pars,
                 replicates = replicates,
                 mainland_params = mainland_params,
                 divdepmodel = divdepmodel,
                 ddmodel = ddmodel,
                 island_type = island_type,
                 nonoceanic = nonoceanic,
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
