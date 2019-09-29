#' Plot per capita STT and per-capita rate change thourgh simulation
#'
#' @param simulation_outputs Set of simulation outputs as created by \code{\link{DAISIE_sim}}
#' @param totaltime Numeric; total time to simulate
#' @param pars Contains the model parameters: \cr \cr \code{pars[1]}
#' corresponds to lambda^c (cladogenesis rate) \cr \code{pars[2]} corresponds
#' to mu (extinction rate) \cr \code{pars[3]} corresponds to K (clade-level
#' carrying capacity). Set K=Inf for non-diversity dependence.\cr
#' \code{pars[4]} corresponds to gamma (immigration rate) \cr \code{pars[5]}
#' corresponds to lambda^a (anagenesis rate) \cr \code{pars[6]} corresponds to
#' lambda^c (cladogenesis rate) for type 2 species \cr \code{pars[7]}
#' corresponds to mu (extinction rate) for type 2 species\cr \code{pars[8]}
#' corresponds to K (clade-level carrying capacity) for type 2 species.  Set
#' K=Inf for non-diversity dependence.\cr \code{pars[9]} corresponds to gamma
#' (immigration rate) for type 2 species\cr \code{pars[10]} corresponds to
#' lambda^a (anagenesis rate) for type 2 species\cr The elements 6:10 are
#' optional and are required only when type 2 species are included.
#' @param island_ontogeny a string describing the type of island ontogeny. Can be \code{"const"},
#' \code{"beta"} for a beta function describing area through time,
#'  or \code{"linear"} for a linear function
#' @param Apars A numeric vector:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#' }
#' @param Epars a numeric vector:
#' \itemize{
#'   \item{[1]: minimum extinction when area is at peak}
#'   \item{[2]: extinction rate when current area is 0.10 of maximum area}
#' }
#' @param mainland_n The size of the mainland pool, i.e the number of species that can
#' potentially colonize the island
#' @param resolution numeric indicating resolution of plot. Should be < 0.
#'
#' @return A composite plot containing the STT plot, and variation in area,
#' extinction rate, cladogenesis rate and immigration rate
#' @export
#'
#' @author Pedro Neves
DAISIE_plot_pc_rates <- function(simulation_outputs,
                                 totaltime,
                                 pars,
                                 island_ontogeny,
                                 Apars,
                                 Epars,
                                 mainland_n = 1000,
                                 resolution = 0.001) {
  testit::assert(DAISIE::is_island_ontogeny_input(island_ontogeny))
  testit::assert(DAISIE::are_area_params(Apars))
  testit::assert(resolution < 0 && is.numeric(resolution) && !is.infinite(resolution))
  if (!requireNamespace("ggplot2", quietly = TRUE) && !requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Packages \"ggplot2\" and \"gridExtra\" needed for this function to work. Please install them.",
         call. = FALSE)
  }

  lac <- pars[1]
  K <- pars[3]
  gam <- pars[4]

  area_plot <- DAISIE_plot_area(
    totaltime = totaltime,
    Apars = Apars,
    island_ontogeny = island_ontogeny,
    resolution = resolution
  )

  clado_plot <- DAISIE_plot_cladogenesis(
    totaltime = totaltime,
    K = K,
    Apars = Apars,
    lac = lac,
    island_ontogeny = island_ontogeny,
    removed_timepoints = 1,
    resolution = resolution
  )

  ext_plot <- DAISIE_plot_extinction(
    totaltime = totaltime,
    K = K,
    Apars = Apars,
    Epars = Epars,
    island_ontogeny = island_ontogeny,
    removed_timepoints = 1,
    resolution = resolution)

  immig_plot <- DAISIE_plot_immigration(
    totaltime = totaltime,
    K = K,
    Apars = Apars,
    gam = gam,
    mainland_n = mainland_n,
    island_ontogeny = island_ontogeny,
    removed_timepoints = 1,
    resolution = resolution
  )

  if (!is.null(simulation_outputs)) {
    sims_plot <- DAISIE_plot_sims(simulation_outputs)
  }

  gridExtra::grid.arrange(area_plot, clado_plot, immig_plot, ext_plot, nrow = 2)
}
