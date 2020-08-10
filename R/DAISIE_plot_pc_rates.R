#' Plot per capita STT and per-capita rate change thourgh simulation
#'
#' @inheritParams default_params_doc
#'
#' @return A composite plot containing the STT plot, and variation in area,
#' extinction rate, cladogenesis rate and immigration rate
#' @keywords internal
#'
#' @author Pedro Neves
DAISIE_plot_pc_rates <- function(totaltime,
                                 pars,
                                 island_ontogeny,
                                 sea_level = "const",
                                 area_pars,
                                 peak,
                                 hyper_pars,
                                 mainland_n = 1000,
                                 resolution = 0.001,
                                 simulation_outputs = NULL) {
  testit::assert(DAISIE::is_island_ontogeny_input(island_ontogeny))
  testit::assert(DAISIE::are_area_pars(area_pars))
  testit::assert(resolution < 1 && is.numeric(resolution) && !is.infinite(resolution))
  if (!requireNamespace("ggplot2", quietly = TRUE) && !requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Packages \"ggplot2\" and \"gridExtra\" needed for this function to work. Please install them.",
         call. = FALSE)
  }

  lac <- pars[1]
  mu <- pars[2]
  K <- pars[3]
  gam <- pars[4]

  area_plot <- DAISIE_plot_area(
    totaltime = totaltime,
    area_pars = area_pars,
    peak = peak,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    resolution = resolution
  )

  clado_plot <- DAISIE_plot_cladogenesis(
    totaltime = totaltime,
    K = K,
    area_pars = area_pars,
    peak = peak,
    hyper_pars = hyper_pars,
    lac = lac,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    removed_timepoints = 1,
    resolution = resolution
  )

  ext_plot <- DAISIE_plot_extinction(
    totaltime = totaltime,
    area_pars = area_pars,
    peak = peak,
    mu = mu,
    hyper_pars = hyper_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    removed_timepoints = 1,
    resolution = resolution)

  immig_plot <- DAISIE_plot_immigration(
    totaltime = totaltime,
    K = K,
    area_pars = area_pars,
    peak = peak,
    gam = gam,
    mainland_n = 1,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    removed_timepoints = 1,
    resolution = resolution
  )

  if (!is.null(simulation_outputs)) {
    sims_plot <- DAISIE_plot_sims(simulation_outputs)
  }

  gridExtra::grid.arrange(area_plot, clado_plot, immig_plot, ext_plot, nrow = 2)
}
