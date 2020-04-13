#### THESE FUNCTIONS NEED MORE TESTING


#' Plots island area function through time
#'
#' @inheritParams default_params_doc
#'
#' @family rate calculations
#'
#' @return a plot with the area size through time
#' @export
DAISIE_plot_area <- function(totaltime,
                             area_pars,
                             island_ontogeny = "beta",
                             resolution,
                             sea_level = "const") {
  testit::assert(DAISIE::is_island_ontogeny_input(island_ontogeny))
  island_ontogeny <- translate_island_ontogeny(
    island_ontogeny = island_ontogeny
  )
  sea_level <- translate_sea_level(
    sea_level = sea_level
  )

  testit::assert(are_area_pars(area_pars))
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work.
         Please install it.",
         call. = FALSE)
  }
  axis <- seq(0, totaltime, by = resolution)
  area <- c()
  for (i in seq_along(axis)) {
    testit::assert(are_area_pars(area_pars))
    area[i] <- DAISIE::island_area(timeval = axis[i],
                                   area_pars = area_pars,
                                   island_ontogeny = island_ontogeny,
                                   sea_level = sea_level
    )
  }
  island_area_time <- data.frame(Area = area, Time = axis, Totaltime = totaltime)
  Time <- NULL; rm(Time) # nolint, fixes warning: no visible binding for global variable
  Area <- NULL; rm(Area) # nolint, fixes warning: no visible binding for global variable
  area_plot <- ggplot2::ggplot(
    data = island_area_time,
    ggplot2::aes(x = Time, y = Area)) +
    ggplot2::ggtitle("Variation of island area during simulation")  +
    ggplot2::theme_classic() +
    ggplot2::geom_line(size = 1.5, color = "darkgreen")
  area_plot
}

#' Plots extinction rate function through time
#'
#' @inheritParams default_params_doc
#'
#' @author Pedro Neves
#' @return per capita extinction rate through time plot and dataframe with extinction
#' at corresponding time
#' @export
DAISIE_plot_extinction <- function(totaltime,
                                   K,
                                   area_pars,
                                   mu,
                                   hyper_pars,
                                   island_ontogeny = "beta",
                                   removed_timepoints,
                                   resolution) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  testit::assert(DAISIE::is_island_ontogeny_input(island_ontogeny))
  island_ontogeny <- translate_island_ontogeny(
    island_ontogeny = island_ontogeny
  )
  axis <- seq(0, totaltime, by = resolution)
  ext_rate <- c()
  for (i in seq_along(axis)) {
    ext_rate[i] <- DAISIE::get_ext_rate(
      timeval = axis[i],
      mu = mu,
      hyper_pars = hyper_pars,
      area_pars = area_pars,
      K = K,
      extcutoff = 100,
      num_spec = 1,
      island_ontogeny = island_ontogeny
    )
  }
  ext_rate_time <- data.frame(Extinction = ext_rate[removed_timepoints:length(ext_rate)], Time = axis[removed_timepoints:length(axis)])
  Time <- NULL; rm(Time) # nolint, fixes warning: no visible binding for global variable
  Extinction <- NULL; rm(Extinction) # nolint, fixes warning: no visible binding for global variable
  ext_plot <- ggplot2::ggplot(
    data = ext_rate_time,
    ggplot2::aes(x = Time, y = Extinction)) +
    ggplot2::ggtitle("Variation of per-capita extinction rate")  +
    ggplot2::theme_classic() +
    ggplot2::geom_line(size = 1, color = "red4") + ggplot2::ylim(0.075, min(1, max(ext_rate_time$Extinction)))
  ext_plot
}

#' Plot immigration rate through time
#'
#' @inheritParams default_params_doc
#'
#' @author Pedro Neves
#' @return a plot with per capita immigration rate through time and dataframe with immigration
#' at corresponding time
#' @export
DAISIE_plot_immigration <- function(totaltime,
                                    K,
                                    area_pars,
                                    gam,
                                    mainland_n,
                                    hyper_pars = NULL,
                                    island_ontogeny = "beta",
                                    removed_timepoints,
                                    resolution) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  testit::assert(DAISIE::is_island_ontogeny_input(island_ontogeny))
  island_ontogeny <- translate_island_ontogeny(
    island_ontogeny = island_ontogeny
  )
  axis <- seq(0, totaltime, by = resolution)
  immig_rate <- c()
  for (i in seq_along(axis)) {
    immig_rate[i] <- get_immig_rate(
      timeval = axis[i],
      totaltime = totaltime,
      area_pars = area_pars,
      gam = gam,
      K = K,
      mainland_n = 1,
      num_spec = 1,
      island_ontogeny = island_ontogeny
    )
  }
  immig_rate_time <- data.frame(Immigration = immig_rate[removed_timepoints:length(immig_rate)], Time = axis[removed_timepoints:length(axis)])
  Time <- NULL; rm(Time) # nolint, fixes warning: no visible binding for global variable
  Immigration <- NULL; rm(Immigration) # nolint, fixes warning: no visible binding for global variable
  immig_plot <- graphics::plot(ggplot2::ggplot(data = immig_rate_time, ggplot2::aes(x = Time, y = Immigration)) +
                                 ggplot2::ggtitle("Variation of per-capita immigration rate") +
                                 ggplot2::theme_classic() +
                                 ggplot2::geom_line(size = 1.5, color = "blue4") +
                                 ggplot2::ylim(0, 0.002))
  immig_plot
}


#' Plot cladogenesis rate through time
#'
#' @inheritParams default_params_doc
#'
#' @return a plot with per capita cladogenesis rate through time and dataframe with immigration
#' at corresponding time
#' @export
#'
#' @author Pedro Neves
DAISIE_plot_cladogenesis <- function(totaltime,
                                     K,
                                     area_pars,
                                     lac,
                                     island_ontogeny = "beta",
                                     hyper_pars = NULL,
                                     removed_timepoints,
                                     resolution) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  testit::assert(DAISIE::is_island_ontogeny_input(island_ontogeny))
  island_ontogeny <- translate_island_ontogeny(
    island_ontogeny = island_ontogeny
  )
  axis <- seq(0, totaltime, by = resolution)
  clado_rate <- c()
  for (i in seq_along(axis)) {
    clado_rate[i] <- get_clado_rate(timeval = axis[i],
                                    area_pars = area_pars,
                                    lac = lac,
                                    K = K,
                                    hyper_pars = hyper_pars,
                                    num_spec = 1,
                                    island_ontogeny = island_ontogeny)
  }
  clado_rate_time <- data.frame(Cladogenesis = clado_rate[removed_timepoints:length(clado_rate)],
                                Time = axis[removed_timepoints:length(axis)])
  Time <- NULL; rm(Time) # nolint, fixes warning: no visible binding for global variable
  Cladogenesis <- NULL; rm(Cladogenesis) # nolint, fixes warning: no visible binding for global variable
  clado_plot <- ggplot2::ggplot(data = clado_rate_time,
                                ggplot2::aes(x = Time, y = Cladogenesis)) +
    ggplot2::ggtitle("Variation of per-capita cladogenesis rate") +
    ggplot2::theme_classic() +
    ggplot2::geom_line(size = 1, color = "darkorchid4")
  clado_plot
}
