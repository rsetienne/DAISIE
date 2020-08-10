#### THESE FUNCTIONS NEED MORE TESTING


#' Plots island area function through time
#'
#' @inheritParams default_params_doc
#'
#' @family rate calculations
#'
#' @return a plot with the area size through time
#' @keywords internal
DAISIE_plot_area <- function(totaltime,
                             area_pars,
                             peak,
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
    area[i] <- island_area(timeval = axis[i],
                                   totaltime = totaltime,
                                   area_pars = area_pars,
                                   peak = peak,
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
#' @keywords internal
DAISIE_plot_extinction <- function(totaltime,
                                   area_pars,
                                   peak,
                                   mu,
                                   hyper_pars,
                                   island_ontogeny = "beta",
                                   sea_level = "const",
                                   removed_timepoints,
                                   resolution,
                                   extcutoff = 1000) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  testit::assert(DAISIE::is_island_ontogeny_input(island_ontogeny))
  island_ontogeny <- translate_island_ontogeny(
    island_ontogeny = island_ontogeny
  )
  sea_level <- translate_sea_level(
    sea_level = sea_level
  )
  axis <- seq(0, totaltime, by = resolution)
  ext_rate <- c()
  A_vector <- sapply(
    X = axis,
    FUN = island_area,
    totaltime = totaltime,
    area_pars = area_pars,
    peak = peak,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level
  )

  ext_rates <- sapply(
    X = A_vector,
    FUN = get_ext_rate,
    extcutoff = extcutoff,
    mu = mu,
    hyper_pars = hyper_pars,
    num_spec = 1
  )

  ext_rate_time <- data.frame(
    Extinction = ext_rates[removed_timepoints:length(ext_rates)],
    Time = axis[removed_timepoints:length(axis)]
  )
  Time <- NULL; rm(Time) # nolint, fixes warning: no visible binding for global variable
  Extinction <- NULL; rm(Extinction) # nolint, fixes warning: no visible binding for global variable
  ext_plot <- ggplot2::ggplot(
    data = ext_rate_time,
    ggplot2::aes(x = Time, y = Extinction)) +
    ggplot2::ggtitle("Variation of per-capita extinction rate")  +
    ggplot2::theme_classic() +
    ggplot2::geom_line(size = 1, color = "red4") +
    ggplot2::ylim(0.2, 2)
  ext_plot
}

#' Plot immigration rate through time
#'
#' @inheritParams default_params_doc
#'
#' @author Pedro Neves
#' @return a plot with per capita immigration rate through time and dataframe with immigration
#' at corresponding time
#' @keywords internal
DAISIE_plot_immigration <- function(totaltime,
                                    K,
                                    area_pars,
                                    gam,
                                    peak,
                                    mainland_n,
                                    hyper_pars = NULL,
                                    island_ontogeny = "beta",
                                    sea_level = "const",
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
  sea_level <- translate_sea_level(
    sea_level = sea_level
  )
  axis <- seq(0, totaltime, by = resolution)
  immig_rate <- c()
  A_vector <- sapply(
    X = axis,
    FUN = island_area,
    totaltime = totaltime,
    area_pars = area_pars,
    peak = peak,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level
  )
  immig_rates <- sapply(
    X = A_vector,
    FUN = get_immig_rate,
    gam = gam,
    num_spec = 1,
    mainland_n = mainland_n,
    K = K
  )

  immig_rate_time <- data.frame(Immigration = immig_rates[removed_timepoints:length(immig_rates)], Time = axis[removed_timepoints:length(axis)])
  Time <- NULL; rm(Time) # nolint, fixes warning: no visible binding for global variable
  Immigration <- NULL; rm(Immigration) # nolint, fixes warning: no visible binding for global variable
  immig_plot <- graphics::plot(ggplot2::ggplot(data = immig_rate_time, ggplot2::aes(x = Time, y = Immigration)) +
                                 ggplot2::ggtitle("Variation of per-capita immigration rate") +
                                 ggplot2::theme_classic() +
                                 ggplot2::geom_line(size = 1.5, color = "blue4"))
  immig_plot
}


#' Plot cladogenesis rate through time
#'
#' @inheritParams default_params_doc
#'
#' @return a plot with per capita cladogenesis rate through time and dataframe with immigration
#' at corresponding time
#' @keywords internal
#'
#' @author Pedro Neves
DAISIE_plot_cladogenesis <- function(totaltime,
                                     K,
                                     area_pars,
                                     peak,
                                     lac,
                                     island_ontogeny = "beta",
                                     sea_level = "const",
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
  sea_level <- translate_sea_level(
    sea_level = sea_level
  )
  axis <- seq(0, totaltime, by = resolution)
  clado_rate <- c()
  A_vector <- sapply(
    X = axis,
    FUN = island_area,
    totaltime = totaltime,
    area_pars = area_pars,
    peak = peak,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level
  )

  clado_rate <- sapply(
    X = A_vector,
    FUN = get_clado_rate,
    lac = lac,
    hyper_pars = hyper_pars,
    num_spec = 1,
    K = K
  )

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
