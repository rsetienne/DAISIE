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
#' @param totaltime total time of simulation
#' @param K K (clade-level carrying capacity)
#' @param area_pars a named list containing area and sea level parameters as
#' created by \code{\link{create_area_pars}}:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#'   \item{[5]: amplitude of area fluctuation from sea level}
#'   \item{[6]: frequency of sine wave of area change from sea level}
#' }
#' @param ext_pars a numeric vector:
#' \itemize{
#'   \item{[1]: minimum extinction when area is at peak}
#'   \item{[2]: extinction rate when current area is 0.10 of maximum area}
#' }
#' @param island_ontogeny a string describing the type of island ontogeny. Can be \code{NULL},
#' \code{beta} for a beta function describing area through time.
#' @param removed_timepoints starting position of time vector
#' @param resolution resolution of time axis
#'
#' @author Pedro Neves
#' @return per capita extinction rate through time plot and dataframe with extinction
#' at corresponding time
#' @export
DAISIE_plot_extinction <- function(totaltime,
                                   K,
                                   area_pars,
                                   ext_pars,
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
      area_pars = area_pars,
      ext_pars = ext_pars,
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
#' @param totaltime total time of simulation
#' @param K K (clade-level carrying capacity)
#' @param area_pars a named list containing area and sea level parameters as
#' created by \code{\link{create_area_pars}}:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#'   \item{[5]: amplitude of area fluctuation from sea level}
#'   \item{[6]: frequency of sine wave of area change from sea level}
#' }
#' @param gam minimum per capita immigration rate
#' @param mainland_n number of mainland species. Set as 1 for clade-specific
#' diversity dependence
#' @param island_ontogeny a string describing the type of island ontogeny. Can be \code{NULL},
#' \code{beta} for a beta function describing area through time.
#' @param removed_timepoints starting position of time vector
#' @param resolution resolution of time axis
#' @param hyper_pars A numeric vector for hyperparameters for the rate
#' calculations:
#' \itemize{
#' \item{[1]: is d_0 the scaling parameter for exponent for calculating
#' cladogenesis rate}
#' \item{[2]: is x the exponent for calculating extinction rate}
#' \item{[3]: is alpha, the exponent for calculating the immigration rate}
#' \item{[4]: is beta the exponent for calculating the anagenesis rate.}
#' }
#' @param dist_pars a numeric for the distance from the mainland.
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
                                    dist_pars = NULL,
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
      hyper_pars = hyper_pars,
      dist_pars = dist_pars,
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
#' @param totaltime total time of simulation
#' @param K K (clade-level carrying capacity)
#' @param area_pars a named list containing area and sea level parameters as
#' created by \code{\link{create_area_pars}}:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#'   \item{[5]: amplitude of area fluctuation from sea level}
#'   \item{[6]: frequency of sine wave of area change from sea level}
#' }
#' @param lac minimum per capita cladogenesis rate
#' @param island_ontogeny a string describing the type of island ontogeny. Can be \code{NULL},
#' \code{beta} for a beta function describing area through time.
#' @param removed_timepoints starting position of time vector
#' @param resolution resolution of time axis
#' @param hyper_pars A numeric vector for hyperparameters for the rate
#' calculations:
#' \itemize{
#' \item{[1]: is d_0 the scaling parameter for exponent for calculating
#' cladogenesis rate}
#' \item{[2]: is x the exponent for calculating extinction rate}
#' \item{[3]: is alpha, the exponent for calculating the immigration rate}
#' \item{[4]: is beta the exponent for calculating the anagenesis rate.}
#' }
#' @param dist_pars a numeric for the distance from the mainland.
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
                                     dist_pars = NULL,
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
                                    dist_pars = dist_pars,
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
