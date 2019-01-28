#### THESE FUNCTIONS NEED MORE TESTING


#' Plots island area function through time
#'
#' @param totaltime total time of simulation
#' @param Apars a named list containing area parameters:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the 
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#' }
#' @param island_ontogeny a string describing the type of island ontogeny. Can be \code{NULL},
#' \code{"quadratic"} for a beta function describing area through time,
#'  or \code{"linear"} for a linear function
#' @param resolution numeric indicating resolution of plot. Should be < 0.
#' @family rates calculation
#'
#' @return a plot with the area size through time
#' @export
DAISIE_plot_area <- function(totaltime,
                             Apars,
                             island_ontogeny = "quadratic",
                             resolution) {
  
  testit::assert(are_area_params(Apars))
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  axis <- seq(0, totaltime, by = resolution)
  area <- c()
  for (i in seq_along(axis)) {
    testit::assert(are_area_params(Apars))
    area[i] <- DAISIE::island_area(timeval = axis[i],
                                   Apars = Apars,
                                   island_ontogeny = island_ontogeny
    )
    
  }
  island_area_time <- data.frame(Area = area, Time = axis, Totaltime = totaltime)
  
  Time <- NULL; rm(Time) # nolint, fixes warning: no visible binding for global variable
  Area <- NULL; rm(Area) # nolint, fixes warning: no visible binding for global variable
  plot(ggplot2::ggplot(data = island_area_time, ggplot2::aes(x = Time, y = Area)) +
    ggplot2::geom_line(size = 1.5))
  invisible(island_area_time)
}

#' Plots extinction rate function through time
#'
#' @param totaltime total time of simulation
#' @param K K (clade-level carrying capacity)
#' @param Apars a named list containing area parameters:
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
#' @param island_ontogeny a string describing the type of island ontogeny. Can be \code{NULL},
#' \code{quadratic} for a beta function describing area through time,
#'  or \code{linear} for a linear function
#' @param removed_timepoints starting position of time vector
#' @param resolution resolution of time axis
#'
#' @author Pedro Neves
#' @return per capita extinction rate through time plot and dataframe with extinction 
#' at corresponding time
#' @export
DAISIE_plot_extinction <- function(totaltime,
                                   K, 
                                   Apars, 
                                   Epars, 
                                   island_ontogeny = "quadratic", 
                                   removed_timepoints,
                                   resolution) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  
  axis <- seq(0, totaltime, by = resolution)

  
  
  ext_rate <- c()
  for (i in seq_along(axis)) {
    ext_rate[i] <- DAISIE::get_ext_rate(
      timeval = axis[i],
      Apars = Apars,
      Epars = Epars,
      mu = NA, 
      K = K, 
      extcutoff = 1100, 
      island_spec = matrix(ncol = 1),
      island_ontogeny = island_ontogeny
    )
  }
  
  ext_rate_time <- data.frame(Extinction = ext_rate[removed_timepoints:length(ext_rate)], Time = axis[removed_timepoints:length(axis)])
  
  Time <- NULL; rm(Time) # nolint, fixes warning: no visible binding for global variable
  Extinction <- NULL; rm(Extinction) # nolint, fixes warning: no visible binding for global variable
  plot(ggplot2::ggplot(data = ext_rate_time, ggplot2::aes(x = Time, y = Extinction)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::ylim(0, 5))
  invisible(ext_rate_time)
}

#' Plot immigration rate through time
#'
#' @param totaltime total time of simulation
#' @param K K (clade-level carrying capacity)
#' @param Apars a named list containing area parameters as created by create_area_params:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the 
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#' }
#' @param gam minimum per capita immigration rate
#' @param mainland_n number of mainland species. Set as 1 for clade-specific 
#' diversity dependence
#' @param island_ontogeny a string describing the type of island ontogeny. Can be \code{NULL},
#' \code{quadratic} for a beta function describing area through time,
#'  or \code{linear} for a linear function
#' @param removed_timepoints starting position of time vector
#' @param resolution resolution of time axis
#' 
#' @author Pedro Neves
#' @return a plot with per capita immigration rate through time and dataframe with immigration 
#' at corresponding time
#' @export
DAISIE_plot_immigration <- function(totaltime,
                                    K, 
                                    Apars, 
                                    gam,
                                    mainland_n,
                                    island_ontogeny = "quadratic", 
                                    removed_timepoints,
                                    resolution) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  
  axis <- seq(0, totaltime, by = resolution)

  immig_rate <- c()
  for (i in seq_along(axis)) {
    immig_rate[i] <- get_immig_rate(timeval = axis[i],
                                    totaltime = totaltime,
                                    Apars = Apars,
                                    gam = gam, 
                                    K = K, 
                                    mainland_n = 1, 
                                    island_spec = matrix(ncol = 1),
                                    island_ontogeny = island_ontogeny
    )
  }
  
  immig_rate_time <- data.frame(Immigration = immig_rate[removed_timepoints:length(immig_rate)], Time = axis[removed_timepoints:length(axis)])
  
  Time <- NULL; rm(Time) # nolint, fixes warning: no visible binding for global variable
  Immigration <- NULL; rm(Immigration) # nolint, fixes warning: no visible binding for global variable
  plot(ggplot2::ggplot(data = immig_rate_time, ggplot2::aes(x = Time, y = Immigration)) +
    ggplot2::geom_line(size = 1) )
    # ggplot2::ylim(0, 0.01))
  invisible(immig_rate_time)
}


#' Plot cladogenesis rate through time
#'
#' @param totaltime total time of simulation
#' @param K K (clade-level carrying capacity)
#' @param Apars a named list containing area parameters as created by create_area_params:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the 
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#' }
#' @param lac minimum per capita cladogenesis rate
#' @param island_ontogeny a string describing the type of island ontogeny. Can be \code{NULL},
#' \code{quadratic} for a beta function describing area through time,
#'  or \code{linear} for a linear function
#' @param removed_timepoints starting position of time vector
#' @param resolution resolution of time axis 
#'
#' @return a plot with per capita cladogenesis rate through time and dataframe with immigration 
#' at corresponding time
#' @export
#'
#' @author Pedro Neves
DAISIE_plot_cladogenesis <- function(totaltime,
                                     K, 
                                     Apars, 
                                     lac,
                                     island_ontogeny = "quadratic", 
                                     removed_timepoints,
                                     resolution) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  
  axis <- seq(0, totaltime, by = resolution)
  
  clado_rate <- c()
  for (i in seq_along(axis)) {
    clado_rate[i] <- get_clado_rate(timeval = axis[i],
                                    totaltime = totaltime,
                                    Apars = Apars,
                                    lac = lac, 
                                    K = K,
                                    island_spec = matrix(ncol = 1),
                                    island_ontogeny = island_ontogeny)
  }
  
  clado_rate_time <- data.frame(Cladogenesis = clado_rate[removed_timepoints:length(clado_rate)],
                                Time = axis[removed_timepoints:length(axis)])
  
  Time <- NULL; rm(Time) # nolint, fixes warning: no visible binding for global variable
  Cladogenesis <- NULL; rm(Cladogenesis) # nolint, fixes warning: no visible binding for global variable
  plot(ggplot2::ggplot(data = clado_rate_time, ggplot2::aes(x = Time, y = Cladogenesis)) +
         ggplot2::geom_line(size = 1))
  invisible(clado_rate_time)
}