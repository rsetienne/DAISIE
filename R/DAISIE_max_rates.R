#' Calculates the maximum rates for a Gillespie simulation
#' @description Internal function that updates the all the max rates at time t.
#' @family rate calculations
#'
#' @inheritParams default_params_doc
#'
#' @seealso \code{\link{update_rates}()}
#'
#' @return a named list with the updated effective rates.
#' @keywords internal
update_max_rates <- function(gam,
                             laa,
                             lac,
                             mu,
                             hyper_pars = NULL,
                             extcutoff,
                             K,
                             num_spec,
                             num_immigrants,
                             mainland_n,
                             Amin,
                             Amax) {

  immig_max_rate <- get_immig_rate(
    gam = gam,
    A = Amax,
    num_spec = num_spec,
    K = K,
    mainland_n = mainland_n
  )

  # testit::assert(is.numeric(immig_max_rate))
  clado_max_rate <- get_clado_rate(
    lac = lac,
    hyper_pars = hyper_pars,
    num_spec = num_spec,
    K = K,
    A = Amax
  )
  # testit::assert(is.numeric(clado_max_rate))

  ext_max_rate <- get_ext_rate(
    mu = mu,
    hyper_pars = hyper_pars,
    extcutoff = extcutoff,
    num_spec = num_spec,
    A = Amin
  )
  # testit::assert(is.numeric(ext_max_rate) && ext_max_rate >= 0.0)

  ana_max_rate <- get_ana_rate(
    laa = laa,
    num_immigrants = num_immigrants
  )
  # testit::assert(is.numeric(ana_max_rate) && ana_max_rate >= 0.0)

  max_rates <- list(
    ext_max_rate = ext_max_rate,
    immig_max_rate = immig_max_rate,
    ana_max_rate = ana_max_rate,
    clado_max_rate = clado_max_rate
  )
  return(max_rates)
}


#' Get the maximum area
#'
#' @inheritParams default_params_doc
#'
#' @return Numeric maximum area during the simulation.
#'
#' @noRd
#' @author Pedro Neves, Joshua Lambert, Shu Xie
get_global_max_area <- function(total_time,
                                area_pars,
                                peak,
                                island_ontogeny,
                                sea_level) {

  max <- stats::optimize(
    f = island_area,
    interval = c(0, total_time),
    total_time = total_time,
    area_pars = area_pars,
    peak = peak,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level, # Fixed at no sea_level for the moment
    maximum = TRUE,
    tol = .Machine$double.eps
  )

  global_max_area_time <- max$maximum

  # testit::assert(is.numeric((global_max_area_time)))
  global_max_area_time <- DDD::roundn(global_max_area_time, 14)

  Amax <- island_area(
    timeval = global_max_area_time,
    total_time = total_time,
    area_pars = area_pars,
    peak = peak,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level
  )
  return(Amax)
}

#' Get the minimum area
#'
#' @inheritParams default_params_doc
#'
#' @return Numeric with time at which area is minimum during the simulation
#' @noRd
#'
#' @author Pedro Neves, Joshua Lambert, Shu Xie
get_global_min_area <- function(total_time,
                                area_pars,
                                peak,
                                island_ontogeny,
                                sea_level) {
  fx <- function(timeval) {
    y <- island_area(
      timeval,
      total_time = total_time,
      area_pars = area_pars,
      peak = peak,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level
    )
    if (is.nan(y)) {
      return(Inf)
    } else {
      return(y)
    }
  }
  global_min_area_time <- subplex::subplex(par = 0, fn = fx)$par
  # testit::assert(is.numeric((global_min_area_time)))
  global_min_area_time <- DDD::roundn(global_min_area_time, 14)

  Amin <- island_area(
    timeval = global_min_area_time,
    total_time = total_time,
    area_pars = area_pars,
    peak = peak,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level
  )
  return(Amin)
}

#' Piecewise area bounds for the thinning of the time-dependent simulation
#'
#' @description Internal function. The time-dependent simulation draws
#' candidate events from an upper bound of the total event rate and accepts
#' them with probability rate / bound (thinning). Using the global minimum and
#' maximum area over the whole simulation for this bound is very wasteful, in
#' particular for island ontogeny, where the area is zero at the start of the
#' island so that the extinction bound is always `extcutoff`. This function
#' splits the simulation time in windows and returns the minimum and maximum
#' area within each window, evaluated on a grid and widened by a small safety
#' margin and clipped to the global bounds.
#'
#' @inheritParams default_params_doc
#' @param n_grid Integer, number of grid intervals on which the area is
#'   evaluated.
#' @param ratio Numeric, a window is extended as long as its maximum area is
#'   at most `ratio` times its minimum area.
#' @param max_windows Integer, upper limit of the number of windows.
#' @param margin Numeric, fraction of the area range within a window by which
#'   the bounds are widened.
#'
#' @return A list with the vectors `end` (end time of each window), `Amin` and
#' `Amax` (area bounds in each window).
#' @keywords internal
#' @noRd
get_area_windows <- function(total_time,
                             area_pars,
                             peak,
                             island_ontogeny,
                             sea_level,
                             Amin,
                             Amax,
                             n_grid = 1000,
                             ratio = 1.5,
                             max_windows = 200,
                             margin = 0.02) {
  freq <- area_pars$sea_level_frequency
  if (is.null(freq) || !is.finite(freq)) freq <- 0
  cycles <- freq * total_time / area_pars$total_island_age
  n_grid <- max(n_grid, ceiling(50 * cycles))
  tt <- seq(0, total_time, length.out = n_grid + 1)
  a <- vapply(tt, island_area, numeric(1),
              total_time = total_time, area_pars = area_pars,
              peak = peak, island_ontogeny = island_ontogeny,
              sea_level = sea_level)
  a[!is.finite(a)] <- NA_real_
  tiny <- 1e-3 * Amax

  build <- function(ratio) {
    starts <- 1L
    i <- 1L
    n <- length(tt)
    while (i < n) {
      j <- i + 1L
      lo <- hi <- a[i]
      while (j <= n) {
        lo2 <- min(lo, a[j], na.rm = TRUE)
        hi2 <- max(hi, a[j], na.rm = TRUE)
        if (hi2 > max(ratio * lo2, tiny)) break
        lo <- lo2
        hi <- hi2
        j <- j + 1L
      }
      i <- max(i + 1L, j - 1L)
      starts <- c(starts, i)
    }
    unique(c(starts, n))
  }
  idx <- build(ratio)
  while (length(idx) - 1 > max_windows) {
    ratio <- ratio * 1.5
    idx <- build(ratio)
  }
  nw <- length(idx) - 1
  amin <- amax <- numeric(nw)
  for (w in seq_len(nw)) {
    aw <- a[idx[w]:idx[w + 1]]
    aw <- aw[!is.na(aw)]
    if (length(aw) == 0) {
      amin[w] <- Amin
      amax[w] <- Amax
    } else {
      pad <- margin * (max(aw) - min(aw))
      amin[w] <- max(Amin, min(aw) - pad)
      amax[w] <- min(Amax, max(aw) + pad)
    }
  }
  list(end = tt[idx[-1]], Amin = amin, Amax = amax)
}
