#' Internal function of the DAISIE simulation
#'
#' @param time Simulated amount of time
#' @param mainland_n A numeric stating the number of mainland species, that
#' is the number of species that can potentially colonize the island.
#' If \code{\link{DAISIE_sim}} uses a clade-specific diversity dependence,
#' this value is set to 1.
#' If \code{\link{DAISIE_sim}} uses an island-wide diversity dependence,
#' this value is set to the number of mainland species.
#' @param pars A numeric vector:
#' \itemize{
#'   \item{[1]: cladogenesis rate}
#'   \item{[2]: extinction rate}
#'   \item{[3]: carrying capacity}
#'   \item{[4]: immigration rate}
#'   \item{[5]: anagenesis rate}
#' }
#' @param ddmodel_sim The a numeric vector to determined which parameters should
#' be diversity dependent.
#' @param island_type Option island_type = 'oceanic' is a model equal to Valente
#' et al., 2015. island_type = 'nonoceanic' is a nonoceanic model where initial
#' species richness is non-zero determined by the nonoceanic parameters.
#' @param nonoceanic_pars A vector of length three with: the island area as a
#' proportion of the mainland, the probability of native species being
#' nonendemic and the size of the mainland pool.
#' @param area_pars A named list containing area parameters as created by
#' create_area_pars:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#' }
#' @param ext_pars A numeric vector:
#' \itemize{
#'   \item{[1]: minimum extinction when area is at peak}
#'   \item{[2]: extinction rate when current area is 0.10 of maximum area}
#' }
#' @param k_dist_pars Vector of two elements to define gamma
#' distribution for sampling carrying capacities. \code{k_dist_pars[1]}
#' is the shape parameter, \code{k_dist_pars[2]} is the rate parameter
#' @param island_ontogeny a numeric describing the type of island ontogeny.
#' Can be \code{0} for constant, \code{1} for a beta function describing area.
#' @param ext_multiplier Numeric between 0 and 1 reducing effective extinction
#' rate for the calculation of t_hor. Default is 0.5.
#' @param sea_level a string describing the type of sea level.
#' Can be \code{"const"} or \code{"sine"} for a sine function describing area
#' @param shift_times a numeric vector specifying when the rate shifts occur
#' before the present.
#' @param hyper_pars A numeric vector for hyperparameters for the rate
#' calculations, \code{hyper_pars[1]} is d_0 the scaling parameter for
#' exponent for calculating cladogenesis rate, \code{hyper_pars[2]}
#' is x the exponent for calculating extinction rate,
#' \code{hyper_pars[3]} is alpha the exponent for calculating the
#' immigration rate, \code{hyper_pars[4]} is beta the exponent for
#' calculating the anagenesis rate.
#' @param dist_pars a numeric for the distance from the mainland.
DAISIE_sim_core <- function(
  time,
  mainland_n,
  pars,
  ddmodel_sim = 11,
  island_type = "oceanic",
  nonoceanic_pars = NULL,
  k_dist_pars = NULL,
  island_ontogeny = 0,
  sea_level = 0,
  hyper_pars = NULL,
  area_pars = NULL,
  dist_pars = NULL,
  ext_pars = NULL,
  shift_times = NULL,
  ext_multiplier = 0.5
) {
  timeval <- 0
  totaltime <- time
  island_ontogeny <- translate_island_ontogeny(island_ontogeny)
  sea_level <- translate_sea_level(sea_level)

  testit::assert(length(pars) == 5 || length(pars) == 10)
  if (!is.null(area_pars) &&
      (island_ontogeny == 0 && sea_level == 0)) {
    stop("area_pars specified for constant island_ontogeny and sea_level.
         Set area_pars to NULL.")
  }
  testit::assert(is.null(area_pars) || are_area_pars(area_pars))
  if (pars[4] == 0 && island_type == "oceanic") {
    stop("Island has no species and the rate of
    colonisation is zero. Island cannot be colonised.")
  }
  if ((is.null(ext_pars) || is.null(area_pars)) &&
      (island_ontogeny != 0 || sea_level != 0)) {
    stop("Island ontogeny and/or sea level specified but area parameters
    and/or extinction parameters not available. Please either set
    island_ontogeny and sea_level to NULL, or specify area_pars and ext_pars.")
  }

  if (island_ontogeny == 0 && sea_level == 0) {
    area_pars <- create_area_pars(
      max_area = 1,
      proportional_peak_t = 0,
      peak_sharpness = 0,
      total_island_age = totaltime,
      sea_level_amplitude = 0,
      sea_level_frequency = 0
    )
  }

  testit::assert(are_area_pars(area_pars))
  testit::assert((totaltime <= area_pars$total_island_age) ||
                   is.null(area_pars))

  if (island_type == "nonoceanic") {
    nonoceanic_sample <- DAISIE_nonoceanic_spec(
      prob_samp = nonoceanic_pars[1],
      prob_nonend = nonoceanic_pars[2],
      mainland_n = mainland_n)
    init_nonend_spec_vec <- nonoceanic_sample[[1]]
    init_end_spec_vec <- nonoceanic_sample[[2]]
    mainland_spec <- nonoceanic_sample[[3]]
  }
  if (island_type == "oceanic") {
    mainland_spec <- seq(1, mainland_n, 1)
    init_nonend_spec <- 0
    init_end_spec <- 0
  }
  maxspecID <- mainland_n

  if (!is.null(k_dist_pars)) {
    K <- stats::rgamma(1, shape = k_dist_pars[[1]], rate = k_dist_pars[[2]])
  }
  island_spec <- c()
  stt_table <- matrix(ncol = 4)
  colnames(stt_table) <- c("Time", "nI", "nA", "nC")
  if (island_type == "oceanic") {
    stt_table[1, ] <- c(totaltime, 0, 0, 0)
  } else {
    nonoceanic_tables <- DAISIE_nonoceanic_stt_table(stt_table,
                                                     totaltime,
                                                     timeval,
                                                     init_nonend_spec_vec,
                                                     init_end_spec_vec,
                                                     mainland_spec,
                                                     island_spec
    )
    stt_table <- nonoceanic_tables$stt_table
    init_nonend_spec <- nonoceanic_tables$init_nonend_spec
    init_end_spec <- nonoceanic_tables$init_end_spec
    mainland_spec <- nonoceanic_tables$mainland_spec
    island_spec <- nonoceanic_tables$island_spec
  }

  lac <- pars[1]
  mu <- pars[2]
  K <- pars[3]
  gam <- pars[4]
  laa <- pars[5]
  extcutoff <- max(1000, 1000 * (laa + lac + gam))
  testit::assert(is.numeric(extcutoff))
  land_bridge <- land_bridge_periods(0,
                                     totaltime,
                                     shift_times)
  if (land_bridge$present == TRUE) {
    lac <- pars[6]
    mu <- pars[7]
    K <- pars[8]
    gam <- pars[9]
    laa <- pars[10]
    extcutoff <- max(1000, 1000 * (laa + lac + gam))
    testit::assert(is.numeric(extcutoff))
  }

  num_spec <- length(island_spec[, 1])
  num_immigrants <- length(which(island_spec[, 4] == "I"))

  max_rates <- update_max_rates(
    timeval = timeval,
    totaltime = totaltime,
    gam = gam,
    mu = mu,
    laa = laa,
    lac = lac,
    ddmodel_sim = ddmodel_sim,
    hyper_pars = hyper_pars,
    area_pars = area_pars,
    dist_pars = dist_pars,
    ext_pars = ext_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    extcutoff = extcutoff,
    K = K,
    num_spec = num_spec,
    num_immigrants = num_immigrants,
    mainland_n = mainland_n
  )
  timeval_and_dt <- calc_next_timeval(
    max_rates = max_rates,
    timeval = timeval
  )
  timeval <- timeval_and_dt$timeval

  #### Start Gillespie ####
  while (timeval < totaltime) {

    land_bridge <- land_bridge_periods(timeval,
                                       totaltime,
                                       shift_times)

    if (land_bridge$present == TRUE) {
      lac <- pars[6]
      mu <- pars[7]
      K <- pars[8]
      gam <- pars[9]
      laa <- pars[10]
      extcutoff <- max(1000, 1000 * (laa + lac + gam))
      testit::assert(is.numeric(extcutoff))
    }

    rates <- update_rates(
      timeval = timeval,
      totaltime = totaltime,
      gam = gam,
      mu = mu,
      laa = laa,
      lac = lac,
      ddmodel_sim = ddmodel_sim,
      hyper_pars = hyper_pars,
      area_pars = area_pars,
      dist_pars = dist_pars,
      ext_pars = ext_pars,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      extcutoff = extcutoff,
      K = K,
      num_spec = num_spec,
      num_immigrants = num_immigrants,
      mainland_n = mainland_n
    )
    testit::assert(are_rates(rates))

    possible_event <- DAISIE_sample_event(
      rates = rates,
      max_rates = max_rates
    )

    updated_state <- DAISIE_sim_update_state(
      timeval = timeval,
      totaltime = totaltime,
      possible_event = possible_event,
      maxspecID = maxspecID,
      mainland_spec = mainland_spec,
      island_spec = island_spec,
      stt_table = stt_table
    )

    island_spec <- updated_state$island_spec
    maxspecID <- updated_state$maxspecID
    stt_table <- updated_state$stt_table
    num_spec <- length(island_spec[, 1])
    num_immigrants <- length(which(island_spec[, 4] == "I"))

    max_rates <- update_max_rates(
      timeval = timeval,
      totaltime = totaltime,
      gam = gam,
      mu = mu,
      laa = laa,
      lac = lac,
      ddmodel_sim = ddmodel_sim,
      hyper_pars = hyper_pars,
      area_pars = area_pars,
      dist_pars = dist_pars,
      ext_pars = ext_pars,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      extcutoff = extcutoff,
      K = K,
      num_spec = num_spec,
      num_immigrants = num_immigrants,
      mainland_n = mainland_n
    )
    timeval_and_dt <- calc_next_timeval(max_rates = max_rates, timeval = timeval)
    timeval <- timeval_and_dt$timeval

    land_bridge_at_next_dt <- land_bridge_periods(timeval,
                                                  totaltime,
                                                  shift_times)

    if (land_bridge$shift_time != land_bridge_at_next_dt$shift_time) {
      timeval <- land_bridge$shift_time
    }
  }
  #### Finalize STT ####
  stt_table <- rbind(
    stt_table,
    c(
      0,
      stt_table[nrow(stt_table), 2],
      stt_table[nrow(stt_table), 3],
      stt_table[nrow(stt_table), 4]
    )
  )
  island <- DAISIE_create_island(
    stt_table = stt_table,
    totaltime = totaltime,
    island_spec = island_spec,
    mainland_n = mainland_n,
    init_nonend_spec = init_nonend_spec,
    init_end_spec = init_end_spec,
    carrying_capacity = K)
  return(island)
}
