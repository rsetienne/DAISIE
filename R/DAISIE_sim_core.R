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
#' @param nonoceanic_pars A vector of length three with: the island area as a
#' proportion of the mainland, the probability of native species being
#' nonendemic and the size of the mainland pool.
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
#' @param ext_pars A numeric vector:
#' \itemize{
#'   \item{[1]: minimum extinction when area is at peak}
#'   \item{[2]: extinction rate when current area is 0.10 of maximum area}
#' }
#' @param island_ontogeny a numeric describing the type of island ontogeny.
#' Can be \code{0} for constant, \code{1} for a beta function describing area.
#' rate for the calculation of t_hor. Default is 0.5.
#' @param sea_level a numeric describing the type of sea level.
#' Can be \code{0} or \code{1} for a sine function describing area
#' @param shift_times a numeric vector specifying when the rate shifts occur
#' before the present.
#' @param hyper_pars A numeric vector for hyperparameters for the rate
#' calculations:
#' \itemize{
#' \item{[1]: is d_0 the scaling parameter for exponent for calculating
#' cladogenesis rate}
#' \item{[2]: is x the exponent for calculating extinction rate}
#' \item{[3]: is alpha, the exponent for calculating the immigration rate}
#' \item{[4]: is beta the exponent for calculating the anagenesis rate.}
#' }
#' @param extcutoff the maximum per capita extinction rate.
#' @param dist_pars a numeric for the distance from the mainland.
DAISIE_sim_core <- function(
  time,
  mainland_n,
  pars,
  nonoceanic_pars = c(0, 0),
  island_ontogeny = 0,
  sea_level = 0,
  hyper_pars = NULL,
  area_pars = NULL,
  dist_pars = NULL,
  ext_pars = NULL,
  extcutoff = 1000,
  shift_times = NULL
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
  if (pars[4] == 0 && nonoceanic_pars[1] == 0) {
    stop("Island has no species and the rate of
    colonisation is zero. Island cannot be colonised.")
  }
  if ((is.null(ext_pars) || is.null(area_pars)) &&
      (island_ontogeny != 0 || sea_level != 0)) {
    stop("Island ontogeny and/or sea level specified but area parameters
    and/or extinction parameters not available. Please either set
    island_ontogeny and sea_level to NULL, or specify area_pars and ext_pars.")
  }
  default_metapars <- create_default_pars(
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    area_pars = area_pars,
    hyper_pars = hyper_pars,
    dist_pars = dist_pars,
    ext_pars = ext_pars,
    totaltime = totaltime,
    pars = pars)
  hyper_pars <- default_metapars$hyper_pars
  dist_pars <- default_metapars$dist_pars
  ext_pars <- default_metapars$ext_pars
  area_pars <- default_metapars$area_pars

  testit::assert(are_hyper_pars(hyper_pars = hyper_pars))
  testit::assert(are_area_pars(area_pars = area_pars))
  testit::assert(are_dist_pars(dist_pars = dist_pars))
  testit::assert((totaltime <= area_pars$total_island_age) ||
                   is.null(area_pars))
  nonoceanic_sample <- DAISIE_nonoceanic_spec(
    prob_samp = nonoceanic_pars[1],
    prob_nonend = nonoceanic_pars[2],
    mainland_n = mainland_n)
  maxspecID <- mainland_n
  island_spec <- c()
  stt_table <- matrix(ncol = 4)
  colnames(stt_table) <- c("Time", "nI", "nA", "nC")
  spec_tables <- DAISIE_spec_tables(stt_table,
                                    totaltime,
                                    timeval,
                                    nonoceanic_sample,
                                    island_spec)
    stt_table <- spec_tables$stt_table
    mainland_spec <- spec_tables$mainland_spec
    island_spec <- spec_tables$island_spec

  initial_land_bridge <- land_bridge_periods(0,
                                             totaltime,
                                             shift_times)

  if (initial_land_bridge$present == FALSE) {
    lac <- pars[1]
    mu <- pars[2]
    K <- pars[3]
    gam <- pars[4]
    laa <- pars[5]
  } else {
    lac <- pars[6]
    mu <- pars[7]
    K <- pars[8]
    gam <- pars[9]
    laa <- pars[10]
  }
  testit::assert(is.numeric(extcutoff))

  num_spec <- length(island_spec[, 1])
  num_immigrants <- length(which(island_spec[, 4] == "I"))

  max_rates <- update_max_rates(
    timeval = timeval,
    totaltime = totaltime,
    gam = gam,
    laa = laa,
    lac = lac,
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

  if (timeval < totaltime) {
    initial_land_bridge_plus_dt <- land_bridge_periods(timeval,
                                                       totaltime,
                                                       shift_times)

    if (!is.character(initial_land_bridge$shift_time) &&
        (initial_land_bridge$shift_time !=
         initial_land_bridge_plus_dt$shift_time)) {
      timeval <- initial_land_bridge_plus_dt$shift_time

      if (initial_land_bridge$present == FALSE) {
        lac <- pars[1]
        mu <- pars[2]
        K <- pars[3]
        gam <- pars[4]
        laa <- pars[5]
      } else {
        lac <- pars[6]
        mu <- pars[7]
        K <- pars[8]
        gam <- pars[9]
        laa <- pars[10]
      }

      max_rates <- update_max_rates(
        timeval = timeval,
        totaltime = totaltime,
        gam = gam,
        laa = laa,
        lac = lac,
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
    }
  }
  #### Start Gillespie ####
  while (timeval < totaltime) {
    rates <- update_rates(
      timeval = timeval,
      totaltime = totaltime,
      gam = gam,
      laa = laa,
      lac = lac,
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
      laa = laa,
      lac = lac,
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

    land_bridge_before_dt <- land_bridge_periods(timeval,
                                                 totaltime,
                                                 shift_times)

    timeval_and_dt <- calc_next_timeval(max_rates = max_rates,
                                        timeval = timeval)
    timeval <- timeval_and_dt$timeval

    if (timeval < totaltime) {
      land_bridge <- land_bridge_periods(timeval,
                                         totaltime,
                                         shift_times)

      if (!is.character(land_bridge$shift_time) &&
          (land_bridge_before_dt$shift_time != land_bridge$shift_time)) {
        timeval <- land_bridge$shift_time
        if (land_bridge$present == FALSE) {
          lac <- pars[1]
          mu <- pars[2]
          K <- pars[3]
          gam <- pars[4]
          laa <- pars[5]
        } else {
          lac <- pars[6]
          mu <- pars[7]
          K <- pars[8]
          gam <- pars[9]
          laa <- pars[10]
        }
        max_rates <- update_max_rates(
          timeval = timeval,
          totaltime = totaltime,
          gam = gam,
          laa = laa,
          lac = lac,
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
      }
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
    mainland_n = mainland_n)
  return(island)
}
