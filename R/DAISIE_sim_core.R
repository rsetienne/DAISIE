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
#' @param nonoceanic_params A vector of length three with: the island area as a
#' proportion of the mainland, the probability of native species being
#' nonendemic and the size of the mainland pool.
#' @param Apars A named list containing area parameters as created by
#' create_area_params:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#' }
#' @param Epars A numeric vector:
#' \itemize{
#'   \item{[1]: minimum extinction when area is at peak}
#'   \item{[2]: extinction rate when current area is 0.10 of maximum area}
#' }
#' @param keep_final_state Logical indicating if final state of simulation
#' should be returned. Default is \code{FALSE}.
#' @param island_spec A matrix with species on island (state of system
#' at each time point).
DAISIE_sim_core <- function(
  time,
  mainland_n,
  pars,
  ddmodel_sim = 11,
  island_type = "oceanic",
  nonoceanic_params = NULL,
  k_dist_params = NULL,
  island_ontogeny = 0,
  Apars = NULL,
  Epars = NULL,
  keep_final_state = FALSE,
  island_spec = NULL
) {
  testit::assert(is.logical(keep_final_state))
  testit::assert(length(pars) == 5)
  testit::assert(is.null(Apars) || are_area_params(Apars))
  # testit::assert(is.null(island_spec) || is.matrix(island_spec))
  if (pars[4] == 0 && island_type == "oceanic") {
    stop("Island has no species and the rate of
    colonisation is zero. Island cannot be colonised.")
  }
  if (!is.null(Apars) && island_ontogeny == "const") {
    stop("Apars specified for constant island_ontogeny. Set Apars to NULL.")
  }
  timeval <- 0
  totaltime <- time
  lac <- pars[1]
  mu <- pars[2]
  K <- pars[3]
  gam <- pars[4]
  laa <- pars[5]
  extcutoff <- max(1000, 1000 * (laa + lac + gam))
  testit::assert(is.numeric(extcutoff))
  ext_multiplier <- 0.5
  testit::assert((totaltime <= Apars$total_island_age) || is.null(Apars))
  # Make island_ontogeny be numeric
  island_ontogeny <- translate_island_ontogeny(island_ontogeny)
  if ((is.null(Epars) || is.null(Apars)) && (island_ontogeny != 0)) {
    stop ("Island ontogeny specified but Area parameters and/or extinction
         parameters not available. Please either set island_ontogeny to NULL, or
         specify Apars and Epars.")
  }
  if (island_type == "nonoceanic") {
    nonoceanic_sample <- DAISIE_nonoceanic_spec(prob_samp = nonoceanic_params[1],
                                                prob_nonend = nonoceanic_params[2],
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

  if (!is.null(k_dist_params)) {
      K <- rgamma(1, shape = k_dist_params[[1]], rate = k_dist_params[[2]])
  }

  #### Start Gillespie ####
  # Start output and tracking objects
  if (is.null(island_spec)) {
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
                                                       island_spec)
     stt_table <- nonoceanic_tables$stt_table
     init_nonend_spec <- nonoceanic_tables$init_nonend_spec
     init_end_spec <- nonoceanic_tables$init_end_spec
     mainland_spec <- nonoceanic_tables$mainland_spec
     island_spec <- nonoceanic_tables$island_spec
    }
  }
    #if starting using keep_final_state
    if (keep_final_state == TRUE) {
      # stt_table <- matrix(stt_table[nrow(stt_table), ], nrow = 1, ncol = 4)
      stt_table <- matrix(ncol = 4)
      colnames(stt_table) <- c("Time", "nI", "nA", "nC")
      stt_table[1, 1] <- totaltime
      stt_table[1, 2] <- length(which(island_spec[, 4] == "I"))
      stt_table[1, 3] <- length(which(island_spec[, 4] == "A"))
      stt_table[1, 4] <- length(which(island_spec[, 4] == "C"))
    }
  testit::assert(is.null(Apars) || are_area_params(Apars))
  # Pick t_hor (before timeval, to set Amax t_hor)
  t_hor <- get_t_hor(
    timeval = 0,
    totaltime = totaltime,
    Apars = Apars,
    ext = 0,
    ext_multiplier = ext_multiplier,
    island_ontogeny = island_ontogeny,
    t_hor = NULL)
  while (timeval < totaltime) {
    # Calculate rates
    rates <- update_rates(
      timeval = timeval,
      totaltime = totaltime,
      gam = gam,
      mu = mu,
      laa = laa,
      lac = lac,
      ddmodel_sim = ddmodel_sim,
      Apars = Apars,
      Epars = Epars,
      island_ontogeny = island_ontogeny,
      extcutoff = extcutoff,
      K = K,
      island_spec = island_spec,
      mainland_n = mainland_n,
      t_hor = t_hor)
    timeval_and_dt <- calc_next_timeval(rates, timeval)
    timeval <- timeval_and_dt$timeval
    dt <- timeval_and_dt$dt
    if (timeval <= t_hor) {
      testit::assert(are_rates(rates))
      # Determine event
      possible_event <- DAISIE_sample_event(
        rates = rates,
        island_ontogeny = island_ontogeny)

      updated_state <- DAISIE_sim_update_state(
        timeval = timeval,
        totaltime = totaltime,
        possible_event = possible_event,
        maxspecID = maxspecID,
        mainland_spec = mainland_spec,
        island_spec = island_spec,
        stt_table = stt_table)
      island_spec <- updated_state$island_spec
      maxspecID <- updated_state$maxspecID
      stt_table <- updated_state$stt_table
    } else {
      #### After t_hor is reached ####
      timeval <- t_hor
      t_hor <- get_t_hor(
        timeval = timeval,
        totaltime = totaltime,
        Apars = Apars,
        ext = rates$ext_rate,
        ext_multiplier = ext_multiplier,
        island_ontogeny = island_ontogeny,
        t_hor = t_hor)
    }
    # TODO Check if this is redundant, or a good idea
    if (rates$ext_rate_max >= extcutoff && length(island_spec[, 1]) == 0) {
      timeval <- totaltime
    }
  }
  # Finalize stt_table
  stt_table <- rbind(
    stt_table,
    c(0,
      stt_table[nrow(stt_table), 2],
      stt_table[nrow(stt_table), 3],
      stt_table[nrow(stt_table), 4])
    )
    island <- DAISIE_create_island(
      stt_table = stt_table,
      totaltime = totaltime,
      island_spec = island_spec,
      mainland_n = mainland_n,
      keep_final_state = keep_final_state,
      init_nonend_spec = init_nonend_spec,
      init_end_spec = init_end_spec,
      carrying_capacity = K)
    return(island)
}
