#' @title Simulate (non-)oceanic islands with given parameters under a
#'   rate-shift regime
#'
#' @description
#' This function simulates islands with given cladogenesis,
#' extinction, Kprime, immigration and anagenesis parameters, all of which
#' modelled as time-constant parameters, which can be switched to a different
#' diversification regime (i.e., different set of parameters) at one or more
#' set times before the present. Further, it allows for the simulation of
#' non-oceanic islands, generating islands for which the starting condition
#' includes potential endemic and non-endemic species.
#'
#' @inheritParams default_params_doc
#'
#' @return
#' A list. The highest level of the least corresponds to each individual
#' replciate. The first element of each replicate is composed of island
#' information containing:
#' \itemize{
#'   \item{\code{$island_age}: A numeric with the island age.}
#'   \item{\code{$not_present}: A numeric with the number of mainland lineages
#'   that are not present on the island.}
#'   \item{\code{$stt_all}: STT table for all species on the island
#'     (nI - number of non-endemic species; nA - number of anagenetic species,
#'     nC - number of cladogenetic species, present - number of independent
#'     colonisations present)}
#'   \item{\code{$brts_table}: Only for simulations under \code{"IW"}. Table
#' containing information on order of events in the data, for use in maximum
#' likelihood optimization.).}
#' }
#' The subsequent elements of the list pertaining to each replcate contain
#' information on a single colonist lineage on the island and have 4 components:
#' \itemize{
#'   \item{\code{$branching_times}: island age and stem age of the
#'     population/species in the case of Non-endemic, Non-endemic_MaxAge and
#'     Endemic anagenetic species. For cladogenetic species these should
#'     be island age and branching times of the radiation including the
#'     stem age of the radiation.}
#'   \item{\code{$stac}: An integer ranging from 1 to 4
#'   indicating the status of the colonist:}
#'   \enumerate{
#'     \item Non_endemic_MaxAge
#'     \item Endemic
#'     \item Endemic&Non_Endemic
#'     \item Non_endemic_MaxAge
#' }
#' \item{\code{$missing_species}: number of island species that were
#' not sampled for particular clade (only applicable for endemic clades)}
#' \item{\code{$type_1or2}: whether the colonist belongs to type 1 or type 2}
#' }
#' @author Luis Valente, Albert Phillimore, Torsten Hauffe
#' @seealso \code{\link{DAISIE_plot_sims}()} for plotting STT of simulation
#' outputs.
#' @family simulation models
#' @references Valente, L.M., A.B. Phillimore and R.S. Etienne (2015).
#' Equilibrium and non-equilibrium dynamics simultaneously operate in the
#' Galapagos islands. Ecology Letters 18: 844-852.
#'
#' Hauffe, T., D. Delicado, R.S. Etienne and L. Valente (2020).
#' Lake expansion elevates equilibrium diversity via increasing colonization.
#' @keywords models
#' @export
DAISIE_sim_cr_shift <- function(
  time,
  M,
  pars,
  replicates,
  shift_times,
  divdepmodel = "CS",
  nonoceanic_pars = c(0, 0),
  num_guilds = NULL,
  sample_freq = 25,
  plot_sims = TRUE,
  hyper_pars = create_hyper_pars(d = 0, x = 0),
  area_pars = DAISIE::create_area_pars(
    max_area = 1,
    current_area = 1,
    proportional_peak_t = 0,
    total_island_age = 0,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0),
  cond = 0,
  verbose = TRUE,
  ...
) {
  testit::assert(
    "length(pars) is not ten, set ten parameters",
    length(pars) == 10
  )
  testit::assert(are_hyper_pars(hyper_pars = hyper_pars))
  testit::assert(are_area_pars(area_pars = area_pars))

  total_time <- time
  island_replicates <- list()
  if (divdepmodel == "IW") {
    for (rep in 1:replicates) {
      island_replicates[[rep]] <- DAISIE_sim_core_cr_shift(
        time = total_time,
        mainland_n = M,
        pars = pars,
        nonoceanic_pars = nonoceanic_pars,
        hyper_pars = hyper_pars,
        area_pars = area_pars,
        shift_times = shift_times
      )
      if (verbose == TRUE) {
        message("Island replicate ", rep)
      }
    }
    island_replicates <- DAISIE_format_IW(island_replicates = island_replicates,
                                          time = total_time,
                                          M = M,
                                          sample_freq = sample_freq,
                                          verbose = verbose)
  }
  if (divdepmodel == "CS") {
    for (rep in 1:replicates) {
      island_replicates[[rep]] <- list()
      full_list <- list()
      if (cond == 0) {
        number_present <- -1
      } else {
        number_present <- 0
      }
      while (number_present < cond) {
        for (m_spec in 1:M) {
          full_list[[m_spec]] <- DAISIE_sim_core_cr_shift(
            time = total_time,
            mainland_n = 1,
            pars = pars,
            nonoceanic_pars = nonoceanic_pars,
            hyper_pars = hyper_pars,
            area_pars = area_pars,
            shift_times = shift_times
          )
        }
        stac_vec <- unlist(full_list)[which(names(unlist(full_list)) == "stac")]
        present <- which(stac_vec != 0)
        number_present <- length(present)
      }
      island_replicates[[rep]] <- full_list
      if (verbose == TRUE) {
        message("Island replicate ", rep)
      }
    }
    island_replicates <- DAISIE_format_CS(
      island_replicates = island_replicates,
      time = total_time,
      M = M,
      sample_freq = sample_freq,
      verbose = verbose
    )
  }

  if (divdepmodel == "GW") {
    if (!is.numeric(num_guilds)) {
      stop("num_guilds must be numeric")
    }
    guild_size <- M / num_guilds
    testit::assert(num_guilds < M)
    testit::assert(M %% num_guilds == 0)
    for (rep in 1:replicates) {
      island_replicates[[rep]] <- list()
      full_list <- list()
      for (m_spec in 1:num_guilds) {
        full_list[[m_spec]]  <- DAISIE_sim_core_cr_shift(
          time = total_time,
          mainland_n = guild_size,
          pars = pars,
          nonoceanic_pars = nonoceanic_pars,
          hyper_pars = hyper_pars,
          area_pars = area_pars,
          shift_times = shift_times
        )
      }
      island_replicates[[rep]] <- full_list
      if (verbose == TRUE) {
        message("Island replicate ", rep)
      }
    }
    island_replicates <- DAISIE_format_GW(island_replicates = island_replicates,
                                          time = total_time,
                                          M = M,
                                          sample_freq = sample_freq,
                                          num_guilds = num_guilds,
                                          verbose = verbose)
  }
  if (plot_sims == TRUE) {
    DAISIE_plot_sims(
      island_replicates = island_replicates,
      sample_freq = sample_freq
    )
  }
  return(island_replicates)
}
