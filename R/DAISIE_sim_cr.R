#' @title Simulate (non-)oceanic islands with given parameters under
#'   time-constant rates
#' @name DAISIE_sim
#' @aliases DAISIE_sim_cr DAISIE_sim
#'
#' @description
#' This function simulates islands with given cladogenesis,
#' extinction, Kprime, immigration and anagenesis parameters, all of which
#' modelled as time-constant parameters. If a single
#' parameter set is provided (5 parameters) it simulates islands where all
#' species have the same macro-evolutionary process. If two paramater sets
#' (10 parameters) are provided, it simulates islands where two different
#' macro-evolutionary processes operate, one applying to type 1 species and
#' other to type 2 species. Further, it allows for the simulation of
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
#'   \item{\code{$not_present}: the number of mainland lineages that are not
#'     present on the island. It is only present if only 1 type of species is
#'     simulated. Becomes \code{$not_present_type1}: the number of mainland
#'     lineages of type 1 that are not present on the island and
#'     \code{$not_present_type2}: the number of mainland lineages of type 2
#'     that are not present on the island, if two types are simulated.}
#'   \item{\code{$stt_all}: STT table for all species on the island
#'     (nI - number of non-endemic species; nA - number of anagenetic species,
#'     nC - number of cladogenetic species, present - number of independent
#'     colonisations present)}
#'   \item{\code{$stt_stt_type1}: STT table for type 1 species on the island -
#'     only if 2 types of species were simulated (nI - number of non-endemic
#'     species; nA - number of anagenetic species, nC - number of cladogenetic
#'     species, present - number of independent colonisations present).}
#'   \item{\code{$stt_stt_type2}: STT table for type 2 species on the island
#'      - only if 2 types of species were simulated (nI - number of non-endemic
#'      species; nA - number of anagenetic species, nC - number of cladogenetic
#'      species, present - number of independent colonisations present ).}
#'   \item{\code{$brts_table}: Only for simulations under \code{"IW"}. Table
#' containing information on order of events in the data, for use in maximum
#' likelihood optimization.).}
#' }
#' The subsequent elements of the list pertaining to each replcate contain
#' information on a single colonist lineage on the island and have 4 components:
#' \itemize{
#'   \item{\code{$branching_times}: island age and stem age of the
#'     population/species in the case of Non-endemic, Non-endemic_MaxAge and
#'     Endemic anagenetic species.
#'
#'     For cladogenetic species these should
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
#' @author Luis Valente, Albert Phillimore, Joshua Lambert, Shu Xie, Pedro
#' Neves, Richèl J. C. Bilderbeek, Rampal Etienne
#' @seealso \code{\link{DAISIE_plot_sims}()} for plotting STT of simulation
#' outputs.
#' @family simulation models
#' @references Valente, L.M., A.B. Phillimore and R.S. Etienne (2015).
#' Equilibrium and non-equilibrium dynamics simultaneously operate in the
#' Galapagos islands. Ecology Letters 18: 844-852.
#' @keywords models
#' @examples
#' ## Simulate 2 islands for 1 million years, where all species have equal
#' ## rates. Pool size 100.
#'
#' clado_rate <- 0.5
#' ext_rate <- 0.2
#' carr_cap <- Inf
#' immig_rate <- 0.05
#' ana_rate <- 1
#' sim_pars <- c(clado_rate, ext_rate, carr_cap, immig_rate, ana_rate)
#' set.seed(1)
#' island_replicates <- DAISIE_sim_cr(
#'   time = 1,
#'   M = 100,
#'   pars = sim_pars,
#'   replicates = 2,
#'   plot_sims = FALSE,
#'   verbose = FALSE
#' )
#'
#' ## Simulate 2 islands for 1 million years with two types of species (type1
#' ## and type 2). Pool size 100
#' ## Fraction of type 2 species in source pool is 0.15. Function will
#' ## simulate until number of islands where type 2 species has colonised is
#' ## equal to number specified in replicates.
#'
#' clado_rate <- 0.5
#' ext_rate <- 0.2
#' carr_cap <- Inf
#' immig_rate <- 0.005
#' ana_rate <- 1
#' sim_pars_type1 <- c(clado_rate, ext_rate, carr_cap, immig_rate, ana_rate)
#' sim_pars_type2 <- sim_pars_type1 * 2
#' set.seed(1)
#' island_replicates_2types <- DAISIE_sim_cr(
#'   time = 1,
#'   M = 100,
#'   pars = c(sim_pars_type1, sim_pars_type2),
#'   replicates = 2,
#'   prop_type2_pool = 0.15,
#'   plot_sims = FALSE,
#'   verbose = FALSE
#' )
#' ## Simulate two non-oceanic island for 1 million years.
#' ## Pool size 500. Island area as a proportion
#' ## of mainland is 0.1, proportion of native species is 0.9.
#' clado_rate <- 0.5
#' ext_rate <- 0.2
#' carr_cap <- Inf
#' immig_rate <- 0.005
#' ana_rate <- 1
#' sim_pars <- c(clado_rate, ext_rate, carr_cap, immig_rate, ana_rate)
#' set.seed(1)
#' island_replicates <- DAISIE_sim_cr(
#'   time = 1,
#'   M = 500,
#'   pars = sim_pars,
#'   replicates = 2,
#'   nonoceanic_pars = c(0.1, 0.9),
#'   plot_sims = FALSE,
#'   verbose = FALSE
#' )
#'
#' ## Simulate 2 islands for 1 million years with a shift in immigration rate
#' ## at 0.195 Ma, and plot the species-through-time plot. Pool size 296.
#'
#' pars_before_shift <- c(0.079, 0.973, Inf, 0.136, 0.413)
#' pars_after_shift <- c(0.079, 0.973, Inf, 0.652, 0.413)
#' tshift <- 0.195
#' set.seed(1)
#' island_shift_replicates <- DAISIE_sim_cr_shift(
#'   time = 1,
#'   M = 296,
#'   pars = c(pars_before_shift, pars_after_shift),
#'   replicates = 2,
#'   shift_times = tshift,
#'   plot_sims = FALSE,
#'   verbose = FALSE
#' )
#' @export DAISIE_sim_cr
#' @export DAISIE_sim
DAISIE_sim_cr <- DAISIE_sim <- function(
  time,
  M,
  pars,
  replicates,
  divdepmodel = "CS",
  nonoceanic_pars = c(0, 0),
  num_guilds = NULL,
  prop_type2_pool = NA,
  replicates_apply_type2 = TRUE,
  sample_freq = 25,
  plot_sims = TRUE,
  hyper_pars = create_hyper_pars(d = 0, x = 0),
  area_pars = create_area_pars(
    max_area = 1,
    current_area = 1,
    proportional_peak_t = 0,
    total_island_age = 0,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0),
  cond = 0,
  verbose = TRUE,
  files_to_write = FALSE,
  ...
) {
  testit::assert(
    "length(pars) is not five and/or shift_times is not null, set
    five parameters with no shift_times or ten parameters with
    non-null shift_times",
    length(pars) == 5 || (length(pars) == 10 && !is.na(prop_type2_pool))
  )
  testit::assert(
    "2 type islands cannot have species on the island initially",
    is.na(prop_type2_pool) || !is.na(prop_type2_pool) && nonoceanic_pars[1] == 0
  )
  testit::assert(
    "prop_type2_pool should either be NA for no type 2 species or value between
    0 and 1",
    is.na(prop_type2_pool) || (prop_type2_pool >= 0 && prop_type2_pool <= 1)
  )
  testit::assert(are_hyper_pars(hyper_pars = hyper_pars))
  testit::assert(are_area_pars(area_pars = area_pars))

  total_time <- time

  if (divdepmodel == "IW") {
    island_replicates <- DAISIE_sim_cr_iw(
      total_time = total_time,
      M = M,
      pars = pars,
      replicates = replicates,
      nonoceanic_pars = nonoceanic_pars,
      sample_freq = sample_freq,
      hyper_pars = hyper_pars,
      area_pars = area_pars,
      cond = cond,
      verbose = verbose)
  }
  if (divdepmodel == "CS") {
    island_replicates <- DAISIE_sim_cr_cs(
      total_time = total_time,
      M = M,
      pars = pars,
      replicates = replicates,
      nonoceanic_pars = nonoceanic_pars,
      prop_type2_pool = prop_type2_pool,
      replicates_apply_type2 = replicates_apply_type2,
      sample_freq = sample_freq,
      hyper_pars = hyper_pars,
      area_pars = area_pars,
      cond = cond,
      verbose = verbose,
      files_to_write = files_to_write)
  }
  if (divdepmodel == "GW") {
    island_replicates <- DAISIE_sim_cr_gw(
      total_time = total_time,
      M = M,
      pars = pars,
      replicates = replicates,
      nonoceanic_pars = nonoceanic_pars,
      num_guilds = num_guilds,
      sample_freq = sample_freq,
      hyper_pars = hyper_pars,
      area_pars = area_pars,
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
