#' @title Simulate (non-)oceanic islands with given parameters under a
#'   time-dependent regime
#'
#' @description
#' This function simulates islands with given cladogenesis,
#' extinction, Kprime, immigration and anagenesis parameters, all of which
#' modelled as time-dependent parameters.
#'
#' Time dependency aims to capture the
#' effect of area changes islands undego from their emergence until subsidence.
#' Thus, oceanic, volcanic island ontogeny scenarios can be modelled (by a
#' beta function), as well as the effect of sea level fluctuations (modelled
#' through a sine function). See paramter entry \code{area_pars} for details.
#' Both island ontogeny and sea level fluctuations are allowed to operate
#' simultaneuosly.
#'
#' This function also allows for the simulation of
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
#' @author Luis Valente and Albert Phillimore
#' @seealso \code{\link{DAISIE_plot_sims}()} for plotting STT of simulation
#' outputs.
#' @family simulation models
#' @references Valente, L.M., A.B. Phillimore and R.S. Etienne (2015).
#' Equilibrium and non-equilibrium dynamics simultaneously operate in the
#' Galapagos islands. Ecology Letters 18: 844-852.
#'
#' Valente, L.M., Etienne, R.S. and Phillimore, A.B. (2014). The effects of
#' island ontogeny on species diversity and phylogeny.
#' Proceedings of the Royal Society B: Biological Sciences 281(1784),
#' p.20133227.
#' @keywords models
#' @export
DAISIE_sim_time_dependent <- function(
  time,
  M,
  pars,
  replicates,
  area_pars,
  hyper_pars,
  divdepmodel = "CS",
  nonoceanic_pars = c(0, 0),
  num_guilds = NULL,
  sample_freq = 25,
  plot_sims = TRUE,
  island_ontogeny = "const",
  sea_level = "const",
  extcutoff = 1000,
  verbose = TRUE,
  ...
) {
  testit::assert(
    "island_ontogeny is not valid input. Specify 'const',\n
    or  'beta'", is_island_ontogeny_input(island_ontogeny)
  )
  testit::assert(
    "sea_level is not valid input. Specify 'const, \n or 'sine'",
    is_sea_level_input(sea_level)
  )

  testit::assert(
    "length(pars) is not five",
    length(pars) == 5
  )
  testit::assert(
    "one hyper parameter must be non-zero for time-dependency",
    hyper_pars$d != 0 || hyper_pars$x != 0)

  totaltime <- time

  testit::assert(are_hyper_pars(hyper_pars = hyper_pars))
  testit::assert(are_area_pars(area_pars = area_pars))
  testit::assert(totaltime <= area_pars$total_island_age)

  island_replicates <- list()
  island_ontogeny <- translate_island_ontogeny(island_ontogeny)
  sea_level <- translate_sea_level(sea_level)

  if (island_ontogeny == 1) {
    peak <- calc_peak(totaltime,
                      area_pars)
  } else {
    peak <- NULL
  }

  Amax <- get_global_max_area(
    totaltime = totaltime,
    area_pars = area_pars,
    peak = peak,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level
  )
  Amin <- get_global_min_area(
    totaltime = totaltime,
    area_pars = area_pars,
    peak = peak,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level
  )

  testit::assert(is.numeric(Amax))
  testit::assert(is.finite(Amax))
  testit::assert(is.numeric(Amin))
  testit::assert(is.finite(Amin))

  #### IW ####
  if (divdepmodel == "IW") {
    for (rep in 1:replicates) {
      island_replicates[[rep]] <- DAISIE_sim_core_time_dependent(
        time = totaltime,
        mainland_n = M,
        pars = pars,
        nonoceanic_pars = nonoceanic_pars,
        island_ontogeny = island_ontogeny,
        sea_level = sea_level,
        hyper_pars = hyper_pars,
        area_pars = area_pars,
        peak = peak,
        Amax = Amax,
        Amin = Amin,
        extcutoff = extcutoff
      )
      if (verbose == TRUE) {
        print(paste("Island replicate ", rep, sep = ""))
      }
    }
    island_replicates <- DAISIE_format_IW(island_replicates = island_replicates,
                                          time = totaltime,
                                          M = M,
                                          sample_freq = sample_freq,
                                          verbose = verbose)
  }

  #### CS ####
  if (divdepmodel == "CS") {
  for (rep in 1:replicates) {
    island_replicates[[rep]] <- list()
    full_list <- list()
    for (m_spec in 1:M) {
      full_list[[m_spec]] <- DAISIE_sim_core_time_dependent(
        time = totaltime,
        mainland_n = 1,
        pars = pars,
        nonoceanic_pars = nonoceanic_pars,
        island_ontogeny = island_ontogeny,
        sea_level = sea_level,
        hyper_pars = hyper_pars,
        area_pars = area_pars,
        peak = peak,
        Amax = Amax,
        Amin = Amin,
        extcutoff = extcutoff
      )
    }
    island_replicates[[rep]] <- full_list
    if (verbose == TRUE) {
      print(paste("Island replicate ", rep, sep = ""))
    }
  }
  island_replicates <- DAISIE_format_CS(
    island_replicates = island_replicates,
    time = totaltime,
    M = M,
    sample_freq = sample_freq,
    verbose = verbose
  )
  }

  #### GW ####
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
        full_list[[m_spec]]  <- DAISIE_sim_core_time_dependent(
          time = totaltime,
          mainland_n = guild_size,
          pars = pars,
          nonoceanic_pars = nonoceanic_pars,
          island_ontogeny = island_ontogeny,
          sea_level = sea_level,
          hyper_pars = hyper_pars,
          area_pars = area_pars,
          peak = peak,
          Amax = Amax,
          Amin = Amin,
          extcutoff = extcutoff
        )
      }
      island_replicates[[rep]] <- full_list
      if (verbose == TRUE) {
        print(paste("Island replicate ", rep, sep = ""))
      }
    }
    island_replicates <- DAISIE_format_GW(island_replicates = island_replicates,
                                          time = totaltime,
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
