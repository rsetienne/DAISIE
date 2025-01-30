#' @title Simulate islands with given trait-dependent parameters.
#' @description This function simulates islands with given cladogenesis,
#' extinction, K, immigration and anagenesis parameters for binary trait states.
#'
#' Returns R list object that contains the simulated islands
#'
#' @inheritParams default_params_doc
#'
#' @return
#' A list. The highest level of the least corresponds to each individual
#' replicate. The first element of each replicate is composed of island
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
#' \item{\code{$missing_species}: number of island species that were
#' not sampled for particular clade (only applicable for endemic clades)}
#' \item{\code{$type_1or2}: whether the colonist belongs to type 1 or type 2}
#' }
#' @author Luis Valente and Albert Phillimore
#' @seealso \code{\link{DAISIE_format_CS}} \code{\link{DAISIE_plot_sims}}
#' @references Valente, L.M., A.B. Phillimore and R.S. Etienne (2015).
#' Equilibrium and non-equilibrium dynamics simultaneously operate in the
#' Galapagos islands. Ecology Letters 18: 844-852.
#' Hauffe, T., D. Delicado, R.S. Etienne and L. Valente (submitted).
#' Lake expansion increases equilibrium diversity via the target effect of
#' island biogeography.
#' @keywords models
#' @export
DAISIE_sim_trait_dep <- function(
  time,
  M,
  pars,
  replicates,
  divdepmodel = "CS",
  sample_freq = 25,
  plot_sims = TRUE,
  island_ontogeny = "const",
  sea_level = "const",
  hyper_pars = create_hyper_pars(d = 0, x = 0),
  area_pars = DAISIE::create_area_pars(
    max_area = 1,
    current_area = 1,
    proportional_peak_t = 0,
    total_island_age = 0,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0),
  extcutoff = 1000,
  cond = 0,
  verbose = TRUE,
  trait_pars = NULL,
  ...
) {
  total_time <- time
  island_replicates <- list()
  island_ontogeny <- translate_island_ontogeny(island_ontogeny)
  sea_level <- translate_sea_level(sea_level)

  #### IW ####
  if (divdepmodel == "IW") {
    for (rep in 1:replicates) {
      island_replicates[[rep]] <- DAISIE_sim_core_trait_dep(
        time = total_time,
        mainland_n = M,
        pars = pars,
        island_ontogeny = island_ontogeny,
        sea_level = sea_level,
        hyper_pars = hyper_pars,
        area_pars = area_pars,
        extcutoff = extcutoff,
        trait_pars = trait_pars
      )
      if (verbose == TRUE) {
        message("Island replicate ", rep)
      }
    }
    island_replicates <- DAISIE_format_IW(
      island_replicates = island_replicates,
      time = total_time,
      M = M,
      sample_freq = sample_freq,
      verbose = verbose,
      trait_pars = trait_pars)
  }

  #### CS ####
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
        if(M == 0 || is.null(trait_pars)){
          stop("One state exists on mainland, should use constant rate DAISIE.")
        }else{
          for (m_spec in 1:M) {
            ### M1 = 1, M2 = 0
            trait_pars_onecolonize <- create_trait_pars(
              trans_rate = trait_pars$trans_rate,
              immig_rate2 = trait_pars$immig_rate2,
              ext_rate2 = trait_pars$ext_rate2,
              ana_rate2 = trait_pars$ana_rate2,
              clado_rate2 = trait_pars$clado_rate2,
              trans_rate2 = trait_pars$trans_rate2,
              M2 = 0)
            full_list[[m_spec]] <- DAISIE_sim_core_trait_dep(
              time = total_time,
              mainland_n = 1,
              pars = pars,
              island_ontogeny = island_ontogeny,
              sea_level = sea_level,
              hyper_pars = hyper_pars,
              area_pars = area_pars,
              extcutoff = extcutoff,
              trait_pars = trait_pars_onecolonize
            )
          }
          for(m_spec in (M + 1):(M + trait_pars$M2)) {
            ### M1 = 0, M2 = 1
            trait_pars_onecolonize <- create_trait_pars(
              trans_rate = trait_pars$trans_rate,
              immig_rate2 = trait_pars$immig_rate2,
              ext_rate2 = trait_pars$ext_rate2,
              ana_rate2 = trait_pars$ana_rate2,
              clado_rate2 = trait_pars$clado_rate2,
              trans_rate2 = trait_pars$trans_rate2,
              M2 = 1)
            full_list[[m_spec]] <- DAISIE_sim_core_trait_dep(
              time = total_time,
              mainland_n = 0,
              pars = pars,
              island_ontogeny = island_ontogeny,
              sea_level = sea_level,
              hyper_pars = hyper_pars,
              area_pars = area_pars,
              extcutoff = extcutoff,
              trait_pars = trait_pars_onecolonize
            )
          }
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
      verbose = verbose,
      trait_pars = trait_pars
    )
  }

  if (plot_sims == TRUE) {
    DAISIE_plot_sims(
      island_replicates = island_replicates,
      sample_freq = sample_freq,
      trait_pars = trait_pars
    )
  }
  return(island_replicates)
}
