#' @title Simulate islands with given parameters.
#' @description This function simulates islands with given cladogenesis,
#' extinction, Kprime, immigration and anagenesis parameters. If a single
#' parameter set is provided (5 parameters) it simulates islands where all
#' species have the same macro-evolutionary process. If two paramater sets
#' (10 parameters) are provided, it simulates islands where two different
#' macro-evolutionary processes operate, one applying to type 1 species and
#' other to type 2 species. If two parameter sets and a time shift (11
#' parameters) are provided, it simulates islands where at the given time
#' a shift between the parameter sets will occur.
#'
#' Returns R list object that contains the simulated islands
#'
#' @inheritParams default_params_doc
#'
#' @return Each simulated dataset is an element of the list, which can be
#' called using [[x]]. For example if the object is called island_replicates,
#' the last replicates is a list in itself. The first (e.g.
#' island_replicates[[x]][[1]]) element of that list has the following
#' components: \cr \code{$island_age} - the island age \cr Then, depending on
#' whether a distinction between types is made, we have:\cr \code{$not_present}
#' - the number of mainland lineages that are not present on the island \cr
#' or:\cr \code{$not_present_type1} - the number of mainland lineages of type 1
#' that are not present on the island \cr \code{$not_present_type2} - the
#' number of mainland lineages of type 2 that are not present on the island \cr
#' \code{$stt_all} - STT table for all species on the island (nI - number of
#' non-endemic species; nA - number of anagenetic species, nC - number of
#' cladogenetic species, present - number of independent colonisations present
#' )\cr \code{$stt_stt_type1} - STT table for type 1 species on the island -
#' only if 2 types of species were simulated (nI - number of non-endemic
#' species; nA - number of anagenetic species, nC - number of cladogenetic
#' species, present - number of independent colonisations present )\cr
#' \code{$stt_stt_type2} - STT table for type 2 species on the island - only if
#' 2 types of species were simulated (nI - number of non-endemic species; nA -
#' number of anagenetic species, nC - number of cladogenetic species, present -
#' number of independent colonisations present )\cr \code{$brts_table} - Only
#' for simulations under 'IW'. Table containing information on order of events
#' in the data, for use in maximum likelihood optimization.)\cr
#'
#' The subsequent elements of the list each contain information on a single
#' colonist lineage on the island and has 4 components:\cr
#' \code{$branching_times} - island age and stem age of the population/species
#' in the case of Non-endemic, Non-endemic_MaxAge and Endemic anagenetic
#' species. For cladogenetic species these should be island age and branching
#' times of the radiation including the stem age of the radiation.\cr
#' \code{$stac} - the status of the colonist \cr * Non_endemic_MaxAge: 1 \cr *
#' ndemic: 2 \cr * Endemic&Non_Endemic: 3 \cr * Non_endemic: 4 \cr
#' \code{$missing_species} - number of island species that were not sampled for
#' particular clade (only applicable for endemic clades) \cr \code{$type_1or2}
#' - whether the colonist belongs to type 1 or type 2 \cr
#' @author Luis Valente and Albert Phillimore
#' @seealso \code{\link{DAISIE_format_CS}} \code{\link{DAISIE_plot_sims}}
#' @references Valente, L.M., A.B. Phillimore and R.S. Etienne (2015).
#' Equilibrium and non-equilibrium dynamics simultaneously operate in the
#' Galapagos islands. Ecology Letters 18: 844-852.
#' Hauffe, T., D. Delicado, R.S. Etienne and L. Valente (submitted).
#' Lake expansion increases equilibrium diversity via the target effect of
#' island biogeography.
#' @keywords models
#' @examples
#'
#' cat("
#' ## Simulate 40 islands for 4 million years, where all species have equal
#' ## rates, and plot the species-through-time plot. Pool size 1000.
#'
#' pars_equal = c(2.550687345,2.683454548,Inf,0.00933207,1.010073119)
#' island_replicates_equal = DAISIE_sim_constant_rate(
#'    time = 4,
#'    M = 1000,
#'    pars = pars_equal,
#'    replicates = 40
#'    )
#'
#' ## Simulate 15 islands for 4 million years with two types of species (type1
#' ## and type 2), and plot the species-through-time plot. Pool size 1000.
#' ## Fraction of type 2 species in source pool is 0.163. Function will
#' ## simulate until number of islands where type 2 species has colonised is
#' ## equal to number specified in replicates.
#'
#' pars_type1 = c(0.195442017,0.087959583,Inf,0.002247364,0.873605049)
#' pars_type2 = c(3755.202241,8.909285094,14.99999923,0.002247364,0.873605049)
#' island_replicates_2types = DAISIE_sim_constant_rate(
#'    time = 4,
#'    M = 1000,
#'    pars = c(pars_type1,pars_type2),
#'    replicates = 15,
#'    prop_type2_pool = 0.163
#'    )
#' ## Simulate a non-oceanic island for 4 million years, and plot the
#' ## species-through-time plot. Pool size 1000. Island area as a proportion
#' ## of mainland is 0.1, proportion of native species is 0.9.
#'  pars = c(2.550687345,2.683454548,Inf,0.00933207,1.010073119)
#'  island_replicates = DAISIE_sim_constant_rate(
#'    time = 4,
#'    M = 1000,
#'    pars = pars,
#'    replicates = 40,
#'    nonoceanic_pars = c(0.1, 0.9)
#'  )
#' ## Simulate 15 islands for 4 million years with a shift in immigration rate
#' ## at 0.195 Ma, and plot the species-through-time plot. Pool size 296.
#'
#' pars_before_shift = c(0.079, 0.973, Inf, 0.136, 0.413)
#' pars_after_shift = c(0.079, 0.973, Inf, 0.652, 0.413)
#' tshift = 0.195
#' island_shift_replicates = DAISIE_sim_constant_rate_shift(
#'    time = 4,
#'    M = 296,
#'    pars = c(pars_before_shift, pars_after_shift),
#'    replicates = 15,
#'    shift_times = tshift
#'  )
#' ")
#'
#' @export
DAISIE_sim_constant_rate <- function(
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
  hyper_pars = NULL,
  area_pars = NULL,
  dist_pars = NULL,
  verbose = TRUE,
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

  totaltime <- time
  island_replicates <- list()
  if (divdepmodel == "IW") {
    for (rep in 1:replicates) {
      island_replicates[[rep]] <- DAISIE_sim_core_constant_rate(
        time = totaltime,
        mainland_n = M,
        pars = pars,
        nonoceanic_pars = nonoceanic_pars,
        hyper_pars = hyper_pars,
        area_pars = area_pars,
        dist_pars = dist_pars
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
  if (divdepmodel == "CS") {
    if (length(pars) == 5) {
      for (rep in 1:replicates) {
        island_replicates[[rep]] <- list()
        full_list <- list()
        for (m_spec in 1:M) {
          full_list[[m_spec]] <- DAISIE_sim_core_constant_rate(
            time = totaltime,
            mainland_n = 1,
            pars = pars,
            nonoceanic_pars = nonoceanic_pars,
            hyper_pars = hyper_pars,
            area_pars = area_pars,
            dist_pars = dist_pars
          )
        }
        island_replicates[[rep]] <- full_list
        if (verbose == TRUE) {
          print(paste("Island replicate ", rep, sep = ""))
        }
      }
    } else if (length(pars) == 10) {
      if (replicates_apply_type2 == TRUE) {
        island_replicates <- DAISIE_sim_min_type2(
          time = totaltime,
          M = M,
          pars = pars,
          replicates = replicates,
          prop_type2_pool = prop_type2_pool,
          verbose = verbose)
      } else {
        for (rep in 1:replicates) {
          pool2 <- DDD::roundn(M * prop_type2_pool)
          pool1 <- M - pool2
          lac_1 <- pars[1]
          mu_1 <- pars[2]
          K_1 <- pars[3]
          gam_1 <- pars[4]
          laa_1 <- pars[5]
          lac_2 <- pars[6]
          mu_2 <- pars[7]
          K_2 <- pars[8]
          gam_2 <- pars[9]
          laa_2 <- pars[10]
          full_list <- list()
          #### species of pool1
          for (m_spec in 1:pool1) {
            full_list[[m_spec]] <- DAISIE_sim_core_constant_rate(
              time = totaltime,
              mainland_n = 1,
              pars = c(lac_1,
                       mu_1,
                       K_1,
                       gam_1,
                       laa_1)
            )
            full_list[[m_spec]]$type1or2  <- 1
          }
          #### species of pool2
          for (m_spec in (pool1 + 1):(pool1 + pool2)) {
            full_list[[m_spec]] <- DAISIE_sim_core_constant_rate(
              time = totaltime,
              mainland_n = 1,
              pars = c(lac_2,
                       mu_2,
                       K_2,
                       gam_2,
                       laa_2)
            )
            full_list[[m_spec]]$type1or2 <- 2
          }
          island_replicates[[rep]] <- full_list
          if (verbose == TRUE) {
            print(paste("Island replicate ", rep, sep = ""))
          }
        }
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
        full_list[[m_spec]]  <- DAISIE_sim_core_constant_rate(
          time = totaltime,
          mainland_n = guild_size,
          pars = pars,
          nonoceanic_pars = nonoceanic_pars,
          hyper_pars = hyper_pars,
          area_pars = area_pars,
          dist_pars = dist_pars
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
