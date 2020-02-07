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
#' @param time Length of the simulation in time units. For examples, if an
#' island is known to be 4 million years old, setting time = 4 will simulate
#' the entire life span of the island; setting time = 2 will stop the
#' simulation at the mid-life of the island.
#' @param M The size of mainland pool, i.e. the number of species that can
#' potentially colonize the island.
#' @param pars Contains the model parameters: \cr \cr \code{pars[1]}
#' corresponds to lambda^c (cladogenesis rate) \cr \code{pars[2]} corresponds
#' to mu (extinction rate) \cr \code{pars[3]} corresponds to K (clade-level
#' carrying capacity). Set K=Inf for non-diversity dependence.\cr
#' \code{pars[4]} corresponds to gamma (immigration rate) \cr \code{pars[5]}
#' corresponds to lambda^a (anagenesis rate) \cr \code{pars[6]} corresponds to
#' lambda^c (cladogenesis rate) for type 2 species \cr \code{pars[7]}
#' corresponds to mu (extinction rate) for type 2 species\cr \code{pars[8]}
#' corresponds to K (clade-level carrying capacity) for type 2 species.  Set
#' K=Inf for non-diversity dependence.\cr \code{pars[9]} corresponds to gamma
#' (immigration rate) for type 2 species\cr \code{pars[10]} corresponds to
#' lambda^a (anagenesis rate) for type 2 species\cr The elements 6:10 are
#' optional and are required only when type 2 species are included.
#' @param replicates Number of island replicates to be simulated. When
#' pars_shift = TRUE, the pars[6:10] corresponds to lambda^c (cladogenesis
#' rate, mu (extinction rate), K (carrying capacity), gamma (immigration rate),
#' and lambda^a (anagensis rate) respectively. \code{pars[11]} is defined as
#' the time before the end of the simulation. For example, setting time = 4
#' and \code{pars[11]} = 1.5 will simulate with \code{pars[1:5]} from 4 to 1.5
#' and with \code{pars[6:10]} from 1.5 to 0.
#' @param divdepmodel Option divdepmodel = 'CS' runs a model with clade-specific
#' carrying capacity, where diversity-dependence operates only within single
#' clades, i.e. only among species originating from the same mainland colonist.
#' Option divdepmodel = 'IW' runs a model with island-wide carrying capacity,
#' where diversity-dependence operates within and among clades. Option
#' divdepmodel = 'GW' runs a model with diversity-dependence operates within
#' a guild.
#' @param nonoceanic_pars A vector of length three with: the island area as a
#' proportion of the mainland, the probability of native species being
#' nonendemic and the size of the mainland pool.
#' @param num_guilds The number of guilds on the mainland. The number of
#' mainland species is divided by the number of guilds when \code{divdepmodel =
#' "GW"}
#' @param prop_type2_pool Fraction of mainland species that belongs to the
#' second subset of species (type 2). Applies only when two types of species
#' are simulated (length(pars) = 10).
#' @param replicates_apply_type2 Applies only when two types of species are
#' being simulated. Default replicates_apply_type2 = TRUE runs simulations
#' until the number of islands where a type 2 species has colonised is equal to
#' the specified number of replicates. This is recommended if prop_type2_pool
#' is small of if the rate of immigration of type two species (pars[9]) is low,
#' meaning that more replicates are needed to achieved an adequate sample size
#' of islands with type 2 species. Setting replicates_apply_type2 = FALSE
#' simulates islands up to the specified number of replicates regardless of
#' whether type 2 species have colonised or not.
#' @param sample_freq Specifies the number of units times should be divided by
#' for plotting purposes. Larger values will lead to plots with higher
#' resolution, but will also run slower.
#' @param plot_sims Default = TRUE plots species-through-time (STT) plots. It
#' detects how many types of species are present. If only one type of species
#' is present, STT is plotted for all species. If two types are present, three
#' plots are produced: STT for all, STT for type 1 and STT for type 2.
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
#' @param verbose \code{Default = TRUE} Give intermediate output, also if
#' everything goes ok.
#' @param ... Any arguments to pass on to plotting functions.
#' @param hyper_pars A numeric vector for hyperparameters for the rate
#' calculations:
#' \itemize{
#' \item{[1]: is d_0 the scaling parameter for exponent for calculating
#' cladogenesis rate}
#' \item{[2]: is x the exponent for calculating extinction rate}
#' \item{[3]: is alpha, the exponent for calculating the immigration rate}
#' \item{[4]: is beta the exponent for calculating the anagenesis rate.}
#' }
#' @param dist_pars a numeric for the distance from the mainland.
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
