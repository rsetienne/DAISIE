#' @title DEPRECATED - Simulate islands with given parameters.
#' @description This function simulates islands with given cladogenesis, extinction, Kprime,
#' immigration and anagenesis parameters that shift at some time.
#' Runs the function with clade-specific
#' carrying capacity, where diversity-dependence operates only within single
#' clades, i.e. only among species originating from the same mainland
#' colonist.
#'
#' Returns R list object that contains the simulated islands.
#'
#' @details This function's use has been deprecated in favor
#' of \code{\link{DAISIE_sim_cr_shift}()}. Please use that
#' function instead.
#'
#' @param time Length of the simulation in time units. For example, if an
#' island is know to be 4 million years old, setting time = 4 will simulate
#' entire life span of the island; setting time = 2 will stop the simulation at
#' the mid-life of the island.
#' @param M The size of the mainland pool, i.e the number of species that can
#' potentially colonize the island
#' @param pars Contains the model parameters: \cr \cr \code{pars[1]}
#' corresponds to lambda^c (cladogenesis rate) before the shift \cr \code{pars[2]} corresponds
#' to mu (extinction rate) before the shift \cr \code{pars[3]} corresponds to K (clade-level
#' carrying capacity) before the shift. Set K=Inf for non-diversity dependence.\cr
#' \code{pars[4]} corresponds to gamma (immigration rate) before the shift \cr \code{pars[5]}
#' corresponds to lambda^a (anagenesis rate) before the shift \cr \code{pars[6]} corresponds to
#' lambda^c (cladogenesis rate) after the shift \cr \code{pars[7]}
#' corresponds to mu (extinction rate) after the shift \cr \code{pars[8]}
#' corresponds to K (clade-level carrying capacity) after the shift.  Set
#' K=Inf for non-diversity dependence. \cr \code{pars[9]} corresponds to gamma
#' (immigration rate)  after the shift \cr \code{pars[10]} corresponds to
#' lambda^a (anagenesis rate) after the shift \cr \code{pars[11]} corresponds to
#' the time of shift. This is defined as time before the end of the simulation. For example,
#' setting time = 4 and pars[11] = 1.5 will simulate with pars[1:5] from 4 to 1.5 and
#' with pars[6:10] from 1.5 to 0.
#' @param replicates Number of island replicates to be simulated.
#' @param sample_freq Specifies the number of units time should be divided by
#' for plotting purposes. Larger values will lead to plots with higher
#' resolution, but will also run slower.
#' @param plot_sims Default=TRUE plots species-through-time (STT) plots. It
#' detects how many types of species are present. If only one type of species
#' is present, STT is plotted for all species. If two types are present, three
#' plots are produced: STT for all, STT for type 1 and STT for type 2.
#' @param ddep = diversity-dependent model,mode of diversity-dependence \cr \cr
#' \code{ddep == 0} no diversity-dependence
#' \code{ddep == 1} linear dependence in speciation rate (anagenesis and cladogenesis)
#' \code{ddep == 11} linear dependence in speciation rate and immigration rate
#' \code{ddep == 3} linear dependence in extinction rate
#' @param verbose \code{Default=TRUE} Give intermediate output, also if
#'   everything goes OK.
#' @param ...  Any arguments to pass on to plotting functions.
#' @return Each simulated dataset is an element of the list, which can be
#' called using [[x]]. For example if the object is called island_replicates,
#' the 1st replicate can be called using island_replicates[[1]] Each of the
#' island replicates is a list in itself. The first (e.g.
#' island_replicates[[x]][[1]]) element of that list has the following
#' components: \cr \code{$island_age} - the island age \cr Then, we have:\cr \code{$not_present}
#' - the number of mainland lineages that are not present on the island \cr
#' \code{$stt_all} - STT table for all species on the island (nI - number of
#' non-endemic species; nA - number of anagenetic species, nC - number of
#' cladogenetic species, present - number of independent colonisations present)\cr
#' The subsequent elements of the list each contain information on a single
#' colonist lineage on the island and has 4 components:\cr
#' \code{$branching_times} - island age and stem age of the population/species
#' in the case of Non-endemic, Non-endemic_MaxAge and Endemic anagenetic
#' species. For cladogenetic species these should be island age and branching
#' times of the radiation including the stem age of the radiation.\cr
#' \code{$stac} - the status of the colonist \cr * Non_endemic_MaxAge: 1 \cr *
#' Endemic: 2 \cr * Endemic&Non_Endemic: 3 \cr * Non_endemic: 4 \cr
#' \code{$missing_species} - number of island species that were not sampled for
#' particular clade (only applicable for endemic clades) \cr \code{$type_1or2}
#' - whether the colonist belongs to type 1 or type 2 \cr
#' @author Luis Valente, Albert Phillimore, and Torsten Hauffe
#' @seealso \code{\link{DAISIE_plot_sims}()}
#' @references Hauffe, T., D. Delicado, R.S. Etienne and L. Valente (2020).
#' Lake expansion increases equilibrium diversity via the target effect of
#' island biogeography
#' @noRd
DAISIE_SR_sim <- function(time,
                          M,
                          pars,
                          replicates,
                          sample_freq = 25,
                          plot_sims = TRUE,
                          ddep = 11,
                          verbose = FALSE,
                          ...)
{
  message("This is a deprecated function to replicate the Hauffe et al. 2020
           paper.\nFor a general purpose split rates model, please use
          DAISIE_sim_cr_shift().")
  if (length(pars) != 11)
  {
    stop("Shift in rates requires 11 parameters")
  }
  island_replicates = list()
  for (rep in 1:replicates)
  {
    island_replicates[[rep]] = list()
    full_list = list()
    parstmp <- c(pars[1:10], time - pars[11])
    for (m_spec in 1:M)
    {
      full_list[[m_spec]] = DAISIE_SR_sim_core(time = time, mainland_n = 1, parstmp)
    }
    island_replicates[[rep]] = full_list
    if (verbose)
      message("Island replicate ", rep)
  }
  island_replicates = DAISIE_format_CS(
    island_replicates = island_replicates,
    time = time,
    M = M,
    sample_freq = sample_freq,
    verbose = verbose
  )
  if (plot_sims == TRUE)
  {
    DAISIE_plot_sims(island_replicates)
  }
  return(island_replicates)
}
