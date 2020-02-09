#' Defailt parameter documentation
#'
#' @param time Numeric defining the length of the simulation in time units.
#' For examples, if an island is known to be 4 million years old, setting
#' time = 4 will simulate the entire life span of the island; setting time = 2
#' will stop the simulation at the mid-life of the island.
#' @param M Numeric defining the size of mainland pool, i.e. the number of
#' species that can potentially colonize the island.
#' @param pars A numeric vector containing the model parameters:
#' \itemize {
#'   \item{pars[1]: lambda^c (cladogenesis rate)}
#'   \item{pars[2]: mu (extinction rate)}
#'   \item{pars[3]: K (carrying capacity), set K=Inf for diversity
#'   independence.}
#'   \item{pars[4]: gamma (immigration rate)}
#'   \item{pars[5]: lambda^a (anagenesis rate)}
#'   \item{pars[6]: lambda^c (cladogenesis rate) for either type 2 species
#'   or rate set 2 in rate shift model}
#'   \item{pars[7]: mu (extinction rate) for either type 2 species or rate
#'   set 2 in rate shift model}
#'   \item{pars[8]: K (carrying capacity) for either type 2 species or rate
#'   set 2 in rate shift model, set K=Inf for diversity independence.}
#'   \item{pars[9]: gamma (immigration rate) for either type 2 species or
#'   rate set 2 in rate shift model}
#'   \item{pars[10]: lambda^a (anagenesis rate) for either type 2 species
#'   or rate set 2 in rate shift model}
#' }
#' \cr The elements 6:10 are required only when type 2 species are included
#' or in the rate shift model.
#' @param replicates Number of island replicates to be simulated.
#' @param divdepmodel Option divdepmodel = 'CS' runs a model with clade-specific
#' carrying capacity, where diversity-dependence operates only within single
#' clades, i.e. only among species originating from the same mainland colonist.
#' Option divdepmodel = 'IW' runs a model with island-wide carrying capacity,
#' where diversity-dependence operates within and among clades. Option
#' divdepmodel = 'GW' runs a model with diversity-dependence operates within
#' a guild.
#' @param nonoceanic_pars A vector of length two with:.
#' #' \itemize{
#'   \item{[1]: the probability of sampling a species from the mainland}
#'   \item{[2]: the probability of the species sampled from the mainland
#'   being nonendemic}
#'   }
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
#' @param sample_freq Numeric specifing the number of units times should be
#' divided by for plotting purposes. Larger values will lead to plots with
#' higher resolution, but will also run slower.
#' @param plot_sims Default = TRUE plots species-through-time (STT) plots. It
#' detects how many types of species are present. If only one type of species
#' is present, STT is plotted for all species. If two types are present, three
#' plots are produced: STT for all, STT for type 1 and STT for type 2.
#' @param verbose Logical, \code{Default = TRUE} Give intermediate output,
#' also if everything goes ok.
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
#'   \item{[7]: angle of the slope of the island}
#' }
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
#' @param ext_pars A numeric vector:
#' \itemize{
#'   \item{[1]: minimum extinction when area is at peak}
#'   \item{[2]: extinction rate when current area is 0.10 of maximum area}
#' }
#' @param island_ontogeny In \code{\link{DAISIE_sim_constant_rate}} a string
#' describing the type of island ontogeny. Can be \code{"const"}, \code{"beta"}
#' for a beta function describing area through time.
#' \cr In \code{\link{DAISIE_sim_core_constant_rate}} a numeric describing the
#' type of island ontogeny. Can be \code{0} for constant, \code{1} for a beta
#' function describing area through time.
#' @param sea_level In In \code{\link{DAISIE_sim_constant_rate}} a string
#' describing the type of sea level. Can be \code{"const"} or \code{"sine"}
#' for a sine function describing area through time.
#' \cr In \code{\link{DAISIE_sim_core_constant_rate}} a numeric describing the
#' type of sea level. Can be \code{0} for constant, \code{1} for a sine
#' function describing area through time.
#' @param extcutoff the maximum extinction rate.
#' @param shift_times a numeric vector specifying when the rate shifts occur
#' before the present.
#' @param mainland_n A numeric stating the number of mainland species, that
#' is the number of species that can potentially colonize the island.
#' If using a clade-specific diversity dependence, this value is set to 1.
#' If using an island-wide diversity dependence, this value is set to the
#' number of mainland species.
#' @param island_replicates List output from
#' \code{\link{DAISIE_sim_core_constant_rate}},
#' \code{\link{DAISIE_sim_core_time_dependent}}, or
#' \code{\link{DAISIE_sim_core_constant_rate_shift}} functions.
#' @param island_spec matrix with current state of simulation.
#' @param stt_table matrix with number of species at each time step.
#' @param ... Any arguments to pass on to plotting functions.
#'
#'
#' @return Nothing
#'
default_params_doc <- function(
  time,
  M,
  pars,
  replicates,
  divdepmodel,
  nonoceanic_pars,
  num_guilds,
  prop_type2_pool,
  replicates_apply_type2,
  sample_freq,
  plot_sims,
  verbose,
  area_pars,
  hyper_pars,
  dist_pars,
  ext_pars,
  island_ontogeny,
  sea_level,
  extcufoff,
  shift_times,
  mainland_n,
  island_replicates,
  island_spec,
  ...,
) {
  # Nothing
}

