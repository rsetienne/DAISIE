#' Default parameter documentation
#'
#' @param time Numeric defining the length of the simulation in time units.
#'   For example, if an island is known to be 4 million years old, setting
#'   time = 4 will simulate the entire life span of the island; setting time = 2
#'   will stop the simulation at the mid-life of the island.
#' @param M Numeric defining the size of mainland pool, i.e. the number of
#'   species that can potentially colonize the island.
#' @param pars A numeric vector containing the model parameters:
#'   \itemize{
#'     \item{\code{pars[1]}: lambda^c (cladogenesis rate)}
#'     \item{\code{pars[2]}: mu (extinction rate)}
#'     \item{\code{pars[3]}: K (carrying capacity), set K=Inf for diversity
#'     independence.}
#'     \item{\code{pars[4]}: gamma (immigration rate)}
#'     \item{\code{pars[5]}: lambda^a (anagenesis rate)}
#'     \item{\code{pars[6]}: lambda^c (cladogenesis rate) for either type 2 species
#'     or rate set 2 in rate shift model}
#'     \item{\code{pars[7]}: mu (extinction rate) for either type 2 species or rate
#'     set 2 in rate shift model}
#'     \item{\code{pars[8]}: K (carrying capacity) for either type 2 species or rate
#'     set 2 in rate shift model, set K=Inf for diversity independence.}
#'     \item{\code{pars[9]}: gamma (immigration rate) for either type 2 species
#'     or rate set 2 in rate shift model}
#'     \item{\code{pars[10]}: lambda^a (anagenesis rate) for either type 2
#'     species or rate set 2 in rate shift model}
#'   }
#'   Elements 6:10 are required only when type 2 species are included
#'   or in the rate shift model. For \code{\link{DAISIE_sim_relaxed_rate}()}
#'   \code{pars[6]} is the standard deviation of the gamma distribution for the
#'   relaxed parameter and the parameter chosen by the \code{relaxed_par}
#'   argument is the mean of the gamma distribution for the relaxed parameter.
#' @param replicates Integer specifying number of island replicates to be
#'   simulated.
#' @param divdepmodel Option divdepmodel = 'CS' runs a model with clade-specific
#'   carrying capacity, where diversity-dependence operates only within single
#'   clades, i.e. only among species originating from the same mainland
#'   colonist. Option divdepmodel = 'IW' runs a model with island-wide
#'   carrying capacity, where diversity-dependence operates within and among
#'   clades. Option divdepmodel = 'GW' runs a model with diversity-dependence
#'   operates within a guild.
#' @param nonoceanic_pars A vector of length two with:
#'   \itemize{
#'     \item{[1]: the probability of sampling a species from the mainland}
#'     \item{[2]: the probability of the species sampled from the mainland
#'     being nonendemic}
#'   }
#' @param num_guilds The number of guilds on the mainland. The number of
#'   mainland species is divided by the number of guilds when
#'   \code{divdepmodel = "GW"}
#' @param prop_type2_pool Fraction of mainland species that belongs to the
#'   second subset of species (type 2). Applies only when two types of species
#'   are simulated (length(pars) = 10). For \code{\link{DAISIE_dataprep}()}
#'   applies only if number_clade_types = 2.  In \code{\link{DAISIE_dataprep}()}
#'   the default \code{"proportional"} sets the fraction to be proportional to
#'   the number of clades of distinct macroevolutionary process that have
#'   colonised the island.
#' @param replicates_apply_type2 Applies only when two types of species are
#'   being simulated. Default replicates_apply_type2 = TRUE runs simulations
#'   until the number of islands where a type 2 species has colonised is equal
#'   to the specified number of replicates. This is recommended if
#'   \code{prop_type2_pool} is small or if the rate of immigration of type two
#'   species (\code{pars[9]}) is low, meaning that more replicates are needed to
#'   achieved an adequate sample size of islands with type 2 species. Setting
#'   \code{replicates_apply_type2 = FALSE} simulates islands up to the
#'   specified number of replicates regardless of whether type 2 species have
#'   colonised or not.
#' @param sample_freq Numeric specifing the number of units times should be
#'   divided by for plotting purposes. Larger values will lead to plots with
#'   higher resolution, but will also run slower.
#' @param plot_sims \code{Default = TRUE} plots species-through-time (STT)
#'   plots. It detects how many types of species are present. If only one type
#'   of species is present, STT is plotted for all species. If two types are
#'   present, three plots are produced: STT for all, STT for type 1 and STT for
#'   type 2.
#' @param verbose In simulation and dataprep functions a logical,
#'   \code{Default = TRUE} gives intermediate output should be printed.
#'   For ML functions a numeric determining if intermediate output should be
#'   printed, \code{Default = 0} does not print, \code{verbose = 1} prints
#'   intermediate output of the parameters and loglikelihood, \code{verbose = 2}
#'   means also intermediate progress during loglikelihood computation is shown.
#' @param area_pars A named list containing area and sea level parameters as
#'   created by \code{\link{create_area_pars}()}:
#'   \itemize{
#'     \item{[1]: maximum area}
#'     \item{[2]: current area}
#'     \item{[3]: value from 0 to 1 indicating where in the island's history the
#'     peak area is achieved}
#'     \item{[4]: total island age}
#'     \item{[5]: amplitude of area fluctuation from sea level}
#'     \item{[6]: frequency of sine wave of area change from sea level}
#'     \item{[7]: angle of the slope of the island}
#'   }
#' @param hyper_pars A named list of numeric hyperparameters for the rate
#'   calculations as returned by \code{\link{create_hyper_pars}()}:
#'   \itemize{
#'     \item{[1]: is d the scaling parameter for exponent for calculating
#'     cladogenesis rate}
#'     \item{[2]: is x the exponent for calculating extinction rate}
#'   }
#' @param island_ontogeny In \code{\link{DAISIE_sim_time_dep}()},
#'   \code{\link{DAISIE_ML_CS}} and plotting a string describing the type of
#'   island ontogeny. Can be \code{"const"}, \code{"beta"} for a beta function
#'   describing area through time. String checked by
#'   \code{\link{is_island_ontogeny_input}()}. \cr In all other functions a
#'   numeric describing the type of island ontogeny. Can be \code{0} for
#'   constant, \code{1} for a beta function describing area through time. In ML
#'   functions \code{island_ontogeny = NA} assumes constant ontogeny.
#' @param sea_level In \code{\link{DAISIE_sim_time_dep}()} and plotting a
#'   string describing the type of sea level. Can be \code{"const"} or
#'   \code{"sine"} for a sine function describing area through time. String
#'   checked by \code{\link{is_sea_level_input}()}.
#'   \cr In all other functions a numeric describing the type of sea level. Can
#'   be \code{0} for constant, \code{1} for a sine function describing area
#'   through time.
#' @param extcutoff A numeric with the cutoff for the the maximum extinction
#'   rate preventing it from being too large and slowing down simulation.
#' @param shift_times a numeric vector specifying when the rate shifts occur
#'   before the present.
#' @param mainland_n A numeric stating the number of mainland species, that
#'   is the number of species that can potentially colonize the island.
#'   If using a clade-specific diversity dependence, this value is set to 1.
#'   If using an island-wide diversity dependence, this value is set to the
#'   number of mainland species.
#' @param island_replicates List output from
#'   \code{\link{DAISIE_sim_core_cr}()},
#'   \code{\link{DAISIE_sim_core_time_dep}()},
#'   \code{\link{DAISIE_sim_core_cr_shift}()} or
#'   \code{\link{DAISIE_sim_min_type2}()} functions. Minimally, this must be a
#'   list that has as many elements as replicates. Each element must be a list
#'   with the elements \code{island_age}, \code{not_present} and \code{stt_all}.
#'   \code{stt_all} must be a data frame with the column names \code{Time},
#'   \code{nI}, \code{nA}, \code{nC} and \code{present}.
#' @param island_spec Matrix with current state of simulation containing number
#'   of species.
#' @param stt_table Matrix with number of species at each time step.
#' @param rates named list of numeric rates as returned by
#'   \code{\link{update_rates}()}.
#' @param max_rates named list of numeric max rates as returned by
#'   \code{\link{update_max_rates}()}.
#' @param timeval Numeric defining current time of simulation.
#' @param total_time Numeric defining the length of the simulation in time
#'   units.
#' @param possible_event Numeric defining what event will happen.
#' @param maxspecID Current species IDs.
#' @param mainland_spec Number of mainland species.
#' @param max_area Numeric defining maximum area.
#' @param proportional_peak_t Numeric value from 0 to 1 indicating
#'   where in the island's history the peak area is achieved.
#' @param total_island_age Numeric defining total island age.
#' @param sea_level_amplitude Numeric defining amplitude of area fluctuation
#'   from sea level.
#' @param sea_level_frequency Numeric defining frequency of sine wave of
#'   area change from sea level.
#' @param island_gradient_angle Numeric defining the angle in degrees
#'   specifying the slope of the island.
#' @param d Numeric defining the scaling parameter for exponent for
#'   calculating cladogenesis rate.
#' @param x Numeric defining the exponent for calculating extinction rate.
#' @param simulation_outputs A list with matrices and vectors of simulation
#'   produced by DAISIE_sim functions.
#' @param plot_plus_one Boolean to indicate to plot all values plus one.
#'   Set to \code{TRUE} for default behavior. Set to \code{FALSE} to plot all
#'   values without adding one. Only works when there is one type of species.
#' @param type String to indicate if stt of all species or all possible stt
#'   should be plotted. Default is \code{"all_species"}, \code{"type1_species"}
#'   or \code{"type2_species"} should be plotted.
#' @param plot_lists List of lists containing average and quantile species
#'   through time.
#' @param ... Any arguments to pass on to plotting functions.
#' @param datalist Data object containing information on colonisation and
#'   branching times. This object can be generated using the DAISIE_dataprep
#'   function, which converts a user-specified data table into a data object,
#'   but the object can of course also be entered directly.
#'   It is an R list object with the following elements.\cr The first element
#'   of the list has two or three components: \cr \cr \code{$island_age} - the
#'   island age \cr Then, depending on whether a distinction between types is
#'   made, we have:\cr \code{$not_present} - the number of mainland lineages
#'   that are not present on the island \cr or:\cr \code{$not_present_type1} -
#'   the number of mainland lineages of type 1 that are not present on the
#'   island \cr \code{$not_present_type2} - the number of mainland lineages of
#'   type 2 that are not present on the island \cr \cr The remaining elements of
#'   the list each contains information on a single colonist lineage on the
#'   island and has 5 components:\cr \cr \code{$colonist_name} - the name of the
#'   species or clade that colonized the island \cr \code{$branching_times} -
#'   island age followed by stem age of the population/species in the case of
#'   Non-endemic, Non-endemic_MaxAge species and Endemic species with no close relatives
#'   on the island. For endemic clades with more than one species on the island
#'   (cladogenetic clades/ radiations) these should be island age followed by the
#'   branching times of the island clade
#'   including the stem age of the clade\cr \code{$stac} - the
#'   status of the colonist \cr \cr * Non_endemic_MaxAge: 1 \cr * Endemic: 2
#'   \cr * Endemic&Non_Endemic: 3 \cr * Non_Endemic: 4 \cr
#'   * Endemic_Singleton_MaxAge: 5 \cr * Endemic_Clade_MaxAge: 6
#'   \cr * Endemic&Non_Endemic_Clade_MaxAge: 7 \cr
#'   \cr \code{$missing_species} - number of island species that were not
#'   sampled for particular clade (only applicable for endemic clades) \cr
#'   \code{$type1or2} - whether the colonist belongs to type 1 or type 2 \cr
#' @param datatype Sets the type of data: 'single' for a single island or
#'   archipelago treated as one, and 'multiple' for multiple archipelagoes
#'   potentially sharing the same parameters.
#' @param initparsopt The initial values of the parameters that must be
#'   optimized, they are all positive.
#' @param idparsopt The ids of the parameters that must be optimized. The ids
#'   are defined as follows: \cr \cr id = 1 corresponds to lambda^c
#'   (cladogenesis rate) \cr id = 2 corresponds to mu (extinction rate) \cr
#'   id = 3 corresponds to K (clade-level carrying capacity) \cr id = 4
#'   corresponds to gamma (immigration rate) \cr id = 5 corresponds to lambda^a
#'   (anagenesis rate) \cr id = 6 corresponds to lambda^c (cladogenesis rate)
#'   for an optional subset of the species \cr id = 7 corresponds to mu
#'   (extinction rate) for an optional subset of the species\cr id = 8
#'   corresponds to K (clade-level carrying capacity) for an optional subset of
#'   the species\cr id = 9 corresponds to gamma (immigration rate) for an
#'   optional subset of the species\cr id = 10 corresponds to lambda^a
#'   (anagenesis rate) for an optional subset of the species\cr id = 11
#'   corresponds to p_f (fraction of mainland species that belongs to the second
#'   subset of species.
#' @param idparsfix The ids of the parameters that should not be optimized,
#'  e.g. c(1,3) if lambda^c and K should not be optimized.
#' @param parsfix The values of the parameters that should not be optimized.
#' @param idparsnoshift For datatype = 'single' only: The ids of the parameters
#'   that should not be different between two groups of species; This can only
#'   apply to ids 6:10, e.g. idparsnoshift = c(6,7) means that lambda^c and mu
#'   have the same values for both groups.
#' @param idparsmat For datatype = 'multiple' only: Matrix containing the ids
#'   of the parameters, linking them to initparsopt and parsfix. Per island
#'   system we use the following order: \cr \cr * lac = (initial) cladogenesis
#'   rate \cr * mu = extinction rate \cr * K = maximum number of species possible
#'   in the clade \cr * gam = (initial) immigration rate \cr * laa = (initial)
#'   anagenesis rate \cr Example:
#'   \code{idparsmat = rbind(c(1, 2, 3, 4, 5), c(1, 2, 3, 6, 7))} has different
#'   rates of immigration and anagenesis for the two islands.
#' @param res Sets the maximum number of species for which a probability must
#'   be computed, must be larger than the size of the largest clade.
#' @param ddmodel Sets the model of diversity-dependence: \cr \cr ddmodel = 0 :
#'   no diversity dependence \cr ddmodel = 1 : linear dependence in speciation
#'   rate \cr ddmodel = 11: linear dependence in speciation rate and in
#'   immigration rate \cr ddmodel = 2 : exponential dependence in speciation
#'   rate\cr ddmodel = 21: exponential dependence in speciation rate and in
#'   immigration rate\cr
#' @param cond cond = 0 : conditioning on island age \cr cond = 1 :
#'   conditioning on island age and non-extinction of the island biota \cr.
#'   cond > 1 : conditioning on island age and having at least cond colonizations
#'   on the island. This last option is not yet available for the IW model \cr
#' @param eqmodel Sets the equilibrium constraint that can be used during the
#'   likelihood optimization. Only available for datatype = 'single'.\cr\cr
#'   eqmodel = 0 : no equilibrium is assumed \cr eqmodel = 13 : near-equilibrium
#'   is assumed on endemics using deterministic equation for endemics and
#'   immigrants. Endemics must be within x_E of the equilibrium value\cr eqmodel
#'   = 15 : near-equilibrium is assumed on endemics and immigrants using
#'   deterministic equation for endemics and immigrants. Endemics must be within
#'   x_E of the equilibrium value, while non-endemics must be within x_I of the
#'   equilibrium value.
#' @param x_E Sets the fraction of the equlibrium endemic diversity above which
#'   the endemics are assumed to be in equilibrium; only active for eqmodel = 13
#'   or 15.
#' @param x_I Sets the fraction of the equlibrium non-endemic diversity above
#'   which the system is assumed to be in equilibrium; only active for eqmodel =
#'   15.
#' @param tol Sets the tolerances in the optimization. Consists of: \cr reltolx
#'   = relative tolerance of parameter values in optimization \cr reltolf =
#'   relative tolerance of function value in optimization \cr abstolx = absolute
#'   tolerance of parameter values in optimization.
#' @param maxiter Sets the maximum number of iterations in the optimization.
#' @param methode Method of the ODE-solver. Supported Boost \code{ODEINT}
#'   solvers (steppers) are:
#'   \code{"odeint::runge_kutta_cash_karp54"}
#'   \code{"odeint::runge_kutta_fehlberg78"}
#'   \code{"odeint::runge_kutta_dopri5"}
#'   \code{"odeint::bulirsch_stoer"}
#'   without \code{odeint::}-prefix, \code{\link[deSolve]{ode}} method is
#'   assumed. The default method overall is
#'   \code{"lsodes"} for \code{\link{DAISIE_ML_CS}()}
#'   and \code{"ode45"} from \code{\link[deSolve]{ode}()} for
#'   \code{\link{DAISIE_ML_IW}()}.
#' @param optimmethod Method used in likelihood optimization. Default is
#'   `subplex` (see `\link[subplex]{subplex}()` for full details).
#'   Alternative is \code{"simplex"} which was the method in previous versions.
#' @param tolint Vector of two elements containing the absolute and relative
#'   tolerance of the integration.
#' @param datatable Data frame (table) with user-specified data. See file
#'   \code{Galapagos_datatable} for a template of an input table. Each row on the
#'   table represents and independent colonisation event. Table has the
#'   following four columns. \cr \cr \code{$Clade_name} - name of independent
#'   colonization event \cr \code{$Status} - One of the following categories:
#'   \cr * "Non_endemic": applies to non-endemic species for cases where both
#'   island and non-island populations of the species have been sampled) \cr *
#'   "Non_endemic_MaxAge": applies to non-endemic species for cases where island
#'   individuals of the species have not been sampled and only the age of the
#'   species is available) \cr * "Endemic": applies to endemic species and is
#'   applicable for both cladogenetic and anagenetic species \cr *
#'   "Endemic_MaxAge": applies to endemic species for cases where island
#'   individuals of the species have not been sampled and only the age of the
#'   species is available. This could apply to endemic species that have
#'   recently gone extinct because of antropogenic causes that are (evidently)
#'   not modelled, and for which no DNA data is available.\cr *
#'   "Endemic&Non_Endemic": when endemic clade is present and its mainland
#'   ancestor has re-colonized \cr \code{$Missing_species} - Number of island
#'   species that were not sampled for particular clade (only applicable for
#'   "Endemic" clades) \cr \code{$Branching_times} - Stem age of the
#'   population/species in the case of "Non-endemic", "Non-endemic_MaxAge" and
#'   "Endemic" anagenetic species. For "Endemic" cladogenetic species these
#'   should be branching times of the radiation including the stem age of the
#'   radiation.
#' @param island_age Age of island in appropriate units. In
#'   \code{\link{DAISIE_plot_age_diversity}()} and
#'   \code{\link{DAISIE_plot_island}()} if island input is in table format,
#'   the age of the island must be specified. If island input is in DAISIE list
#'   format, this option will override the island age specified in the island
#'   list.
#' @param number_clade_types Number of clade types. Default: number_clade_types
#'   = 1 all species are considered to belong to same macroevolutionary process.
#'   If number_clade_types = 2, there are two types of clades with distinct
#'   macroevolutionary processes.
#' @param list_type2_clades If \code{number_clade_types = 2}, list_type2_clades
#'   specifies the names of the clades that have a distinct macroevolutionary
#'   process. The names must match those in the $Clade_name column of the source
#'   data table (e.g. \code{list_type2_clades = "Finches"}).  If
#'   \code{number_clade_types = 1}, then list_type2_clades = NA should be
#'   specified (default).
#' @param epss Default= 1E-5 should be appropriate in most cases. This value
#'   is used to set the maximum age of colonisation of "Non_endemic_MaxAge" and
#'   "Endemic_MaxAge" species to an age that is slightly younger than the island
#'   for cases when the age provided for that species is older than the island.
#'   The new maximum age is then used as an upper bound to integrate over all.
#' @param t The time at which the expectations need to be computed.
#' @param initEI The initial values for the number of endemics and
#'   non-endemics. In \code{\link{DAISIE_probdist}()} or
#' \code{\link{DAISIE_margprobdist}()} either this or initprobs must be NULL. In
#' \code{\link{DAISIE_numcol}()} when it is NULL, it is assumed that the island
#'   is empty.
#' @param data_table data table
#' @param endmc Numeric for how many simulations should run.
#' @param archipelago something
#' @param phylo_data  something
#' @param archipelago_data  something
#' @param gam A numeric with the per capita immigration rate.
#' @param laa A numeric with the per capita anagenesis rate.
#' @param lac A numeric with the per capita cladogenesis rate.
#' @param mu A numeric with the per capita extinction rate.
#' @param K A numeric with carrying capacity.
#' @param num_spec A numeric with the current number of species.
#' @param num_immigrants A numeric with the current number of non-endemic
#' species (a.k.a non-endemic species).
#' @param global_min_area_time stub
#' @param global_max_area_time stub
#' @param distance_type Use 'continent' if the distance to the continent should
#'   be used, use 'nearest_big' if the distance to the nearest big landmass
#'   should be used, and use 'biologically_realistic' if the distance should
#'   take into account some biologically realism, e.g. an average of the
#'   previous two if both are thought to contribute.
#' @param distance_dep Sets what type of distance dependence should be used.
#'   Default is a power law, denoted as 'power'. Alternatives are an exponantial
#'   relationship denoted by 'exp' or sigmoids, either 'sigmoidal_col' for a
#'   sigmoid in the colonization, 'sigmoidal_ana' for sigmoidal anagenesis,
#'   'sigmoidal_clado' for sigmoidal cladogenesis, and 'sigmoidal_col_ana' for
#'   signoids in both colonization and anagenesis.
#' @param parallel Sets whether parallel computation should be used. Use 'no'
#'   if no parallel computing should be used, 'cluster' for parallel computing
#'   on a unix/linux cluster, and 'local' for parallel computation on a local
#'   machine.
#' @param cpus Number of cpus used in parallel computing. Default is 3. Will
#'   not have an effect if parallel = 'no'.
#' @param pars1 Vector of model parameters: \cr \cr \code{pars1[1]} corresponds
#'   to lambda^c (cladogenesis rate) \cr \code{pars1[2]} corresponds to mu
#'   (extinction rate) \cr \code{pars1[3]} corresponds to K (clade-level
#'   carrying capacity) \cr \code{pars1[4]} corresponds to gamma
#'   (immigration rate) \cr \code{pars1[5]} corresponds to lambda^a
#'   (anagenesis rate).
#' @param pars2 Vector of settings: \cr \cr \code{pars2[1]} corresponds to res,
#'   the maximum number of endemics or non-endemics for which the ODE system is
#'   solved; this must be much larger than the actual number for which the
#'   probability needs to be calculated.) \cr \code{pars2[2]} corresponds to M,
#'   size of the mainland pool, i.e the number of species that can potentially
#'   colonize the island.
#' @param tvec The times at which the probabilities need to be computed.
#' @param initprobs The initial probability distribution for the number of
#'   endemics and non-endemics; either this or initEI must be NULL.
#' @param pb Probability distribution in matrix format as output by
#'   \code{\link{DAISIE_probdist}()}.
#' @param island Island data object. Can be in DAISIE list format (see
#'   Galapagos_datalist and DAISIE_data_prep for examples) or in table format
#'   (see Galapagos_datatable for an example).
#' @param title Title of the plot
#' @param plot_lists_simulations List with simulation output after parsing by
#' \code{DAISIE_prepare_data_plotting}.
#' @param plot_lists_simulations_MLE List with simulation output after parsing
#'   by \code{DAISIE_prepare_data_plotting}, but obtained by simulating MLE
#'   output.
#' @param kind_of_plot Character vector stating how STT plot resulting from MLE
#'   based simulations should be plotted. Default is \code{"line"} for multiple
#'   individual lines. Can also be \code{"shade"} for the 5\% quantile.
#' @param resolution numeric indicating resolution of plot. Should be < 0.
#' @param resol numeric for resolution of summary stats calculation. Should be
#'   > 1.
#' @param removed_timepoints Positive integer with number of first datapoints
#'   to be removed from rate plots (to prevent Inf)
#' @param A A numeric value for island area at a given point in time.
#' @param Amin A numeric value for minimum island area during the simulation.
#' @param Amax A numeric value for maximum island area during the simulation.
#' @param peak A numeric value specifying the peakiness (or shaprness) of the
#'   ontogeny curve. Higher values imply peakier ontogeny. This value is
#'   internally calculated by \code{\link{calc_peak}()} given the area at the
#'   present and the \code{area_pars}.
#' @param proptime A numeric from 0 to 1. The proportion of time that has
#'   elapsed in the simulation, in relation to the total island age (NB: not
#'   the simulation time, but island age).
#' @param proptime_max A numeric from 0 to 1. The same as
#'   \code{proportional_peak_t}. Indicates, in proportion to the total island
#'   age when the ontogeny peak should occur (i.e. 0.5 means a peak halfway in
#'   time).
#' @param current_area A numeric with the current island area at present (i.e.,
#'   at the end of the simulation).
#' @param jitter Numeric for \code{\link[DDD]{optimizer}()}. Jitters the
#'   parameters being optimized by the specified amount which should be very
#'   small, e.g. 1e-5. Jitter when \code{link[subplex]{subplex}()} produces
#'   incorrect output due to parameter transformation.
#' @param num_cycles The number of cycles the optimizer will go through.
#'   Default is 1.
#' @param trait_pars A named list containing diversification rates considering
#'   two trait states created by \code{\link{create_trait_pars}}:
#' \itemize{
#'   \item{[1]:A numeric with the per capita transition rate with state1}
#'   \item{[2]:A numeric with the per capita immigration rate with state2}
#'   \item{[3]:A numeric with the per capita extinction rate with state2}
#'   \item{[4]:A numeric with the per capita anagenesis rate with state2}
#'   \item{[5]:A numeric with the per capita cladogenesis rate with state2}
#'   \item{[6]:A numeric with the per capita transition rate with state2}
#'   \item{[7]:A numeric with the number of species with trait state 2 on
#'    mainland}
#' }
#' @param relaxed_par A string determining which parameter is relaxed in a
#' relaxed rate model.
#' @param relaxed_rate_pars A two element list the first element is the mean
#' of the gamma distribution, the second element is the standard deviation of
#' the gamma distribution. List can be created with
#' \code{create_relaxed_rate_pars()}
#' @param brts Numeric vector of branching times
#' @param stac Numeric of Endemicity status
#' @param missnumspec Numeric of missing species
#' @param CS_version a numeric or list. Default is 1 for the standard DAISIE
#' model, for a relaxed-rate model a list with the following elements:
#' \itemize{
#'   \item{model: the CS model to run, options are \code{1} for single rate
#'   DAISIE model, \code{2} for multi-rate DAISIE, or \code{0} for IW test
#'   model.}
#'   \item{relaxed_par: the parameter to relax (integrate over). Options are
#' \code{"cladogenesis"}, \code{"extinction"}, \code{"carrying_capacity"},
#' \code{"immigration"}, or \code{"anagenesis"}.}
#'   }
#' @param DAISIE_par A numeric parameter to evaluate the integral of the
#' function.
#' @param DAISIE_dist_pars A numeric vector of two elements, first is the mean
#' and second the standard deviation of the distribution.
#' @param abstolint Numeric absolute tolerance of the integration
#' @param reltolint Numeric relative tolerance of the integration
#' @param pick Numeric determining which parameter is selected for the
#' relaxed-rate model
#' @param mean Numeric mean of the distribution
#' @param sd Numeric standard deviation of the distribution
#' @keywords internal
#' @note This is an internal function, so it should be marked with
#'   \code{@noRd}. This is not done, as this will disallow all
#'   functions to find the documentation parameters
#' @param clado_rate Numeric rate of cladogenesis
#' @param ext_rate Numeric rate of extinction
#' @param carr_cap Numeric carrying capacity
#' @param immig_rate Numeric rate of immigration
#' @param ana_rate Numeric rate of anagenesis
#'
#'
#' @return Nothing
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
  island_ontogeny,
  sea_level,
  extcutoff,
  shift_times,
  mainland_n,
  island_replicates,
  island_spec,
  stt_table,
  rates,
  max_rates,
  timeval,
  total_time,
  possible_event,
  maxspecID,
  mainland_spec,
  max_area,
  proportional_peak_t,
  total_island_age,
  sea_level_amplitude,
  sea_level_frequency,
  island_gradient_angle,
  d,
  x,
  simulation_outputs,
  plot_plus_one,
  type,
  plot_lists,
  ...,
  datalist,
  datatype,
  initparsopt,
  idparsopt,
  idparsfix,
  parsfix,
  idparsnoshift,
  idparsmat,
  res,
  ddmodel,
  cond,
  eqmodel,
  x_E,
  x_I,
  tol,
  maxiter,
  methode,
  optimmethod,
  CS_version,
  tolint,
  datatable,
  island_age,
  number_clade_types,
  list_type2_clades,
  epss,
  t,
  initEI,
  data_table,
  endmc,
  archipelago,
  phylo_data,
  archipelago_data,
  gam,
  laa,
  lac,
  mu,
  K,
  num_spec,
  num_immigrants,
  global_min_area_time,
  global_max_area_time,
  distance_type,
  distance_dep,
  parallel,
  cpus,
  pars1,
  pars2,
  tvec,
  initprobs,
  pb,
  island,
  title,
  plot_lists_simulations,
  plot_lists_simulations_MLE,
  kind_of_plot,
  resolution,
  resol,
  removed_timepoints,
  A,
  Amin,
  Amax,
  peak,
  proptime,
  proptime_max,
  current_area,
  jitter,
  num_cycles,
  trait_pars,
  relaxed_par,
  relaxed_rate_pars,
  brts,
  stac,
  missnumspec,
  DAISIE_par,
  DAISIE_dist_pars,
  abstolint,
  reltolint,
  pick,
  mean,
  sd,
  clado_rate,
  ext_rate,
  carr_cap,
  immig_rate,
  ana_rate
) {
  # Nothing
}
