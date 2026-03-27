# Default parameter documentation

Default parameter documentation

## Usage

``` r
default_params_doc(
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
  initEI_mat,
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
  function_to_optimize,
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
  ana_rate,
  islands,
  sort_clade_sizes,
  equal_extinction,
  files_to_write,
  use_rcpp
)
```

## Arguments

- time:

  Numeric defining the length of the simulation in time units. For
  example, if an island is known to be 4 million years old, setting time
  = 4 will simulate the entire life span of the island; setting time = 2
  will stop the simulation at the mid-life of the island.

- M:

  Numeric defining the size of mainland pool, i.e. the number of species
  that can potentially colonize the island.

- pars:

  A numeric vector containing the model parameters:

  - `pars[1]`: lambda^c (cladogenesis rate)

  - `pars[2]`: mu (extinction rate)

  - `pars[3]`: K (carrying capacity), set K=Inf for diversity
    independence.

  - `pars[4]`: gamma (immigration rate)

  - `pars[5]`: lambda^a (anagenesis rate)

  - `pars[6]`: lambda^c (cladogenesis rate) for either type 2 species or
    rate set 2 in rate shift model

  - `pars[7]`: mu (extinction rate) for either type 2 species or rate
    set 2 in rate shift model

  - `pars[8]`: K (carrying capacity) for either type 2 species or rate
    set 2 in rate shift model, set K=Inf for diversity independence.

  - `pars[9]`: gamma (immigration rate) for either type 2 species or
    rate set 2 in rate shift model

  - `pars[10]`: lambda^a (anagenesis rate) for either type 2 species or
    rate set 2 in rate shift model

  Elements 6:10 are required only when type 2 species are included or in
  the rate shift model. For
  [`DAISIE_sim_relaxed_rate()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_relaxed_rate.md)
  `pars[6]` is the standard deviation of the gamma distribution for the
  relaxed parameter and the parameter chosen by the `relaxed_par`
  argument is the mean of the gamma distribution for the relaxed
  parameter.

- replicates:

  Integer specifying number of island replicates to be simulated.

- divdepmodel:

  Option divdepmodel = 'CS' runs a model with clade-specific carrying
  capacity, where diversity-dependence operates only within single
  clades, i.e. only among species originating from the same mainland
  colonist. Option divdepmodel = 'IW' runs a model with island-wide
  carrying capacity, where diversity-dependence operates within and
  among clades. Option divdepmodel = 'GW' runs a model with
  diversity-dependence operates within a guild.

- nonoceanic_pars:

  A vector of length two with:

  - \[1\]: the probability of sampling a species from the mainland

  - \[2\]: the probability of the species sampled from the mainland
    being nonendemic

- num_guilds:

  The number of guilds on the mainland. The number of mainland species
  is divided by the number of guilds when `divdepmodel = "GW"`

- prop_type2_pool:

  Fraction of mainland species that belongs to the second subset of
  species (type 2). Applies only when two types of species are simulated
  (length(pars) = 10). For
  [`DAISIE_dataprep()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_dataprep.md)
  applies only if number_clade_types = 2. In
  [`DAISIE_dataprep()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_dataprep.md)
  the default `"proportional"` sets the fraction to be proportional to
  the number of clades of distinct macroevolutionary process that have
  colonised the island.

- replicates_apply_type2:

  Applies only when two types of species are being simulated. Default
  replicates_apply_type2 = TRUE runs simulations until the number of
  islands where a type 2 species has colonised is equal to the specified
  number of replicates. This is recommended if `prop_type2_pool` is
  small or if the rate of immigration of type two species (`pars[9]`) is
  low, meaning that more replicates are needed to achieved an adequate
  sample size of islands with type 2 species. Setting
  `replicates_apply_type2 = FALSE` simulates islands up to the specified
  number of replicates regardless of whether type 2 species have
  colonised or not.

- sample_freq:

  Numeric specifing the number of units times should be divided by for
  plotting purposes. Larger values will lead to plots with higher
  resolution, but will also run slower.

- plot_sims:

  `Default = TRUE` plots species-through-time (STT) plots. It detects
  how many types of species are present. If only one type of species is
  present, STT is plotted for all species. If two types are present,
  three plots are produced: STT for all, STT for type 1 and STT for type
  2.

- verbose:

  A numeric vector of length 1, which in simulations and
  \`DAISIEdataprep()\` can be \`1\` or \`0\`, where \`1\` gives
  intermediate output should be printed. For ML functions a numeric
  determining if intermediate output should be printed. The default:
  \`0\` does not print, \`1\` prints the initial likelihood and the
  settings that were selected (which parameters are to be optimised,
  fixed or shifted), \`2\` prints the same as \`1 and also the
  intermediate output of the parameters and loglikelihood, while \`3\`
  the same as \`2\` and prints intermediate progress during likelihood
  computation.

- area_pars:

  A named list containing area and sea level parameters as created by
  [`create_area_pars()`](https://rsetienne.github.io/DAISIE/reference/create_area_pars.md):

  - \[1\]: maximum area

  - \[2\]: current area

  - \[3\]: value from 0 to 1 indicating where in the island's history
    the peak area is achieved

  - \[4\]: total island age

  - \[5\]: amplitude of area fluctuation from sea level

  - \[6\]: frequency of sine wave of area change from sea level

  - \[7\]: angle of the slope of the island

- hyper_pars:

  A named list of numeric hyperparameters for the rate calculations as
  returned by
  [`create_hyper_pars()`](https://rsetienne.github.io/DAISIE/reference/create_hyper_pars.md):

  - \[1\]: is d the scaling parameter for exponent for calculating
    cladogenesis rate

  - \[2\]: is x the exponent for calculating extinction rate

- island_ontogeny:

  In
  [`DAISIE_sim_time_dep()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_time_dep.md),
  [`DAISIE_ML_CS`](https://rsetienne.github.io/DAISIE/reference/DAISIE_ML.md)
  and plotting a string describing the type of island ontogeny. Can be
  `"const"`, `"beta"` for a beta function describing area through
  time.  
  In all other functions a numeric describing the type of island
  ontogeny. Can be `0` for constant, `1` for a beta function describing
  area through time. In ML functions `island_ontogeny = NA` assumes
  constant ontogeny. Time dependent estimation is not yet available as
  development is still ongoing. Will return an error if called in that
  case.

- sea_level:

  In
  [`DAISIE_sim_time_dep()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_time_dep.md)
  and plotting a string describing the type of sea level. Can be
  `"const"` or `"sine"` for a sine function describing area through
  time.  
  In all other functions a numeric describing the type of sea level. Can
  be `0` for constant, `1` for a sine function describing area through
  time.

- extcutoff:

  A numeric with the cutoff for the the maximum extinction rate
  preventing it from being too large and slowing down simulation.

- shift_times:

  a numeric vector specifying when the rate shifts occur before the
  present.

- mainland_n:

  A numeric stating the number of mainland species, that is the number
  of species that can potentially colonize the island. If using a
  clade-specific diversity dependence, this value is set to 1. If using
  an island-wide diversity dependence, this value is set to the number
  of mainland species.

- island_replicates:

  List output from
  [`DAISIE_sim_core_cr()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_core_cr.md),
  [`DAISIE_sim_core_time_dep()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_core_time_dep.md),
  [`DAISIE_sim_core_cr_shift()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_core_cr_shift.md)
  or
  [`DAISIE_sim_min_type2()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_min_type2.md)
  functions. Minimally, this must be a list that has as many elements as
  replicates. Each element must be a list with the elements
  `island_age`, `not_present` and `stt_all`. `stt_all` must be a data
  frame with the column names `Time`, `nI`, `nA`, `nC` and `present`.

- island_spec:

  Matrix with current state of simulation containing number of species.

- stt_table:

  Matrix with number of species at each time step.

- rates:

  named list of numeric rates as returned by
  [`update_rates()`](https://rsetienne.github.io/DAISIE/reference/update_rates.md).

- max_rates:

  named list of numeric max rates as returned by
  [`update_max_rates()`](https://rsetienne.github.io/DAISIE/reference/update_max_rates.md).

- timeval:

  Numeric defining current time of simulation.

- total_time:

  Numeric defining the length of the simulation in time units.

- possible_event:

  Numeric defining what event will happen.

- maxspecID:

  Current species IDs.

- mainland_spec:

  Number of mainland species.

- max_area:

  Numeric defining maximum area.

- proportional_peak_t:

  Numeric value from 0 to 1 indicating where in the island's history the
  peak area is achieved.

- total_island_age:

  Numeric defining total island age.

- sea_level_amplitude:

  Numeric defining amplitude of area fluctuation from sea level.

- sea_level_frequency:

  Numeric defining frequency of sine wave of area change from sea level.

- island_gradient_angle:

  Numeric defining the angle in degrees specifying the slope of the
  island.

- d:

  Numeric defining the scaling parameter for exponent for calculating
  cladogenesis rate.

- x:

  Numeric defining the exponent for calculating extinction rate.

- simulation_outputs:

  A list with matrices and vectors of simulation produced by DAISIE_sim
  functions.

- plot_plus_one:

  Boolean to indicate to plot all values plus one. Set to `TRUE` for
  default behavior. Set to `FALSE` to plot all values without adding
  one. Only works when there is one type of species.

- type:

  String to indicate if stt of all species or all possible stt should be
  plotted. Default is `"all_species"`, `"type1_species"` or
  `"type2_species"` should be plotted.

- plot_lists:

  List of lists containing average and quantile species through time.

- ...:

  Any arguments to pass on to plotting functions.

- datalist:

  Data object containing information on colonisation and branching
  times. This object can be generated using the DAISIE_dataprep
  function, which converts a user-specified data table into a data
  object, but the object can of course also be entered directly. It is
  an R list object with the following elements.  
  The first element of the list has two or three components:  
    
  `$island_age` - the island age  
  Then, depending on whether a distinction between types is made, we
  have:  
  `$not_present` - the number of mainland lineages that are not present
  on the island  
  or:  
  `$not_present_type1` - the number of mainland lineages of type 1 that
  are not present on the island  
  `$not_present_type2` - the number of mainland lineages of type 2 that
  are not present on the island  
    
  The remaining elements of the list each contains information on a
  single colonist lineage on the island and has 5 components:  
    
  `$colonist_name` - the name of the species or clade that colonized the
  island  
  `$branching_times` - island age followed by stem age of the
  population/species in the case of Non-endemic, Non-endemic_MaxAge
  species and Endemic species with no close relatives on the island. For
  endemic clades with more than one species on the island (cladogenetic
  clades/ radiations) these should be island age followed by the
  branching times of the island clade including the stem age of the
  clade  
  `$stac` - the status of the colonist  
    
  - Non_endemic_MaxAge: 1  
  - Endemic: 2  
  - Endemic&Non_Endemic: 3  
  - Non_Endemic: 4  
  - Endemic_Singleton_MaxAge: 5  
  - Endemic_Clade_MaxAge: 6  
  - Endemic&Non_Endemic_Clade_MaxAge: 7  
  - Non_endemic_MaxAge_MinAge: 8  
  - Endemic_Singleton_MaxAge_MinAge: 9  
    
  `$missing_species` - number of island species that were not sampled
  for particular clade (only applicable for endemic clades)  
  `$type1or2` - whether the colonist belongs to type 1 or type 2  

- datatype:

  Sets the type of data: 'single' for a single island or archipelago
  treated as one, and 'multiple' for multiple archipelagoes potentially
  sharing the same parameters.

- initparsopt:

  The initial values of the parameters that must be optimized, they are
  all positive.

- idparsopt:

  The ids of the parameters that must be optimized. The ids are defined
  as follows:  
    
  id = 1 corresponds to lambda^c (cladogenesis rate)  
  id = 2 corresponds to mu (extinction rate)  
  id = 3 corresponds to K (clade-level carrying capacity)  
  id = 4 corresponds to gamma (immigration rate)  
  id = 5 corresponds to lambda^a (anagenesis rate)  
  id = 6 corresponds to lambda^c (cladogenesis rate) for an optional
  subset of the species  
  id = 7 corresponds to mu (extinction rate) for an optional subset of
  the species  
  id = 8 corresponds to K (clade-level carrying capacity) for an
  optional subset of the species  
  id = 9 corresponds to gamma (immigration rate) for an optional subset
  of the species  
  id = 10 corresponds to lambda^a (anagenesis rate) for an optional
  subset of the species  
  id = 11 corresponds to p_f (fraction of mainland species that belongs
  to the second subset of species.

- idparsfix:

  The ids of the parameters that should not be optimized, e.g. c(1,3) if
  lambda^c and K should not be optimized.

- parsfix:

  The values of the parameters that should not be optimized.

- idparsnoshift:

  For datatype = 'single' only: The ids of the parameters that should
  not be different between two groups of species; This can only apply to
  ids 6:10, e.g. idparsnoshift = c(6,7) means that lambda^c and mu have
  the same values for both groups.

- idparsmat:

  For datatype = 'multiple' only: Matrix containing the ids of the
  parameters, linking them to initparsopt and parsfix. Per island system
  we use the following order:  
    
  \* lac = (initial) cladogenesis rate  
  \* mu = extinction rate  
  \* K = maximum number of species possible in the clade  
  \* gam = (initial) immigration rate  
  \* laa = (initial) anagenesis rate  
  Example: `idparsmat = rbind(c(1, 2, 3, 4, 5), c(1, 2, 3, 6, 7))` has
  different rates of immigration and anagenesis for the two islands.

- res:

  Sets the maximum number of species for which a probability must be
  computed, must be larger than the size of the largest clade.

- ddmodel:

  Sets the model of diversity-dependence:  
    

  - ddmodel = 0 : no diversity dependence

  - ddmodel = 1 : linear dependence in speciation rate

  - ddmodel = 11: linear dependence in speciation rate and in
    immigration rate

  - ddmodel = 2 : exponential dependence in speciation rate

  - ddmodel = 21: exponential dependence in speciation rate and in
    immigration rate

- cond:

  cond = 0 : conditioning on island age  
  cond = 1 : conditioning on island age and non-extinction of the island
  biota  
  . cond \> 1 : conditioning on island age and having at least cond
  colonizations on the island. This last option is not yet available for
  the IW model  

- eqmodel:

  Sets the equilibrium constraint that can be used during the likelihood
  optimization. Only available for datatype = 'single'.  
    
  eqmodel = 0 : no equilibrium is assumed  
  eqmodel = 13 : near-equilibrium is assumed on endemics using
  deterministic equation for endemics and immigrants. Endemics must be
  within x_E of the equilibrium value  
  eqmodel = 15 : near-equilibrium is assumed on endemics and immigrants
  using deterministic equation for endemics and immigrants. Endemics
  must be within x_E of the equilibrium value, while non-endemics must
  be within x_I of the equilibrium value.

- x_E:

  Sets the fraction of the equlibrium endemic diversity above which the
  endemics are assumed to be in equilibrium; only active for eqmodel =
  13 or 15.

- x_I:

  Sets the fraction of the equlibrium non-endemic diversity above which
  the system is assumed to be in equilibrium; only active for eqmodel =
  15.

- tol:

  Sets the tolerances in the optimization. Consists of:  
  reltolx = relative tolerance of parameter values in optimization  
  reltolf = relative tolerance of function value in optimization  
  abstolx = absolute tolerance of parameter values in optimization.

- maxiter:

  Sets the maximum number of iterations in the optimization.

- methode:

  Method of the ODE-solver. Supported Boost `ODEINT` solvers (steppers)
  are: `"odeint::runge_kutta_cash_karp54"`
  `"odeint::runge_kutta_fehlberg78"` `"odeint::runge_kutta_dopri5"`
  `"odeint::bulirsch_stoer"` without `odeint::`-prefix,
  [`ode`](https://rdrr.io/pkg/deSolve/man/ode.html) method is assumed.
  The default method overall is `"lsodes"` for
  [`DAISIE_ML_CS()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_ML.md)
  and `"ode45"` from [`ode()`](https://rdrr.io/pkg/deSolve/man/ode.html)
  for
  [`DAISIE_ML_IW()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_ML_IW.md).

- optimmethod:

  Method used in likelihood optimization. Default is \`simplex\` in the
  standard Clade Specific scenario. Alternative is \`subplex\` (see
  \`[subplex](https://rdrr.io/pkg/subplex/man/subplex.html)()\` for full
  details) which was the default method in previous versions. In the
  Island Wide, two type scenarios, and split rate scenarios the default
  remains \`subplex\`.

- CS_version:

  a numeric or list. Default is CS_version = list(model = 1,
  function_to_optimize = 'DAISIE'), but for a relaxed-rate model the
  list can contain more elements:

  - model: the CS model to run, options are `1` for single rate DAISIE
    model, `2` for multi-rate DAISIE, or `0` for IW test model

  - function_to_optimize: the DAISIE loglikelihood function that will be
    optimized. Options are: `"DAISIE"`, default, the full DAISIE
    loglikelihood `"DAISIE_approx"`, an approximate loglikelihood
    `"DAISIE_DE"`, an exact loglikelkhood for K = Inf based on the D-E
    approach

  - integration_method: the method used to do integraion in the relaxed
    rate model. Options are: `'standard'` the default numerical
    integration `'MC'` Monte Carlo integration `'stratified'` using
    quantiles of the gamma distribution

  - relaxed_par: the parameter to relax (integrate over) in the relaxed
    rate model. Options are `"cladogenesis"`, `"extinction"`,
    `"carrying_capacity"`, `"immigration"`, or `"anagenesis"`

  - par_sd: standard deviation of the parameter to relax

  - par_upper_bound upper bound of the parameter to relax

  - seed: seed of the random number generator in case of 'MC'

  - sample_size: size of sample in case of 'MC' or 'stratified'

  - parallel: use parallel computing or not in case of 'MC' or
    'stratified'

  - n_cores: number of cores to use when run in parallel

- tolint:

  Vector of two elements containing the absolute and relative tolerance
  of the integration.

- datatable:

  Data frame (table) with user-specified data. See file
  `Galapagos_datatable` for a template of an input table. Each row on
  the table represents and independent colonisation event. Table has the
  following four columns.  
    
  `$Clade_name` - name of independent colonization event  
  `$Status` - One of the following categories:  
  \* "Non_endemic": applies to non-endemic species for cases where both
  island and non-island populations of the species have been sampled)  
  \* "Non_endemic_MaxAge": applies to non-endemic species for cases
  where island individuals of the species have not been sampled and only
  the age of the species is available)  
  \* "Endemic": applies to endemic species and is applicable for both
  cladogenetic and anagenetic species  
  \* "Endemic_MaxAge": applies to endemic species for cases where island
  individuals of the species have not been sampled and only the age of
  the species is available. This could apply to endemic species that
  have recently gone extinct because of antropogenic causes that are
  (evidently) not modelled, and for which no DNA data is available.  
  \* "Endemic&Non_Endemic": when endemic clade is present and its
  mainland ancestor has re-colonized  
  `$Missing_species` - Number of island species that were not sampled
  for particular clade (only applicable for "Endemic" clades)  
  `$Branching_times` - Stem age of the population/species in the case of
  "Non-endemic", "Non-endemic_MaxAge" and "Endemic" anagenetic species.
  For "Endemic" cladogenetic species these should be branching times of
  the radiation including the stem age of the radiation.

- island_age:

  Age of island in appropriate units. In
  [`DAISIE_plot_age_diversity()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_plot_age_diversity.md)
  and
  [`DAISIE_plot_island()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_plot_island.md)
  if island input is in table format, the age of the island must be
  specified. If island input is in DAISIE list format, this option will
  override the island age specified in the island list.

- number_clade_types:

  Number of clade types. Default: number_clade_types = 1 all species are
  considered to belong to same macroevolutionary process. If
  number_clade_types = 2, there are two types of clades with distinct
  macroevolutionary processes.

- list_type2_clades:

  If `number_clade_types = 2`, list_type2_clades specifies the names of
  the clades that have a distinct macroevolutionary process. The names
  must match those in the \$Clade_name column of the source data table
  (e.g. `list_type2_clades = "Finches"`). If `number_clade_types = 1`,
  then list_type2_clades = NA should be specified (default).

- epss:

  Default= 1E-5 should be appropriate in most cases. This value is used
  to set the maximum age of colonisation of "Non_endemic_MaxAge" and
  "Endemic_MaxAge" species to an age that is slightly younger than the
  island for cases when the age provided for that species is older than
  the island. The new maximum age is then used as an upper bound to
  integrate over all.

- initEI:

  The initial values for the number of endemics and non-endemics. In
  [`DAISIE_probdist()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_probdist.md)
  or
  [`DAISIE_margprobdist()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_margprobdist.md)
  either this or initprobs must be NULL. In
  [`DAISIE_numcol()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_numcol.md)
  when it is NULL, it is assumed that the island is empty.

- initEI_mat:

  matrix where each row represents the initial number of endemic and
  non-endemic species per colonizing lineage.

- data_table:

  data table

- endmc:

  Numeric for how many simulations should run.

- archipelago:

  something

- phylo_data:

  something

- archipelago_data:

  something

- gam:

  A numeric with the per capita immigration rate.

- laa:

  A numeric with the per capita anagenesis rate.

- lac:

  A numeric with the per capita cladogenesis rate.

- mu:

  A numeric with the per capita extinction rate.

- K:

  A numeric with carrying capacity.

- num_spec:

  A numeric with the current number of species.

- num_immigrants:

  A numeric with the current number of non-endemic species (a.k.a
  non-endemic species).

- global_min_area_time:

  stub

- global_max_area_time:

  stub

- distance_type:

  Use 'continent' if the distance to the continent should be used, use
  'nearest_big' if the distance to the nearest big landmass should be
  used, and use 'biologically_realistic' if the distance should take
  into account some biologically realism, e.g. an average of the
  previous two if both are thought to contribute.

- distance_dep:

  Sets what type of distance dependence should be used. Default is a
  power law, denoted as 'power'. Alternatives are an exponantial
  relationship denoted by 'exp' or sigmoids, either 'sigmoidal_col' for
  a sigmoid in the colonization, 'sigmoidal_ana' for sigmoidal
  anagenesis, 'sigmoidal_clado' for sigmoidal cladogenesis, and
  'sigmoidal_col_ana' for signoids in both colonization and anagenesis.

- parallel:

  Sets whether parallel computation should be used. Use 'no' if no
  parallel computing should be used, 'cluster' for parallel computing on
  a unix/linux cluster, and 'local' for parallel computation on a local
  machine.

- cpus:

  Number of cpus used in parallel computing. Default is 3. Will not have
  an effect if parallel = 'no'.

- pars1:

  Vector of model parameters:  
    
  `pars1[1]` corresponds to lambda^c (cladogenesis rate)  
  `pars1[2]` corresponds to mu (extinction rate)  
  `pars1[3]` corresponds to K (clade-level carrying capacity)  
  `pars1[4]` corresponds to gamma (immigration rate)  
  `pars1[5]` corresponds to lambda^a (anagenesis rate).

- pars2:

  Vector of settings:  
    
  `pars2[1]` corresponds to res, the maximum number of endemics or
  non-endemics for which the ODE system is solved; this must be much
  larger than the actual number for which the probability needs to be
  calculated.)  
  `pars2[2]` corresponds to M, size of the mainland pool, i.e the number
  of species that can potentially colonize the island.

- tvec:

  The times at which the probabilities need to be computed.

- initprobs:

  The initial probability distribution for the number of endemics and
  non-endemics; either this or initEI must be NULL.

- pb:

  Probability distribution in matrix format as output by
  [`DAISIE_probdist()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_probdist.md).

- island:

  Island data object. Can be in DAISIE list format (see
  Galapagos_datalist and DAISIE_data_prep for examples) or in table
  format (see Galapagos_datatable for an example).

- title:

  Title of the plot

- plot_lists_simulations:

  List with simulation output after parsing by
  `DAISIE_prepare_data_plotting`.

- plot_lists_simulations_MLE:

  List with simulation output after parsing by
  `DAISIE_prepare_data_plotting`, but obtained by simulating MLE output.

- kind_of_plot:

  Character vector stating how STT plot resulting from MLE based
  simulations should be plotted. Default is `"line"` for multiple
  individual lines. Can also be `"shade"` for the 5% quantile.

- resolution:

  numeric indicating resolution of plot. Should be \< 0.

- resol:

  numeric for resolution of summary stats calculation. Should be \> 1.

- removed_timepoints:

  Positive integer with number of first datapoints to be removed from
  rate plots (to prevent Inf)

- A:

  A numeric value for island area at a given point in time.

- Amin:

  A numeric value for minimum island area during the simulation.

- Amax:

  A numeric value for maximum island area during the simulation.

- peak:

  A numeric value specifying the peakiness (or shaprness) of the
  ontogeny curve. Higher values imply peakier ontogeny. This value is
  internally calculated by
  [`calc_peak()`](https://rsetienne.github.io/DAISIE/reference/calc_peak.md)
  given the area at the present and the `area_pars`.

- proptime:

  A numeric from 0 to 1. The proportion of time that has elapsed in the
  simulation, in relation to the total island age (NB: not the
  simulation time, but island age).

- proptime_max:

  A numeric from 0 to 1. The same as `proportional_peak_t`. Indicates,
  in proportion to the total island age when the ontogeny peak should
  occur (i.e. 0.5 means a peak halfway in time).

- current_area:

  A numeric with the current island area at present (i.e., at the end of
  the simulation).

- jitter:

  Numeric for
  [`optimizer()`](https://rsetienne.github.io/DDD/reference/optimizer.html).
  Jitters the parameters being optimized by the specified amount which
  should be very small, e.g. 1e-5. Jitter when
  `link{subplex}{subplex}()` produces incorrect output due to parameter
  transformation.

- num_cycles:

  The number of cycles the optimizer will go through. Default is 1.

- trait_pars:

  A named list containing diversification rates considering two trait
  states created by
  [`create_trait_pars`](https://rsetienne.github.io/DAISIE/reference/create_trait_pars.md):

  - \[1\]:A numeric with the per capita transition rate with state 1

  - \[2\]:A numeric with the per capita immigration rate with state 2

  - \[3\]:A numeric with the per capita extinction rate with state 2

  - \[4\]:A numeric with the per capita anagenesis rate with state 2

  - \[5\]:A numeric with the per capita cladogenesis rate with state 2

  - \[6\]:A numeric with the per capita transition rate with state 2

  - \[7\]:A numeric with the number of species with trait state 2 on
    mainland

- relaxed_par:

  A string determining which parameter is relaxed in a relaxed rate
  model.

- relaxed_rate_pars:

  A list of two numbers, element one is the distribution mean, element
  two is the distribution standard deviation (sd). Currently the
  distribution is the gamma distribution. The list can be created with
  `create_relaxed_rate_pars()`.

- brts:

  Numeric vector of branching times

- stac:

  Numeric of Endemicity status

- missnumspec:

  Numeric of missing species

- DAISIE_par:

  A numeric parameter to evaluate the integral of the function.

- DAISIE_dist_pars:

  A numeric vector of two elements, first is the mean and second the
  standard deviation of the distribution.

- abstolint:

  Numeric absolute tolerance of the integration

- reltolint:

  Numeric relative tolerance of the integration

- pick:

  Numeric determining which parameter is selected for the relaxed-rate
  model

- mean:

  Numeric mean of the distribution

- sd:

  Numeric standard deviation of the distribution

- clado_rate:

  Numeric rate of cladogenesis

- ext_rate:

  Numeric rate of extinction

- carr_cap:

  Numeric carrying capacity

- immig_rate:

  Numeric rate of immigration

- ana_rate:

  Numeric rate of anagenesis

- islands:

  Island datalist or simulated data in DAISIE datalist format. Can be a
  single island (empirical data) generated with DAISIE_dataprep or
  DAISIEprep. Can also be simulated data generated with DAISIE_sim
  function.

- sort_clade_sizes:

  Default sort_clade_sizes = T outputs clade sizes sorted in ascending
  order of number of species. sort_clade_sizes=F outputs clade sizes in
  the same order as they appear in the input datalist.

- equal_extinction:

  If FALSE the extinction rates of endemic and non-endemic species are
  different, otherwise they are set equal in optimization

- files_to_write:

  number of files to write simulations to file

- use_rcpp:

  If TRUE, use Rcpp implementation of DAISIE simulation core. Default is
  FALSE.

## Value

Nothing

## Note

This is an internal function, so it should be marked with `@noRd`. This
is not done, as this will disallow all functions to find the
documentation parameters
