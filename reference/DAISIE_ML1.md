# Computes MLE for single type species under a clade specific scenario

Computes MLE for single type species under a clade specific scenario

## Usage

``` r
DAISIE_ML1(
  datalist,
  initparsopt,
  idparsopt,
  parsfix,
  idparsfix,
  idparsnoshift = 6:10,
  res = 100,
  ddmodel = 0,
  cond = 0,
  eqmodel = 0,
  x_E = 0.95,
  x_I = 0.98,
  tol = c(1e-04, 1e-05, 1e-07),
  maxiter = 1000 * round((1.25)^length(idparsopt)),
  methode = "odeint::runge_kutta_cash_karp54",
  optimmethod = "simplex",
  CS_version = list(model = 1, function_to_optimize = "DAISIE"),
  verbose = 0,
  tolint = c(1e-16, 1e-10),
  island_ontogeny = NA,
  jitter = 0,
  num_cycles = 1,
  equal_extinction = TRUE
)
```

## Arguments

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

- parsfix:

  The values of the parameters that should not be optimized.

- idparsfix:

  The ids of the parameters that should not be optimized, e.g. c(1,3) if
  lambda^c and K should not be optimized.

- idparsnoshift:

  For datatype = 'single' only: The ids of the parameters that should
  not be different between two groups of species; This can only apply to
  ids 6:10, e.g. idparsnoshift = c(6,7) means that lambda^c and mu have
  the same values for both groups.

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

- tolint:

  Vector of two elements containing the absolute and relative tolerance
  of the integration.

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

- jitter:

  Numeric for
  [`optimizer()`](https://rsetienne.github.io/DDD/reference/optimizer.html).
  Jitters the parameters being optimized by the specified amount which
  should be very small, e.g. 1e-5. Jitter when
  `link{subplex}{subplex}()` produces incorrect output due to parameter
  transformation.

- num_cycles:

  The number of cycles the optimizer will go through. Default is 1.

- equal_extinction:

  If FALSE the extinction rates of endemic and non-endemic species are
  different, otherwise they are set equal in optimization

## Value

The output is a dataframe containing estimated parameters and maximum
loglikelihood.

- lambda_c:

  gives the maximum likelihood estimate of lambda^c, the rate of
  cladogenesis

- mu:

  gives the maximum likelihood estimate of mu, the extinction rate

- K:

  gives the maximum likelihood estimate of K, the carrying-capacity

- gamma:

  gives the maximum likelihood estimate of gamma, the immigration rate

- lambda_a:

  gives the maximum likelihood estimate of lambda^a, the rate of
  anagenesis

- loglik:

  gives the maximum loglikelihood

- df:

  gives the number of estimated parameters, i.e. degrees of feedom

- conv:

  gives a message on convergence of optimization; conv = 0 means
  convergence
