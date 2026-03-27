# Maximization of the loglikelihood under the DAISIE model with clade-specific diversity-dependence

This function computes the maximum likelihood estimates of the
parameters of the DAISIE model with clade-specific diversity-dependence
for data from lineages colonizing an island. It also outputs the
corresponding loglikelihood that can be used in model comparisons. The
result of sort(c(idparsopt, idparsfix, idparsnoshift)) should be
identical to c(1:10). If not, an error is reported that the input is
incoherent. The same happens when the length of initparsopt is different
from the length of idparsopt, and the length of parsfix is different
from the length of idparsfix.  
Including the 11th parameter (p_f) in either idparsopt or idparsfix (and
therefore initparsopt or parsfix) is optional. If this parameter is not
specified, then the information in the data is used, otherwise the
information in the data is overruled.

## Usage

``` r
DAISIE_ML_CS(
  datalist,
  datatype = "single",
  initparsopt,
  idparsopt,
  parsfix,
  idparsfix,
  idparsnoshift = 6:10,
  idparsmat = NULL,
  res = 100,
  ddmodel = 0,
  cond = 0,
  island_ontogeny = NA,
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

- island_ontogeny:

  In
  [`DAISIE_sim_time_dep()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_time_dep.md),
  `DAISIE_ML_CS` and plotting a string describing the type of island
  ontogeny. Can be `"const"`, `"beta"` for a beta function describing
  area through time.  
  In all other functions a numeric describing the type of island
  ontogeny. Can be `0` for constant, `1` for a beta function describing
  area through time. In ML functions `island_ontogeny = NA` assumes
  constant ontogeny. Time dependent estimation is not yet available as
  development is still ongoing. Will return an error if called in that
  case.

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
  The default method overall is `"lsodes"` for `DAISIE_ML_CS()` and
  `"ode45"` from [`ode()`](https://rdrr.io/pkg/deSolve/man/ode.html) for
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

- lambda_c2:

  gives the maximum likelihood estimate of lambda^c2, the rate of
  cladogenesis for the optional second group of species

- mu2:

  gives the maximum likelihood estimate of mu2, the extinction rate for
  the optional second group of species

- K2:

  gives the maximum likelihood estimate of K2, the carrying-capacity for
  the optional second group of species

- gamma2:

  gives the maximum likelihood estimate of gamma2, the immigration rate
  for the optional second group of species

- lambda_a2:

  gives the maximum likelihood estimate of lambda^a2, the rate of
  anagenesis for the optional second group of species

- loglik:

  gives the maximum loglikelihood

- df:

  gives the number of estimated parameters, i.e. degrees of feedom

- conv:

  gives a message on convergence of optimization; conv = 0 means
  convergence

## References

Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium and
non-equilibrium dynamics simultaneously operate in the Galapagos
islands. Ecology Letters 18: 844-852. \<doi:10.1111/ele.12461\>.

## See also

[`DAISIE_loglik_all`](https://rsetienne.github.io/DAISIE/reference/DAISIE_loglik_CS.md),
[`DAISIE_sim_cr`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim.md),
[`DAISIE_sim_time_dep`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_time_dep.md),
[`DAISIE_sim_cr_shift`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_cr_shift.md)

## Author

Rampal S. Etienne

## Examples

``` r
cat("
### When all species have the same rates, and we want to optimize all 5 parameters,
# we use:

utils::data(Galapagos_datalist)
DAISIE_ML(
   datalist = Galapagos_datalist,
   initparsopt = c(2.5,2.7,20,0.009,1.01),
   ddmodel = 11,
   idparsopt = 1:5,
   parsfix = NULL,
   idparsfix = NULL
)

### When all species have the same rates, and we want to optimize all parameters
# except K (which we set equal to Inf), we use:

utils::data(Galapagos_datalist)
DAISIE_ML(
   datalist = Galapagos_datalist,
   initparsopt = c(2.5,2.7,0.009,1.01),
   idparsopt = c(1,2,4,5),
   parsfix = Inf,
   idparsfix = 3
   )

### When all species have the same rates except that the finches have a different
# rate of cladogenesis, and we want to optimize all parameters except K (which we
# set equal to Inf), fixing the proportion of finch-type species at 0.163, we use:

utils::data(Galapagos_datalist_2types)
DAISIE_ML(
   datalist = Galapagos_datalist_2types,
   initparsopt = c(0.38,0.55,0.004,1.1,2.28),
   idparsopt = c(1,2,4,5,6),
   parsfix = c(Inf,Inf,0.163),
   idparsfix = c(3,8,11),
   idparsnoshift = c(7,9,10)
   )

### When all species have the same rates except that the finches have a different
# rate of cladogenesis, extinction and a different K, and we want to optimize all
# parameters, fixing the proportion of finch-type species at 0.163, we use:

utils::data(Galapagos_datalist_2types)
DAISIE_ML(
   datalist = Galapagos_datalist_2types,
   ddmodel = 11,
   initparsopt = c(0.19,0.09,0.002,0.87,20,8.9,15),
   idparsopt = c(1,2,4,5,6,7,8),
   parsfix = c(Inf,0.163),
   idparsfix = c(3,11),
   idparsnoshift = c(9,10)
   )


### When all species have the same rates except that the finches have a different
# rate of extinction, and we want to optimize all parameters except K (which we
# set equal to Inf), and we also# want to estimate the fraction of finch species
# in the mainland pool. we use:

utils::data(Galapagos_datalist_2types)
DAISIE_ML(
   datalist = Galapagos_datalist_2types,
   initparsopt = c(2.48,2.7,0.009,1.01,2.25,0.163),
   idparsopt = c(1,2,4,5,7,11),
   parsfix = c(Inf,Inf),
   idparsfix = c(3,8),
   idparsnoshift = c(6,9,10)
   )

### When we have two islands with the same rates except for immigration and anagenesis rate,
# and we want to optimize all parameters, we use:

utils::data(Galapagos_datalist)
DAISIE_ML(
   datalist = list(Galapagos_datalist,Galapagos_datalist),
   datatype = 'multiple',
   initparsopt = c(2.5,2.7,20,0.009,1.01,0.009,1.01),
   idparsmat = rbind(1:5,c(1:3,6,7)),
   idparsopt = 1:7,
   parsfix = NULL,
   idparsfix = NULL
)

### When we consider the four Macaronesia archipelagoes and set all parameters the same
# except for rates of cladogenesis, extinction and immigration for Canary Islands,
# rate of cladogenesis is fixed to 0 for the other archipelagoes,
# diversity-dependence is assumed to be absent
# and we want to optimize all parameters, we use:

utils::data(Macaronesia_datalist)
DAISIE_ML(
   datalist = Macaronesia_datalist,
   datatype = 'multiple',
   initparsopt = c(1.053151832,0.052148979,0.512939011,0.133766934,0.152763179),
   idparsmat = rbind(1:5,c(6,2,3,7,5),1:5,1:5),
   idparsopt = c(2,4,5,6,7),
   parsfix = c(0,Inf),
   idparsfix = c(1,3)
)

")
#> 
#> ### When all species have the same rates, and we want to optimize all 5 parameters,
#> # we use:
#> 
#> utils::data(Galapagos_datalist)
#> DAISIE_ML(
#>    datalist = Galapagos_datalist,
#>    initparsopt = c(2.5,2.7,20,0.009,1.01),
#>    ddmodel = 11,
#>    idparsopt = 1:5,
#>    parsfix = NULL,
#>    idparsfix = NULL
#> )
#> 
#> ### When all species have the same rates, and we want to optimize all parameters
#> # except K (which we set equal to Inf), we use:
#> 
#> utils::data(Galapagos_datalist)
#> DAISIE_ML(
#>    datalist = Galapagos_datalist,
#>    initparsopt = c(2.5,2.7,0.009,1.01),
#>    idparsopt = c(1,2,4,5),
#>    parsfix = Inf,
#>    idparsfix = 3
#>    )
#> 
#> ### When all species have the same rates except that the finches have a different
#> # rate of cladogenesis, and we want to optimize all parameters except K (which we
#> # set equal to Inf), fixing the proportion of finch-type species at 0.163, we use:
#> 
#> utils::data(Galapagos_datalist_2types)
#> DAISIE_ML(
#>    datalist = Galapagos_datalist_2types,
#>    initparsopt = c(0.38,0.55,0.004,1.1,2.28),
#>    idparsopt = c(1,2,4,5,6),
#>    parsfix = c(Inf,Inf,0.163),
#>    idparsfix = c(3,8,11),
#>    idparsnoshift = c(7,9,10)
#>    )
#> 
#> ### When all species have the same rates except that the finches have a different
#> # rate of cladogenesis, extinction and a different K, and we want to optimize all
#> # parameters, fixing the proportion of finch-type species at 0.163, we use:
#> 
#> utils::data(Galapagos_datalist_2types)
#> DAISIE_ML(
#>    datalist = Galapagos_datalist_2types,
#>    ddmodel = 11,
#>    initparsopt = c(0.19,0.09,0.002,0.87,20,8.9,15),
#>    idparsopt = c(1,2,4,5,6,7,8),
#>    parsfix = c(Inf,0.163),
#>    idparsfix = c(3,11),
#>    idparsnoshift = c(9,10)
#>    )
#> 
#> 
#> ### When all species have the same rates except that the finches have a different
#> # rate of extinction, and we want to optimize all parameters except K (which we
#> # set equal to Inf), and we also# want to estimate the fraction of finch species
#> # in the mainland pool. we use:
#> 
#> utils::data(Galapagos_datalist_2types)
#> DAISIE_ML(
#>    datalist = Galapagos_datalist_2types,
#>    initparsopt = c(2.48,2.7,0.009,1.01,2.25,0.163),
#>    idparsopt = c(1,2,4,5,7,11),
#>    parsfix = c(Inf,Inf),
#>    idparsfix = c(3,8),
#>    idparsnoshift = c(6,9,10)
#>    )
#> 
#> ### When we have two islands with the same rates except for immigration and anagenesis rate,
#> # and we want to optimize all parameters, we use:
#> 
#> utils::data(Galapagos_datalist)
#> DAISIE_ML(
#>    datalist = list(Galapagos_datalist,Galapagos_datalist),
#>    datatype = 'multiple',
#>    initparsopt = c(2.5,2.7,20,0.009,1.01,0.009,1.01),
#>    idparsmat = rbind(1:5,c(1:3,6,7)),
#>    idparsopt = 1:7,
#>    parsfix = NULL,
#>    idparsfix = NULL
#> )
#> 
#> ### When we consider the four Macaronesia archipelagoes and set all parameters the same
#> # except for rates of cladogenesis, extinction and immigration for Canary Islands,
#> # rate of cladogenesis is fixed to 0 for the other archipelagoes,
#> # diversity-dependence is assumed to be absent
#> # and we want to optimize all parameters, we use:
#> 
#> utils::data(Macaronesia_datalist)
#> DAISIE_ML(
#>    datalist = Macaronesia_datalist,
#>    datatype = 'multiple',
#>    initparsopt = c(1.053151832,0.052148979,0.512939011,0.133766934,0.152763179),
#>    idparsmat = rbind(1:5,c(6,2,3,7,5),1:5,1:5),
#>    idparsopt = c(2,4,5,6,7),
#>    parsfix = c(0,Inf),
#>    idparsfix = c(1,3)
#> )
#> 
```
