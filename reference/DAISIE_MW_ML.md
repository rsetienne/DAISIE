# Maximization of the loglikelihood under the DAISIE model with clade-specific diversity-dependence and explicit dependencies on island area and isolation as hypothesized by MacArthur & Wilson

This function computes the maximum likelihood estimates of the
parameters of the relationships between parameters of the DAISIE model
(with clade-specific diversity-dependence) and island area and distance
of the island to the mainland for data from lineages colonizing several
islands/archipelagos. It also outputs the corresponding loglikelihood
that can be used in model comparisons.

A note on the sigmoidal functions used in distance_dep: For anagenesis
and cladogenesis, the functional relationship is k \* (d/d0)^x/(1 +
(d/d0)^x); for colonization the relationship is: k - k \* (d/d0)^x/(1 +
(d/d0)^x). The d0 parameter is the 11th parameter entered. In
'sigmoidal_col_ana', the 11th parameter is the d0 for colonization and
the 12th is the d0 for anagenesis.

## Usage

``` r
DAISIE_MW_ML(
  datalist,
  initparsopt,
  idparsopt,
  parsfix,
  idparsfix,
  res = 100,
  ddmodel = 11,
  cond = 0,
  island_ontogeny = NA,
  tol = c(1e-04, 1e-05, 1e-07),
  maxiter = 1000 * round((1.25)^length(idparsopt)),
  methode = "odeint::runge_kutta_cash_karp54",
  optimmethod = "simplex",
  CS_version = list(model = 1, function_to_optimize = "DAISIE"),
  verbose = 0,
  tolint = c(1e-16, 1e-10),
  distance_type = "continent",
  distance_dep = "power",
  parallel = "local",
  cpus = 3,
  num_cycles = 1
)
```

## Arguments

- datalist:

  Data object containing information on colonisation and branching times
  of species for several islands or archipelagos, as well as the area,
  isolation and age of each of the islands/archipelagos. See
  data(archipelagos41) for an example.

- initparsopt:

  The initial values of the parameters that must be optimized; they are
  all positive

- idparsopt:

  The ids of the parameters that must be optimized. The ids are defined
  as follows (see Valente et al 2020 Supplementary Tables 1 and 2 a
  better explanation of the models and parameters):  
    
  id = 1 corresponds to lambda^c0 (cladogenesis rate for unit area)  
  id = 2 corresponds to y (exponent of area for cladogenesis rate)  
  id = 3 corresponds to mu0 (extinction rate for unit area)  
  id = 4 corresponds to x (exponent of 1/area for extinction rate)  
  id = 5 corresponds to K0 (clade-level carrying capacity for unit
  area)  
  id = 6 corresponds to z (exponent of area for clade-level carrying
  capacity)  
  id = 7 corresponds to gamma0 (immigration rate for unit distance)  
  id = 8 corresponds to alpha (exponent of 1/distance for immigration
  rate)  
  id = 9 corresponds to lambda^a0 (anagenesis rate for unit distance)  
  id = 10 corresponds to beta (exponent of 1/distance for anagenesis
  rate)  
  id = 11 corresponds to d0 in models M15 to M19, and models with
  distance_dep = 'sigmoidal_col', 'sigmoidal_ana' or 'sigmoidal_clado';
  or d0 for colonisation (when specifying distance_dep =
  'sigmoidal_col_ana'  
  id = 12 corresponds to d0 for anagenesis when specifying distance_dep
  = 'sigmoidal_col_ana'  

- parsfix:

  The values of the parameters that should not be optimized

- idparsfix:

  The ids of the parameters that should not be optimized, e.g. c(1,3) if
  lambda^c and K should not be optimized.

- res:

  Sets the maximum number of species for which a probability must be
  computed, must be larger than the size of the largest clade

- ddmodel:

  Sets the model of diversity-dependence:  
    
  ddmodel = 0 : no diversity dependence  
  ddmodel = 1 : linear dependence in speciation rate  
  ddmodel = 11: linear dependence in speciation rate and in immigration
  rate  
  ddmodel = 2 : exponential dependence in speciation rate  
  ddmodel = 21: exponential dependence in speciation rate and in
  immigration rate  

- cond:

  cond = 0 : conditioning on island age  
  cond = 1 : conditioning on island age and non-extinction of the island
  biota  

- island_ontogeny:

  type of island ontonogeny. If NA, then constant ontogeny is assumed

- tol:

  Sets the tolerances in the optimization. Consists of:  
  reltolx = relative tolerance of parameter values in optimization  
  reltolf = relative tolerance of function value in optimization  
  abstolx = absolute tolerance of parameter values in optimization

- maxiter:

  Sets the maximum number of iterations in the optimization

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

  Method used in likelihood optimization. Default is "subplex" (see
  subplex package). Alternative is 'simplex' which was the method in
  previous versions.

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

  sets whether parameters and likelihood should be printed (1) or not
  (0)

- tolint:

  Vector of two elements containing the absolute and relative tolerance
  of the integration

- distance_type:

  Use 'continent' if the distance to the continent should be used, use
  'nearest_big' if the distance to the nearest big landmass should be
  used, and use 'biologically_realistic' if the distance should take
  into account some biologically realism, e.g. an average of the
  previous two if both are thought to contribute.

- distance_dep:

  Sets what type of distance dependence should be used. Default is a
  power law, denoted as 'power' (models M1-14 in Valente et al 2020).
  Alternatives are additive or interactive contributions of distance and
  area to the rate of cladogenesis ("area_additive_clado";
  "area_interactive_clado", "area_interactive_clado1" and
  "area_interactive_clado2"). Other alternatives are exponential
  relationship denoted by 'exp'; or sigmoids, either 'sigmoidal_col' for
  a sigmoid in the colonization, 'sigmoidal_ana' for sigmoidal
  anagenesis, 'sigmoidal_clado' for sigmoidal cladogenesis, and
  'sigmoidal_col_ana' for sigmoids in both colonization and
  anagenesis.  
  A key for the different options of distance_dep that should be
  specified to run the models from Valente et al 2020 (Supplementary
  Data Table 1 and 2) is given below:  
  \* M1 to M14 - 'power'  
  \* M15 - 'area_additive_clado'  
  \* M16 and M19 - 'area_interactive_clado' (M19 assumes y = 0)  
  \* M17 - 'area_interactive_clado1'  
  \* M18 - 'area_interactive_clado2'  
  \* M20 and M24 - sigmoidal_col'  
  \* M21, M25 and M28 - sigmoidal_ana'  
  \* M22 and M26 - 'sigmoidal_clado'  
  \* M23 and M27 - 'sigmoidal_col_ana'  

- parallel:

  Sets whether parallel computation should be used. Use 'no' if no
  parallel computing should be used, 'cluster' for parallel computing on
  a unix/linux cluster, and 'local' for parallel computation on a local
  machine.

- cpus:

  Number of cpus used in parallel computing. Default is 3. Will not have
  an effect if parallel = 'no'.

- num_cycles:

  The number of cycles the optimizer will go through. Default is 1.

## Value

The output is a dataframe containing estimated parameters and maximum
loglikelihood.

- lambda_c0:

  gives the maximum likelihood estimate of lambda^c, the rate of
  cladogenesis for unit area

- y:

  gives the maximum likelihood estimate of y, the exponent of area for
  the rate of cladogenesis

- mu0:

  gives the maximum likelihood estimate of mu0, the extinction rate

- x:

  gives the maximum likelihood estimate of x, the exponent of 1/area for
  the extinction rate

- K0:

  gives the maximum likelihood estimate of K0, the carrying-capacity for
  unit area

- z:

  gives the maximum likelihood estimate of z, the exponent of area for
  the carrying capacity

- gamma0:

  gives the maximum likelihood estimate of gamma0, the immigration rate
  for unit distance

- y:

  gives the maximum likelihood estimate of alpha, the exponent of
  1/distance for the rate of colonization

- lambda_a0:

  gives the maximum likelihood estimate of lambda^a0, the rate of
  anagenesis for unit distance

- beta:

  gives the maximum likelihood estimate of beta, the exponent of
  1/distance for the rate of anagenesis

- d0:

  gives the maximum likelihood estimate of d0, the parameter that
  controls the interactive or sigmoidal functions

- loglik:

  gives the maximum loglikelihood

- df:

  gives the number of estimated parameters, i.e. degrees of feedom

- conv:

  gives a message on convergence of optimization; conv = 0 means
  convergence

## References

Valente L, Phillimore AB, Melo M, Warren BH, Clegg SM, Havenstein K,
Tiedemann R, Illera JC, Thébaud C, Aschenbach T, Etienne RS. A simple
dynamic model explains island bird diversity worldwide (2020) Nature,
579, 92-96

## See also

[`DAISIE_ML_CS`](https://rsetienne.github.io/DAISIE/reference/DAISIE_ML.md),

## Author

Rampal S. Etienne & Luis Valente

## Examples

``` r
cat("
### Fit the M19 model as in Valente et al 2020, using the ML
parameters as starting values (see Supplementary Tables 1 and 2).

utils::data(archipelagos41)

DAISIE_MW_ML(
datalist= archipelagos41,
initparsopt =
c(0.040073803,  1.945656546,  0.150429656,
67.25643672,  0.293635061,  0.059096872,  0.382688527,
0.026510781),
idparsopt = c(1,3,4,7,8,9,10,11),
parsfix = c(0,Inf,0) ,
idparsfix = c(2,5,6),
res = 100,
ddmodel = 0,
methode = 'lsodes',
cpus = 4,
parallel = 'local',
optimmethod = 'subplex',
tol = c(1E-4, 1E-5, 1E-7),
distance_type = 'continent',
distance_dep = 'area_interactive_clado'
)
")
#> 
#> ### Fit the M19 model as in Valente et al 2020, using the ML
#> parameters as starting values (see Supplementary Tables 1 and 2).
#> 
#> utils::data(archipelagos41)
#> 
#> DAISIE_MW_ML(
#> datalist= archipelagos41,
#> initparsopt =
#> c(0.040073803,   1.945656546,    0.150429656,
#> 67.25643672, 0.293635061,    0.059096872,    0.382688527,
#> 0.026510781),
#> idparsopt = c(1,3,4,7,8,9,10,11),
#> parsfix = c(0,Inf,0) ,
#> idparsfix = c(2,5,6),
#> res = 100,
#> ddmodel = 0,
#> methode = 'lsodes',
#> cpus = 4,
#> parallel = 'local',
#> optimmethod = 'subplex',
#> tol = c(1E-4, 1E-5, 1E-7),
#> distance_type = 'continent',
#> distance_dep = 'area_interactive_clado'
#> )
```
