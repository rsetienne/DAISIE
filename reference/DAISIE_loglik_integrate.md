# Integrates the loglikelihood of a single clade across a parameter weighted by a given distribution

Integrates the loglikelihood of a single clade across a parameter
weighted by a given distribution

## Usage

``` r
DAISIE_loglik_integrate(
  pars1,
  pars2,
  brts,
  stac,
  missnumspec,
  CS_version,
  methode,
  abstolint,
  reltolint,
  verbose
)
```

## Arguments

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

- brts:

  Numeric vector of branching times

- stac:

  Numeric of Endemicity status

- missnumspec:

  Numeric of missing species

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

- abstolint:

  Numeric absolute tolerance of the integration

- reltolint:

  Numeric relative tolerance of the integration

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

## Value

A loglikelihood value
