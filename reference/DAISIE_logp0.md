# Computes the log probability of no species present under the DAISIE model with clade-specific diversity-dependence

Computes the log probability of no species present under the DAISIE
model with clade-specific diversity-dependence. The output is a log
value.

## Usage

``` r
DAISIE_logp0(
  pars1,
  pars2,
  island_age,
  methode = "odeint::runge_kutta_cash_karp54",
  CS_version = list(model = 1, function_to_optimize = "DAISIE"),
  abstolint = 1e-16,
  reltolint = 1e-10
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

  Contains the model settings  
    
  `pars2[1]` corresponds to lx = length of ODE variable x  
  `pars2[2]` corresponds to ddmodel = diversity-dependent model, model
  of diversity-dependence, which can be one of  
    
  ddmodel = 0 : no diversity dependence  
  ddmodel = 1 : linear dependence in speciation rate  
  ddmodel = 11: linear dependence in speciation rate and in immigration
  rate  
  ddmodel = 2 : exponential dependence in speciation rate  
  ddmodel = 21: exponential dependence in speciation rate and in
  immigration rate  
    

- island_age:

  the island age  

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

- abstolint:

  Numeric absolute tolerance of the integration

- reltolint:

  Numeric relative tolerance of the integration

## Value

The logarithm of the probability

## Author

Rampal S. Etienne & Bart Haegeman
