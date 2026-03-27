# Computes the loglikelihood of the DAISIE model with clade-specific diversity-dependence given data and a set of model parameters that may shift at some time

Computes the loglikelihood of the DAISIE model with clade-specific
diversity-dependence given colonization and branching times for lineages
on an island, and a set of model parameters that may shift at some time

## Usage

``` r
DAISIE_SR_loglik_CS(
  pars1,
  pars2,
  datalist,
  methode = "odeint::runge_kutta_cash_karp54",
  CS_version = list(model = 1, function_to_optimize = "DAISIE"),
  abstolint = 1e-16,
  reltolint = 1e-10
)
```

## Arguments

- pars1:

  Contains the model parameters:  
    
  `pars1[1]` corresponds to lambda^c (cladogenesis rate)  
  `pars1[2]` corresponds to mu (extinction rate)  
  `pars1[3]` corresponds to K (clade-level carrying capacity)  
  `pars1[4]` corresponds to gamma (immigration rate)  
  `pars1[5]` corresponds to lambda^a (anagenesis rate)  
  `pars1[6]` corresponds to lambda^c (cladogenesis rate) after the
  shift  
  `pars1[7]` corresponds to mu (extinction rate) after the shift  
  `pars1[8]` corresponds to K (clade-level carrying capacity) after the
  shift  
  `pars1[9]` corresponds to gamma (immigration rate) after the shift  
  `pars1[10]` corresponds to lambda^a (anagenesis rate) after the
  shift  
  `pars1[11]` corresponds to the time of shift  

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
    
  `pars2[3]` corresponds to cond = setting of conditioning  
    
  cond = 0 : conditioning on island age  
  cond = 1 : conditioning on island age and non-extinction of the island
  biota  
  cond \> 1 : conditioning on island age and having at least cond
  colonizations on the island  
    
  `pars2[4]` sets whether parameters and likelihood should be
  printed (1) or not (0)

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
  The remaining elements of the list each contains information on a
  single colonist lineage on the island and has 5 components:  
    
  `$colonist_name` - the name of the species or clade that colonized the
  island  
  `$branching_times` - island age and stem age of the population/species
  in the case of Non-endemic, Non-endemic_MaxAge and Endemic anagenetic
  species. For cladogenetic species these should be island age and
  branching times of the radiation including the stem age of the
  radiation.  
  `$stac` - the status of the colonist  
    
  - Non_endemic_MaxAge: 1  
  \* Endemic: 2  
  - Endemic&Non_Endemic: 3  
  - Non_endemic: 4  
  - Endemic_MaxAge: 5  
    
  `$missing_species` - number of island species that were not sampled
  for particular clade (only applicable for endemic clades)  

- methode:

  Method of the ODE-solver. See package deSolve for details. Default is
  "lsodes"

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

  Absolute tolerance of the integration

- reltolint:

  Relative tolerance of the integration

## Value

The loglikelihood

## Details

The output is a loglikelihood value

## References

Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium and
non-equilibrium dynamics simultaneously operate in the Galapagos
islands. Ecology Letters 18: 844-852.

## See also

[`DAISIE_ML`](https://rsetienne.github.io/DAISIE/reference/DAISIE_ML.md),
[`DAISIE_sim_cr`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim.md)

## Author

Rampal S. Etienne & Bart Haegeman

## Examples

``` r
utils::data(Galapagos_datalist_2types)
pars1 = c(0.195442017,0.087959583,Inf,0.002247364,0.873605049,
          3755.202241,8.909285094,14.99999923,0.002247364,0.873605049,0.163)
pars2 = c(100,11,0,1)
DAISIE_loglik_all(pars1,pars2,Galapagos_datalist_2types)
#> High lambda detected; approximation used.
#> High lambda detected; approximation used.
#> [1] -61.70372
```
