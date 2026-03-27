# Computes the loglikelihood of the DAISIE model with island-wide diversity-dependence given data and a set of model parameters

Computes the loglikelihood of the DAISIE model given colonization and
branching times for lineages on an island, and a set of model parameters
for the DAISIE model with island-wide diversity-dependence

## Usage

``` r
DAISIE_loglik_IW(
  pars1,
  pars2,
  datalist,
  methode = "lsodes",
  abstolint = 1e-12,
  reltolint = 1e-10,
  verbose = FALSE
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
  `pars1[6]` is optional; it may contain M, the total number of species
  on the mainland  
    

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
  Only ddmodel = 11 is currently implemented  
    
  `pars2[3]` corresponds to cond = setting of conditioning  
    
  cond = 0 : conditioning on island age  
  cond = 1 : conditioning on island age and non-extinction of the island
  biota  
    
  `pars2[4]` Specifies whether intermediate output should be provided,
  because computation may take long. Default is 0, no output. A value of
  1 means the parameters and loglikelihood are printed. A value of 2
  means also intermediate progress during loglikelihood computation is
  shown.

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
    
  \* Non_endemic_MaxAge: 1  
  \* Endemic: 2  
  \* Endemic&Non_Endemic: 3  
  \* Non_endemic: 4  
  \* Endemic_MaxAge: 5  
    
  `$missing_species` - number of island species that were not sampled
  for particular clade (only applicable for endemic clades)  

- methode:

  Method of the ODE-solver. Supported Boost `ODEINT` solvers (steppers)
  are: `'odeint::runge_kutta_cash_karp54'`
  `'odeint::runge_kutta_fehlberg78'` \[default\]
  `'odeint::runge_kutta_dopri5'` `'odeint::bulirsch_stoer'`
  `'odeint::adams_bashforth_[1|2|3|4|5|6|7|8]} \code{'odeint::adams_bashforth_moulton_[1|2|3|4|5|6|7|8]`
  without `odeint::`-prefix,
  [`ode`](https://rdrr.io/pkg/deSolve/man/ode.html) method is assumed.

- abstolint:

  Absolute tolerance of the integration

- reltolint:

  Relative tolerance of the integration

- verbose:

  Logical controling if progress is printed to console.

## Value

The loglikelihood

## Details

The output is a loglikelihood value

## References

Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium and
non-equilibrium dynamics simultaneously operate in the Galapagos
islands. Ecology Letters 18: 844-852.

## See also

[`DAISIE_ML_IW`](https://rsetienne.github.io/DAISIE/reference/DAISIE_ML_IW.md),
[`DAISIE_loglik_CS`](https://rsetienne.github.io/DAISIE/reference/DAISIE_loglik_CS.md),
[`DAISIE_sim_cr`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim.md)

## Author

Rampal S. Etienne & Bart Haegeman
