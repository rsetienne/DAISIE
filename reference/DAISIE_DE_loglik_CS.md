# Computes the loglikelihood of the DAISIE_DE model given data and a set of model parameters

Computes the loglikelihood of the DAISIE_DE model given colonization and
branching times for lineages on an island, and a set of model
parameters. The output is a loglikelihood value

## Usage

``` r
DAISIE_DE_loglik_CS(
  pars1,
  pars2,
  datalist,
  methode = "odeint::runge_kutta_cash_karp54",
  abstolint = 1e-15,
  reltolint = 1e-15,
  equal_extinction = TRUE,
  sampling = "n"
)
```

## Arguments

- pars1:

  Contains the model parameters:  
    
  `pars1[1]` corresponds to lambda^c (cladogenesis rate)  
  `pars1[2]` corresponds to mu (extinction rate of endemic species)  
  `pars1[3]` corresponds to mu2 (extinction rate of non-endemic
  species)  
  `pars1[4]` corresponds to gamma (immigration rate)  
  `pars1[5]` corresponds to lambda^a (anagenesis rate)  

- pars2:

  Contains the model settings  
    
  `pars2[1]` irrelevant for DAISIE_DE  
  `pars2[2]` irrelevant for DAISIE_DE  
  `pars2[3]` corresponds to cond = setting of conditioning  
    
  cond = 0 : conditioning on island age  
  cond = 1 : conditioning on island age and non-extinction of the island
  biota  
    
  cond \> 1 : conditioning on island age and having at least cond
  colonizations on the island  
    
  `pars2[4]` sets the level of verbosity. When equal to 0, no output is
  generated. At higher values (1 or 2) more output will be generated.

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
  `$branching_times` - island age and stem age of the population/species
  in the case of Non-endemic, Non-endemic_MaxAge and Endemic anagenetic
  species. For cladogenetic species these should be island age and
  branching times of the radiation including the stem age of the
  radiation.  
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
  `$type1or2` - whether the colonist belongs to type 1 or type 2.
  Currently not implemented for DAISIE_DE  

- methode:

  Method of the ODE-solver. See package deSolve for details. Default is
  "odeint::runge_kutta_cask_karp54"

- abstolint:

  Absolute tolerance of the integration

- reltolint:

  Relative tolerance of the integration

- equal_extinction:

  If FALSE the extinction rates of endemic and non-endemic species are
  different, otherwise they are set equal in optimization

- sampling:

  Determines whether n-sampling or rho-sampling should be used when
  function_to_optimize = 'DAISIE_DE'.

## Value

The loglikelihood

## References

O.N. Dehayem et al. 2026. Preprint.

## See also

[`DAISIE_ML`](https://rsetienne.github.io/DAISIE/reference/DAISIE_ML.md)

## Author

Rampal S. Etienne & Bart Haegeman
