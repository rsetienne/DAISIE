# The expected number of endemics and non-endemics under the DAISIE model with no diversity-dependence

This function calculates the expected number of endemics, non-endemics
and the sum of these for a given set of parameter values, a given
mainland species pool size and a given time, assuming no
diversity-dependence

## Usage

``` r
DAISIE_ExpEIN(tvec, pars, M, initEI = c(0, 0))
```

## Arguments

- tvec:

  The times at which the probabilities need to be computed.

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

- M:

  Numeric defining the size of mainland pool, i.e. the number of species
  that can potentially colonize the island.

- initEI:

  The initial values for the number of endemics and non-endemics. In
  [`DAISIE_probdist()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_probdist.md)
  or
  [`DAISIE_margprobdist()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_margprobdist.md)
  either this or initprobs must be NULL. In
  [`DAISIE_numcol()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_numcol.md)
  when it is NULL, it is assumed that the island is empty.

## Value

The output is a list with three elements:  
  
`ExpE` The number of endemic species  
`ExpI` The number of non-endemic species  
`ExpN` The sum of the number of endemics and non-endemics

## References

Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium and
non-equilibrium dynamics simultaneously operate in the Galapagos
islands. Ecology Letters 18: 844-852.

## Author

Rampal S. Etienne

## Examples

``` r
### Compute the expected values at t = 4, for a mainland pool size of 1000 potential
# colonists and a vector of 5 parameters (cladogenesis, extinction, clade-level carrying
# capacity, immigration, anagenesis)

DAISIE_ExpEIN(
   tvec = c(2,4),
   pars = c(0.5,0.1,Inf,0.01,0.4),
   M = 1000
   )
#> $ExpE
#> [1] 21.89456 88.46714
#> 
#> $ExpI
#> [1] 8.587570 9.726758
#> 
#> $ExpN
#> [1] 30.48213 98.19390
#> 
```
