# The expectation and marginal distribution of the number of colonizations (lineages) under the DAISIE model

This function calculates expectation and marginal distribution of the
number of colonizations (lineages) for a given set of parameter values,
a given mainland species pool size and a given set of times

## Usage

``` r
DAISIE_numcol(pars1, pars2, tvec, initEI = NULL)
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

- tvec:

  The times at which the probabilities need to be computed.

- initEI:

  The initial values for the number of endemics and non-endemics. In
  [`DAISIE_probdist()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_probdist.md)
  or
  [`DAISIE_margprobdist()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_margprobdist.md)
  either this or initprobs must be NULL. In `DAISIE_numcol()` when it is
  NULL, it is assumed that the island is empty.

## Value

- out:

  A list of three vectors:  
    
  `expC` The expectation of the number of colonizations/lineages at the
  given times  
  `pC` The probability distribution of the number of colonizations
  (lineages) at the given times  

## References

Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium and
non-equilibrium dynamics simultaneously operate in the Galapagos
islands. Ecology Letters 18: 844-852.

## Author

Rampal S. Etienne

## Examples

``` r
### Compute the marginal probability distributions at t = 4 and t = 8, for a mainland
# pool size of 250 potential colonists and a vector of 5 parameters (cladogenesis,
# extinction, clade-level carrying capacity, immigration, anagenesis) starting from
# an empty island

numcol <- DAISIE_numcol(
   pars1 = c(0.3,0.35,Inf,0.75,0.012),
   pars2 = c(100,250),
   tvec = c(4,8),
   initEI = list(c(0,1),c(0,2),c(3,1))
   )
```
