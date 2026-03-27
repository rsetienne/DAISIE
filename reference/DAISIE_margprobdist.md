# The marginal distribution of endemics and non-endemics under the DAISIE model

This function calculates the marginal distribution of the number of
endemics and non-endemics and their sum for a given set of parameter
values, a given mainland species pool size and a given set of times

## Usage

``` r
DAISIE_margprobdist(
  pars1,
  pars2,
  tvec,
  initEI = c(0, 0),
  initprobs = NULL,
  pb = NULL
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

- tvec:

  The times at which the probabilities need to be computed.

- initEI:

  The initial values for the number of endemics and non-endemics. In
  [`DAISIE_probdist()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_probdist.md)
  or `DAISIE_margprobdist()` either this or initprobs must be NULL. In
  [`DAISIE_numcol()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_numcol.md)
  when it is NULL, it is assumed that the island is empty.

- initprobs:

  The initial probability distribution for the number of endemics and
  non-endemics; either this or initEI must be NULL.

- pb:

  Probability distribution in matrix format as output by
  [`DAISIE_probdist()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_probdist.md).

## Value

- out:

  A list of three vectors:  
    
  `pE` The probability distribution of the number of endemic species  
  `pI` The probability distribution of the number of non-endemic
  species  
  `pN` The probability distribution of the sum of the number of endemics
  and non-endemics

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

marg_prob_dists <- DAISIE_margprobdist(
   pars1 = c(0.3,0.35,Inf,0.75,0.012),
   pars2 = c(100,250),
   tvec = c(4,8),
   initEI = c(5,1),
   initprobs = NULL
   )
```
