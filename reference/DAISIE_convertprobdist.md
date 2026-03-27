# Converts the joint distribution of endemics and non-endemics under the DAISIE model to list format

This function converts the joint distribution of the number of endemics
and non-endemics from the matrix format of DAISIE_probdist to a list
format

## Usage

``` r
DAISIE_convertprobdist(pb)
```

## Arguments

- pb:

  Probability distribution in matrix format as output by
  [`DAISIE_probdist()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_probdist.md).

## Value

A list of length nrow(pb) containing matrices of square dimensions of
size sqrt(ncol - 1) containing the joint probabilities with endemics in
the rows and non-endemics in the columns. The last element of the list
is a vector a times at which the joint probability distribution is
evaluated.

## References

Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium and
non-equilibrium dynamics simultaneously operate in the Galapagos
islands. Ecology Letters 18: 844-852.

## Author

Rampal S. Etienne

## Examples

``` r
### Compute the probability distribution at t = 4 and t = 8, for a mainland pool
# size of 250 potential colonists and a vector of 5 parameters (cladogenesis, extinction,
# clade-level carrying capacity, immigration, anagenesis) starting from an empty
# island; store in list format

pb <- DAISIE_probdist(
   pars1 = c(0.3,0.35,Inf,0.75,0.012),
   pars2 = c(100,250),
   tvec = c(4,8),
   initEI = c(0,0),
   initprobs = NULL
   )
prob_dists <- DAISIE_convertprobdist(pb)
```
