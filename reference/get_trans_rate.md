# Calculate transition rate

Internal function. Calculates the transition rate given the current
number of immigrant species and the per capita rate.

## Usage

``` r
get_trans_rate(trait_pars, island_spec)
```

## Arguments

- trait_pars:

  A named list containing diversification rates considering two trait
  states created by
  [`create_trait_pars`](https://rsetienne.github.io/DAISIE/reference/create_trait_pars.md):

  - \[1\]:A numeric with the per capita transition rate with state1

  - \[2\]:A numeric with the per capita immigration rate with state2

  - \[3\]:A numeric with the per capita extinction rate with state2

  - \[4\]:A numeric with the per capita anagenesis rate with state2

  - \[5\]:A numeric with the per capita cladogenesis rate with state2

  - \[6\]:A numeric with the per capita transition rate with state2

  - \[7\]:A numeric with the number of species with trait state 2 on
    mainland

- island_spec:

  Matrix with current state of simulation containing number of species.
