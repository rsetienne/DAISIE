# Calculate cladogenesis rate

Internal function. Calculates the cladogenesis rate given the current
number of species in the system, the carrying capacity and the per
capita cladogenesis rate

## Usage

``` r
get_clado_rate(
  lac,
  hyper_pars,
  num_spec,
  K,
  A,
  trait_pars = NULL,
  island_spec = NULL
)
```

## Arguments

- lac:

  A numeric with the per capita cladogenesis rate.

- hyper_pars:

  A named list of numeric hyperparameters for the rate calculations as
  returned by
  [`create_hyper_pars()`](https://rsetienne.github.io/DAISIE/reference/create_hyper_pars.md):

  - \[1\]: is d the scaling parameter for exponent for calculating
    cladogenesis rate

  - \[2\]: is x the exponent for calculating extinction rate

- num_spec:

  A numeric with the current number of species.

- K:

  A numeric with carrying capacity.

- A:

  A numeric value for island area at a given point in time.

- trait_pars:

  A named list containing diversification rates considering two trait
  states created by
  [`create_trait_pars`](https://rsetienne.github.io/DAISIE/reference/create_trait_pars.md):

  - \[1\]:A numeric with the per capita transition rate with state 1

  - \[2\]:A numeric with the per capita immigration rate with state 2

  - \[3\]:A numeric with the per capita extinction rate with state 2

  - \[4\]:A numeric with the per capita anagenesis rate with state 2

  - \[5\]:A numeric with the per capita cladogenesis rate with state 2

  - \[6\]:A numeric with the per capita transition rate with state 2

  - \[7\]:A numeric with the number of species with trait state 2 on
    mainland

- island_spec:

  Matrix with current state of simulation containing number of species.

## Author

Pedro Neves, Joshua Lambert, Shu Xie
