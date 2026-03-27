# Calculate extinction rate

Calculate extinction rate

## Usage

``` r
get_ext_rate(
  mu,
  hyper_pars,
  extcutoff = 1000,
  num_spec,
  A = 1,
  trait_pars = NULL,
  island_spec = NULL
)
```

## Arguments

- mu:

  A numeric with the per capita extinction rate.

- hyper_pars:

  A named list of numeric hyperparameters for the rate calculations as
  returned by
  [`create_hyper_pars()`](https://rsetienne.github.io/DAISIE/reference/create_hyper_pars.md):

  - \[1\]: is d the scaling parameter for exponent for calculating
    cladogenesis rate

  - \[2\]: is x the exponent for calculating extinction rate

- extcutoff:

  A numeric with the cutoff for the the maximum extinction rate
  preventing it from being too large and slowing down simulation.

- num_spec:

  A numeric with the current number of species.

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

## Value

A numeric, with the extinction rate given the base extinction rate, if
present, the hyperparemeter `x`, A(t) if time-dependent and traits if
running a traint model.

## References

Valente, Luis M., Rampal S. Etienne, and Albert B. Phillimore. "The
effects of island ontogeny on species diversity and phylogeny."
Proceedings of the Royal Society of London B: Biological Sciences
281.1784 (2014): 20133227.

## See also

Other rate calculations:
[`DAISIE_plot_area()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_plot_area.md),
[`get_ana_rate()`](https://rsetienne.github.io/DAISIE/reference/get_ana_rate.md),
[`get_immig_rate()`](https://rsetienne.github.io/DAISIE/reference/get_immig_rate.md),
[`island_area()`](https://rsetienne.github.io/DAISIE/reference/island_area.md),
[`island_area_vector()`](https://rsetienne.github.io/DAISIE/reference/island_area_vector.md),
[`update_max_rates()`](https://rsetienne.github.io/DAISIE/reference/update_max_rates.md),
[`update_rates()`](https://rsetienne.github.io/DAISIE/reference/update_rates.md)

## Author

Pedro Neves, Joshua Lambert, Shu Xie
