# Calculate anagenesis rate

Internal function. Calculates the anagenesis rate given the current
number of immigrant species and the per capita rate.

## Usage

``` r
get_ana_rate(laa, num_immigrants, island_spec = NULL, trait_pars = NULL)
```

## Arguments

- laa:

  A numeric with the per capita anagenesis rate.

- num_immigrants:

  A numeric with the current number of non-endemic species (a.k.a
  non-endemic species).

- island_spec:

  Matrix with current state of simulation containing number of species.

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

## See also

Other rate calculations:
[`DAISIE_plot_area()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_plot_area.md),
[`get_ext_rate()`](https://rsetienne.github.io/DAISIE/reference/get_ext_rate.md),
[`get_immig_rate()`](https://rsetienne.github.io/DAISIE/reference/get_immig_rate.md),
[`island_area()`](https://rsetienne.github.io/DAISIE/reference/island_area.md),
[`island_area_vector()`](https://rsetienne.github.io/DAISIE/reference/island_area_vector.md),
[`update_max_rates()`](https://rsetienne.github.io/DAISIE/reference/update_max_rates.md),
[`update_rates()`](https://rsetienne.github.io/DAISIE/reference/update_rates.md)

## Author

Pedro Neves, Joshua Lambert, Shu Xie
