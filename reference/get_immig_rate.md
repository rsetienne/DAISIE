# Calculate immigration rate

Internal function. Calculates the immigration rate given the current
number of species in the system, the carrying capacity

## Usage

``` r
get_immig_rate(
  gam,
  A = 1,
  num_spec,
  K,
  mainland_n,
  trait_pars = NULL,
  island_spec = NULL
)
```

## Arguments

- gam:

  A numeric with the per capita immigration rate.

- A:

  A numeric value for island area at a given point in time.

- num_spec:

  A numeric with the current number of species.

- K:

  A numeric with carrying capacity.

- mainland_n:

  A numeric stating the number of mainland species, that is the number
  of species that can potentially colonize the island. If using a
  clade-specific diversity dependence, this value is set to 1. If using
  an island-wide diversity dependence, this value is set to the number
  of mainland species.

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

## References

Valente, Luis M., Rampal S. Etienne, and Albert B. Phillimore. "The
effects of island ontogeny on species diversity and phylogeny."
Proceedings of the Royal Society of London B: Biological Sciences
281.1784 (2014): 20133227.

## See also

Other rate calculations:
[`DAISIE_plot_area()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_plot_area.md),
[`get_ana_rate()`](https://rsetienne.github.io/DAISIE/reference/get_ana_rate.md),
[`get_ext_rate()`](https://rsetienne.github.io/DAISIE/reference/get_ext_rate.md),
[`island_area()`](https://rsetienne.github.io/DAISIE/reference/island_area.md),
[`island_area_vector()`](https://rsetienne.github.io/DAISIE/reference/island_area_vector.md),
[`update_max_rates()`](https://rsetienne.github.io/DAISIE/reference/update_max_rates.md),
[`update_rates()`](https://rsetienne.github.io/DAISIE/reference/update_rates.md)

## Author

Pedro Neves, Joshua Lambert
