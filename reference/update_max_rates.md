# Calculates the maximum rates for a Gillespie simulation

Internal function that updates the all the max rates at time t.

## Usage

``` r
update_max_rates(
  gam,
  laa,
  lac,
  mu,
  hyper_pars = NULL,
  extcutoff,
  K,
  num_spec,
  num_immigrants,
  mainland_n,
  Amin,
  Amax
)
```

## Arguments

- gam:

  A numeric with the per capita immigration rate.

- laa:

  A numeric with the per capita anagenesis rate.

- lac:

  A numeric with the per capita cladogenesis rate.

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

- K:

  A numeric with carrying capacity.

- num_spec:

  A numeric with the current number of species.

- num_immigrants:

  A numeric with the current number of non-endemic species (a.k.a
  non-endemic species).

- mainland_n:

  A numeric stating the number of mainland species, that is the number
  of species that can potentially colonize the island. If using a
  clade-specific diversity dependence, this value is set to 1. If using
  an island-wide diversity dependence, this value is set to the number
  of mainland species.

- Amin:

  A numeric value for minimum island area during the simulation.

- Amax:

  A numeric value for maximum island area during the simulation.

## Value

a named list with the updated effective rates.

## See also

[`update_rates()`](https://rsetienne.github.io/DAISIE/reference/update_rates.md)

Other rate calculations:
[`DAISIE_plot_area()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_plot_area.md),
[`get_ana_rate()`](https://rsetienne.github.io/DAISIE/reference/get_ana_rate.md),
[`get_ext_rate()`](https://rsetienne.github.io/DAISIE/reference/get_ext_rate.md),
[`get_immig_rate()`](https://rsetienne.github.io/DAISIE/reference/get_immig_rate.md),
[`island_area()`](https://rsetienne.github.io/DAISIE/reference/island_area.md),
[`island_area_vector()`](https://rsetienne.github.io/DAISIE/reference/island_area_vector.md),
[`update_rates()`](https://rsetienne.github.io/DAISIE/reference/update_rates.md)
