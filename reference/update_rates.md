# Calculates algorithm rates

Internal function that updates the all the rates and max extinction
horizon at time t.

## Usage

``` r
update_rates(
  timeval,
  total_time,
  gam,
  laa,
  lac,
  mu,
  hyper_pars = hyper_pars,
  area_pars = NULL,
  peak = NULL,
  island_ontogeny = NULL,
  sea_level = NULL,
  extcutoff,
  K,
  num_spec,
  num_immigrants,
  mainland_n,
  trait_pars = NULL,
  island_spec = NULL
)
```

## Arguments

- timeval:

  Numeric defining current time of simulation.

- total_time:

  Numeric defining the length of the simulation in time units.

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

- area_pars:

  A named list containing area and sea level parameters as created by
  [`create_area_pars()`](https://rsetienne.github.io/DAISIE/reference/create_area_pars.md):

  - \[1\]: maximum area

  - \[2\]: current area

  - \[3\]: value from 0 to 1 indicating where in the island's history
    the peak area is achieved

  - \[4\]: total island age

  - \[5\]: amplitude of area fluctuation from sea level

  - \[6\]: frequency of sine wave of area change from sea level

  - \[7\]: angle of the slope of the island

- peak:

  A numeric value specifying the peakiness (or shaprness) of the
  ontogeny curve. Higher values imply peakier ontogeny. This value is
  internally calculated by
  [`calc_peak()`](https://rsetienne.github.io/DAISIE/reference/calc_peak.md)
  given the area at the present and the `area_pars`.

- island_ontogeny:

  In
  [`DAISIE_sim_time_dep()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_time_dep.md),
  [`DAISIE_ML_CS`](https://rsetienne.github.io/DAISIE/reference/DAISIE_ML.md)
  and plotting a string describing the type of island ontogeny. Can be
  `"const"`, `"beta"` for a beta function describing area through
  time.  
  In all other functions a numeric describing the type of island
  ontogeny. Can be `0` for constant, `1` for a beta function describing
  area through time. In ML functions `island_ontogeny = NA` assumes
  constant ontogeny. Time dependent estimation is not yet available as
  development is still ongoing. Will return an error if called in that
  case.

- sea_level:

  In
  [`DAISIE_sim_time_dep()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_time_dep.md)
  and plotting a string describing the type of sea level. Can be
  `"const"` or `"sine"` for a sine function describing area through
  time.  
  In all other functions a numeric describing the type of sea level. Can
  be `0` for constant, `1` for a sine function describing area through
  time.

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

a named list with the updated effective rates.

## See also

[`update_max_rates()`](https://rsetienne.github.io/DAISIE/reference/update_max_rates.md)

Other rate calculations:
[`DAISIE_plot_area()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_plot_area.md),
[`get_ana_rate()`](https://rsetienne.github.io/DAISIE/reference/get_ana_rate.md),
[`get_ext_rate()`](https://rsetienne.github.io/DAISIE/reference/get_ext_rate.md),
[`get_immig_rate()`](https://rsetienne.github.io/DAISIE/reference/get_immig_rate.md),
[`island_area()`](https://rsetienne.github.io/DAISIE/reference/island_area.md),
[`island_area_vector()`](https://rsetienne.github.io/DAISIE/reference/island_area_vector.md),
[`update_max_rates()`](https://rsetienne.github.io/DAISIE/reference/update_max_rates.md)
