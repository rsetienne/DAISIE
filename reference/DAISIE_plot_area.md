# Plots island area function through time

Plots island area function through time

## Usage

``` r
DAISIE_plot_area(
  total_time,
  area_pars,
  peak,
  island_ontogeny = "beta",
  resolution,
  sea_level = "const"
)
```

## Arguments

- total_time:

  Numeric defining the length of the simulation in time units.

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

- resolution:

  numeric indicating resolution of plot. Should be \< 0.

- sea_level:

  In
  [`DAISIE_sim_time_dep()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_time_dep.md)
  and plotting a string describing the type of sea level. Can be
  `"const"` or `"sine"` for a sine function describing area through
  time.  
  In all other functions a numeric describing the type of sea level. Can
  be `0` for constant, `1` for a sine function describing area through
  time.

## Value

a plot with the area size through time

## See also

Other rate calculations:
[`get_ana_rate()`](https://rsetienne.github.io/DAISIE/reference/get_ana_rate.md),
[`get_ext_rate()`](https://rsetienne.github.io/DAISIE/reference/get_ext_rate.md),
[`get_immig_rate()`](https://rsetienne.github.io/DAISIE/reference/get_immig_rate.md),
[`island_area()`](https://rsetienne.github.io/DAISIE/reference/island_area.md),
[`island_area_vector()`](https://rsetienne.github.io/DAISIE/reference/island_area_vector.md),
[`update_max_rates()`](https://rsetienne.github.io/DAISIE/reference/update_max_rates.md),
[`update_rates()`](https://rsetienne.github.io/DAISIE/reference/update_rates.md)
