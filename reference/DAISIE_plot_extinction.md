# Plots extinction rate function through time

Plots extinction rate function through time

## Usage

``` r
DAISIE_plot_extinction(
  total_time,
  area_pars,
  peak,
  mu,
  hyper_pars,
  island_ontogeny = "beta",
  sea_level = "const",
  removed_timepoints,
  resolution,
  extcutoff = 1000
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

- mu:

  A numeric with the per capita extinction rate.

- hyper_pars:

  A named list of numeric hyperparameters for the rate calculations as
  returned by
  [`create_hyper_pars()`](https://rsetienne.github.io/DAISIE/reference/create_hyper_pars.md):

  - \[1\]: is d the scaling parameter for exponent for calculating
    cladogenesis rate

  - \[2\]: is x the exponent for calculating extinction rate

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

- removed_timepoints:

  Positive integer with number of first datapoints to be removed from
  rate plots (to prevent Inf)

- resolution:

  numeric indicating resolution of plot. Should be \< 0.

- extcutoff:

  A numeric with the cutoff for the the maximum extinction rate
  preventing it from being too large and slowing down simulation.

## Value

per capita extinction rate through time plot and dataframe with
extinction at corresponding time

## Author

Pedro Neves
