# Calculates the area at a point in time from a beta function

Calculates the area at a point in time from a beta function

## Usage

``` r
calc_Abeta(proptime, proptime_max, peak, Amax)
```

## Arguments

- proptime:

  A numeric from 0 to 1. The proportion of time that has elapsed in the
  simulation, in relation to the total island age (NB: not the
  simulation time, but island age).

- proptime_max:

  A numeric from 0 to 1. The same as `proportional_peak_t`. Indicates,
  in proportion to the total island age when the ontogeny peak should
  occur (i.e. 0.5 means a peak halfway in time).

- peak:

  A numeric value specifying the peakiness (or shaprness) of the
  ontogeny curve. Higher values imply peakier ontogeny. This value is
  internally calculated by
  [`calc_peak()`](https://rsetienne.github.io/DAISIE/reference/calc_peak.md)
  given the area at the present and the `area_pars`.

- Amax:

  A numeric value for maximum island area during the simulation.

## Value

Numeric

## Author

Joshua Lambert, Pedro Neves, Shu Xie
