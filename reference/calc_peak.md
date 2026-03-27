# Calculates the peak of ontogeny curve (beta function)

Calculates the peak of ontogeny curve (beta function)

## Usage

``` r
calc_peak(total_time, area_pars)
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

## Value

numeric
