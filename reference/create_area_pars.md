# Create named list of area parameters

Create named list of area parameters

## Usage

``` r
create_area_pars(
  max_area,
  current_area,
  proportional_peak_t,
  total_island_age,
  sea_level_amplitude,
  sea_level_frequency,
  island_gradient_angle
)
```

## Arguments

- max_area:

  Numeric defining maximum area.

- current_area:

  A numeric with the current island area at present (i.e., at the end of
  the simulation).

- proportional_peak_t:

  Numeric value from 0 to 1 indicating where in the island's history the
  peak area is achieved.

- total_island_age:

  Numeric defining total island age.

- sea_level_amplitude:

  Numeric defining amplitude of area fluctuation from sea level.

- sea_level_frequency:

  Numeric defining frequency of sine wave of area change from sea level.

- island_gradient_angle:

  Numeric defining the angle in degrees specifying the slope of the
  island.

## Value

list of numerical values containing area and sea level parameters for
island ontogeny simulation

## Author

Richel J.C Bilderbeek, Joshua Lambert, Pedro Neves

## Examples

``` r
area_pars <- create_area_pars(
  max_area = 10,
  current_area = 1,
  proportional_peak_t = 0.5,
  total_island_age = 5,
  sea_level_amplitude = 5,
  sea_level_frequency = 10,
  island_gradient_angle = 0


)
```
