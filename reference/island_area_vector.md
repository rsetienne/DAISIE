# Computes island_area, but takes vector as argument (needed by )

Computes island_area, but takes vector as argument (needed by )

## Usage

``` r
island_area_vector(
  timeval,
  area_pars,
  island_ontogeny,
  sea_level,
  total_time,
  peak
)
```

## Arguments

- timeval:

  current time of simulation

- area_pars:

  A vector similar to list produced by create_area_pars

  - \[1\]: maximum area

  - \[2\]: value from 0 to 1 indicating where in the island's history
    the peak area is achieved

  - \[3\]: total island age

- island_ontogeny:

  a string describing the type of island ontogeny. Can be `NULL`,
  `"beta"` for a beta function describing area through time.

- sea_level:

  a numeric describing the type of sea level.

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
[`get_immig_rate()`](https://rsetienne.github.io/DAISIE/reference/get_immig_rate.md),
[`island_area()`](https://rsetienne.github.io/DAISIE/reference/island_area.md),
[`update_max_rates()`](https://rsetienne.github.io/DAISIE/reference/update_max_rates.md),
[`update_rates()`](https://rsetienne.github.io/DAISIE/reference/update_rates.md)

## Author

Pedro Neves
