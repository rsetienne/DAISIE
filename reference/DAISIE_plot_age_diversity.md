# Plot clade age against clade diversity.

Plots clade age against clade diversity for all clades for which
colonisation time is known.

## Usage

``` r
DAISIE_plot_age_diversity(
  island,
  title = "Clade age vs clade diversity",
  island_age = NA
)
```

## Arguments

- island:

  Island data object. Can be in DAISIE list format (see
  Galapagos_datalist and DAISIE_data_prep for examples) or in table
  format (see Galapagos_datatable for an example).

- title:

  Title of the plot

- island_age:

  Age of island in appropriate units. In `DAISIE_plot_age_diversity()`
  and
  [`DAISIE_plot_island()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_plot_island.md)
  if island input is in table format, the age of the island must be
  specified. If island input is in DAISIE list format, this option will
  override the island age specified in the island list.

## Value

R plot.

## Details

R plot showing for each clade in the island object the time of
colonisation and the diversity of the clade. Only clades for which
colonisation time is known are plotted. Blue - endemic; black -
non-endemic. C = number of independent colonisations on island, N =
number of species on the island

## References

Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium and
non-equilibrium dynamics simultaneously operate in the Galapagos
islands. Ecology Letters 18: 844-852.

## See also

[`DAISIE_sim_cr`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim.md),
[`DAISIE_sim_time_dep`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_time_dep.md),
[`DAISIE_sim_cr_shift`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_cr_shift.md),
[`DAISIE_plot_island`](https://rsetienne.github.io/DAISIE/reference/DAISIE_plot_island.md)

## Author

Luis Valente

## Examples

``` r

### Plot Galapagos age-diversity for Galapagos dataset
utils::data(Galapagos_datalist)
DAISIE_plot_age_diversity(Galapagos_datalist)


```
