# Count number of species in DAISIE datalist or simulated data.

Calculates various island diversity metrics from island datasets.

## Usage

``` r
DAISIE_count_species(islands, sort_clade_sizes = TRUE)
```

## Arguments

- islands:

  Island datalist or simulated data in DAISIE datalist format. Can be a
  single island (empirical data) generated with DAISIE_dataprep or
  DAISIEprep. Can also be simulated data generated with DAISIE_sim
  function.

- sort_clade_sizes:

  Default sort_clade_sizes = T outputs clade sizes sorted in ascending
  order of number of species. sort_clade_sizes=F outputs clade sizes in
  the same order as they appear in the input datalist.

## Value

The output is a list containing the following items:

- clade_sizes_sorted:

  List showing the total number of species in each island clade
  (including missing species). Each item \[\[i\]\] on the list gives the
  sizes of all clades for a single island. If option sort_clade_sizes =
  T, the clade sizes for are sorted by increasing number of species. If
  option sort_clade_sizes = F the clade sizes are given in the same
  order as in the input datalist.

- size_largest_clade:

  The total number of species in the largest island clade for each
  island.

- mean_clade_size:

  Mean clade size (average of all island clades)

- number_colonisations:

  The total number of colonisations (clades) on each island.

- total_number_species:

  The total number of species on each island. These are the extant
  species at present, including missing species; in case of simulations,
  this is the number of species present on the island at the end of the
  simulation.

## See also

[`DAISIE_dataprep`](https://rsetienne.github.io/DAISIE/reference/DAISIE_dataprep.md),
[`DAISIE_plot_island`](https://rsetienne.github.io/DAISIE/reference/DAISIE_plot_island.md)

## Author

Luis Valente

## Examples

``` r
# Run function with clade sizes in the order they appear in the input data
data("NewZealand_birds_datalist")
species_count <- DAISIE_count_species(NewZealand_birds_datalist)

# Run function with clade sizes in ascending order
species_count_sorted <- DAISIE_count_species(
  NewZealand_birds_datalist,
  sort_clade_sizes = TRUE
)
```
