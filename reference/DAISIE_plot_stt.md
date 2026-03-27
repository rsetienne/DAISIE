# Create the Species-Through-Time plot. This is used to visualize the output of DAISIE_sim functions

Create the Species-Through-Time plot. This is used to visualize the
output of DAISIE_sim functions

## Usage

``` r
DAISIE_plot_stt(
  plot_plus_one = TRUE,
  time,
  plot_lists = plot_lists,
  type = type
)
```

## Arguments

- plot_plus_one:

  Boolean to indicate to plot all values plus one. Set to `TRUE` for
  default behavior. Set to `FALSE` to plot all values without adding
  one. Only works when there is one type of species.

- time:

  Numeric defining the length of the simulation in time units. For
  example, if an island is known to be 4 million years old, setting time
  = 4 will simulate the entire life span of the island; setting time = 2
  will stop the simulation at the mid-life of the island.

- plot_lists:

  List of lists containing average and quantile species through time.

- type:

  String to indicate if stt of all species or all possible stt should be
  plotted. Default is `"all_species"`, `"type1_species"` or
  `"type2_species"` should be plotted.
