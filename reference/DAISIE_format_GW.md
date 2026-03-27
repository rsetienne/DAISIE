# Formats guild-wide simulation output into standard DAISIE list output

Formats guild-wide simulation output into standard DAISIE list output

## Usage

``` r
DAISIE_format_GW(
  island_replicates,
  time,
  M,
  sample_freq,
  num_guilds,
  verbose = TRUE
)
```

## Arguments

- island_replicates:

  List output from
  [`DAISIE_sim_core_cr()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_core_cr.md),
  [`DAISIE_sim_core_time_dep()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_core_time_dep.md),
  [`DAISIE_sim_core_cr_shift()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_core_cr_shift.md)
  or
  [`DAISIE_sim_min_type2()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_min_type2.md)
  functions. Minimally, this must be a list that has as many elements as
  replicates. Each element must be a list with the elements
  `island_age`, `not_present` and `stt_all`. `stt_all` must be a data
  frame with the column names `Time`, `nI`, `nA`, `nC` and `present`.

- time:

  Numeric defining the length of the simulation in time units. For
  example, if an island is known to be 4 million years old, setting time
  = 4 will simulate the entire life span of the island; setting time = 2
  will stop the simulation at the mid-life of the island.

- M:

  Numeric defining the size of mainland pool, i.e. the number of species
  that can potentially colonize the island.

- sample_freq:

  Numeric specifing the number of units times should be divided by for
  plotting purposes. Larger values will lead to plots with higher
  resolution, but will also run slower.

- num_guilds:

  The number of guilds on the mainland. The number of mainland species
  is divided by the number of guilds when `divdepmodel = "GW"`

- verbose:

  A numeric vector of length 1, which in simulations and
  \`DAISIEdataprep()\` can be \`1\` or \`0\`, where \`1\` gives
  intermediate output should be printed. For ML functions a numeric
  determining if intermediate output should be printed. The default:
  \`0\` does not print, \`1\` prints the initial likelihood and the
  settings that were selected (which parameters are to be optimised,
  fixed or shifted), \`2\` prints the same as \`1 and also the
  intermediate output of the parameters and loglikelihood, while \`3\`
  the same as \`2\` and prints intermediate progress during likelihood
  computation.

## Value

List with GW DAISIE simulation output
