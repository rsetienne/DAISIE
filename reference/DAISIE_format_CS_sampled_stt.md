# Formats clade-specific simulation output into standard DAISIE list output

Formats clade-specific simulation output into standard DAISIE list
output

## Usage

``` r
DAISIE_format_CS_sampled_stt(
  island_replicates,
  time,
  M,
  sample_freq,
  verbose = TRUE,
  trait_pars = NULL
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

- trait_pars:

  A named list containing diversification rates considering two trait
  states created by
  [`create_trait_pars`](https://rsetienne.github.io/DAISIE/reference/create_trait_pars.md):

  - \[1\]:A numeric with the per capita transition rate with state 1

  - \[2\]:A numeric with the per capita immigration rate with state 2

  - \[3\]:A numeric with the per capita extinction rate with state 2

  - \[4\]:A numeric with the per capita anagenesis rate with state 2

  - \[5\]:A numeric with the per capita cladogenesis rate with state 2

  - \[6\]:A numeric with the per capita transition rate with state 2

  - \[7\]:A numeric with the number of species with trait state 2 on
    mainland

## Value

List with CS DAISIE simulation output
