# Converts simulation output into island output

Converts simulation output into island output

## Usage

``` r
DAISIE_create_island(
  stt_table,
  total_time,
  island_spec,
  mainland_n,
  trait_pars = NULL
)
```

## Arguments

- stt_table:

  Matrix with number of species at each time step.

- total_time:

  Numeric defining the length of the simulation in time units.

- island_spec:

  Matrix with current state of simulation containing number of species.

- mainland_n:

  A numeric stating the number of mainland species, that is the number
  of species that can potentially colonize the island. If using a
  clade-specific diversity dependence, this value is set to 1. If using
  an island-wide diversity dependence, this value is set to the number
  of mainland species.

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

list with the island information, composed stt table, branching times of
extant species, status of species on the island and number of missing
species.
