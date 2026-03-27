# Updates state of island given sampled event with two trait states.

Makes the event happen by updating island species matrix and species
IDs. What event happens is determined by the sampling in the algorithm.

## Usage

``` r
DAISIE_sim_update_state_trait_dep(
  timeval,
  total_time,
  possible_event,
  maxspecID,
  mainland_spec,
  island_spec,
  stt_table,
  trait_pars
)
```

## Arguments

- timeval:

  Numeric defining current time of simulation.

- total_time:

  Numeric defining the length of the simulation in time units.

- possible_event:

  Numeric defining what event will happen.

- maxspecID:

  Current species IDs.

- mainland_spec:

  Number of mainland species.

- island_spec:

  Matrix with current state of simulation containing number of species.

- stt_table:

  Matrix with number of species at each time step.

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

The updated state of the system, which is a list with the `island_spec`
matrix, an integer `maxspecID` with the most recent ID of species and
the `stt_table`, a matrix with the current species through time table.

## See also

[DAISIE_sim_core_trait_dep](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_core_trait_dep.md)
