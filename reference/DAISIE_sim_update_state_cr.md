# Updates state of island given sampled event for a constant rate case.

Makes the event happen by updating island species matrix and species
IDs. What event happens is determined by the sampling in the algorithm.

## Usage

``` r
DAISIE_sim_update_state_cr(
  timeval,
  total_time,
  possible_event,
  maxspecID,
  mainland_spec,
  island_spec,
  stt_table
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

## Value

The updated state of the system, which is a list with the `island_spec`
matrix, an integer `maxspecID` with the most recent ID of species and
the `stt_table`, a matrix with the current species through time table.

## See also

[DAISIE_sim_core_cr](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim_core_cr.md),
DAISIE_sim_update_state_cr
