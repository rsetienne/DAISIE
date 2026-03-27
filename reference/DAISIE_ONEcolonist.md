# Convert intermediate output to final simulation output

Convert intermediate output to final simulation output

## Usage

``` r
DAISIE_ONEcolonist(time, island_spec, stt_table)
```

## Arguments

- time:

  Numeric defining the length of the simulation in time units. For
  example, if an island is known to be 4 million years old, setting time
  = 4 will simulate the entire life span of the island; setting time = 2
  will stop the simulation at the mid-life of the island.

- island_spec:

  Matrix with current state of simulation containing number of species.

- stt_table:

  Matrix with number of species at each time step.

## Value

a list with these elements:

- \[1\]: `stt_table`, the same stt_table as put in.

- \[2\]: `branching_times`, a sorted numeric vector, as required by the
  ML estimation functions. The first element always refers to the island
  age. Subsequent elements refer to colonisation, speciation and
  recolonisation times. The most recent recolonisation time, if any is
  always omitted to approximate simulation results to the mathematical
  formulation of the likelihood functions used for MLE.

- \[3\]: `stac`, status of colonist. In this function it can be returned
  as either 2, 4 or 3. If `stac` is 2, then there is only one
  independent colonisation present on the island and the extant species
  are endemic. If stac is 4, then only a singleton endemic is present at
  the present. If stac is 3, then recolonisation occurred, and more than
  one colonising lineage.

- \[4\]: `missing_species`, a numeric value with the number of missing
  species, that is, species not sampled in the phylogeny but present on
  the island. As this code only runs for simulation models, here
  `missing_species` is always set to 0.

- \[5\]: `all_colonisations`, on recolonising lineages only. It is
  comprised of `$event_times` and `$species_type`:

  - `$event_times`:

    ordered numeric vectors containing all events for each extant
    recolonising lineage. This includes all colonisation and branching
    times. Each vector pertains to one colonising lineage.

  - `$species_type`:

    a string. Can be `"A"`, `"C"` or `"I"` depending on whether the
    extant clade is of anagenetic, cladogenetic or immigrant origin,
    respectively.
