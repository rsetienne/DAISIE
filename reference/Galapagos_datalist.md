# Colonization and branching times of 8 terrestrial avifaunal clades in list format, accepted by DAISIE_ML and DAISIE_loglik_all

A list containing the colonization and branching times of the
terrestrial avifauna in the Galapagos where no distinction is made
between types of colonists. This list can be generated using the
DAISIE_dataprep function, which converts a user-specified data table
into a data object, but the object can of course also be entered
directly. It is an R list object with the following elements.  
  
The first element of the list has two components:  
  
`$island_age` - the island age  
`$not_present` - the number of mainland lineages that are not present on
the island  
  
The following 8 elements of the list each contains information on a
single colonist lineage on the island and has 5 components:  
  
`$colonist_name` - the name of the species or clade that colonized the
island  
`$branching_times` - island age followed by stem age of the
population/species in the case of Non-endemic, Non-endemic_MaxAge
species and Endemic species with no close relatives on the island. For
endemic clades with more than one species on the island (cladogenetic
clades/ radiations) these should be island age followed by the branching
times of the island cladeincluding the stem age of the clade.  
`$stac` - the status of the colonist  
  
\* Non_endemic_MaxAge: 1  
\* Endemic: 2  
\* Endemic&Non_Endemic: 3  
\* Non_endemic: 4  
  
`$missing_species` - number of island species that were not sampled for
particular clade (only applicable for endemic clades)  
`$type1or2` - whether the colonist belongs to type 1 or type 2. In this
dataset all are equal to 1.  

## Format

A list with 9 elements the first of which contains 2 elements and the
following 8 containing 5 components.

## Source

Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium and
non-equilibrium dynamics simultaneously operate in the Galapagos
islands. Ecology Letters 18: 844-852.

## See also

[`DAISIE_dataprep`](https://rsetienne.github.io/DAISIE/reference/DAISIE_dataprep.md),
[`DAISIE_ML`](https://rsetienne.github.io/DAISIE/reference/DAISIE_ML.md)
