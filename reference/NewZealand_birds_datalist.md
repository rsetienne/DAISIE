# Colonization and branching times of New Zealand birds.

A list containing the colonization and branching times of the birds of
New Zealand. Main dataset used in Valente, Etienne, Garcia-R (2019)
Current Biology. Island age 52 Myr and mainland pool size of 1000
species.  
The first element of the list has two components:  
  
`$island_age` - the island age  
`$not_present` - the number of mainland lineages that are not present on
the island  
  
The following elements of the list each contain information on a single
colonist lineage on the island and has 5 components:  
  
`$colonist_name` - the name of the species or clade that colonized the
island  
`$branching_times` - island age followed by stem age of the
population/species in the case of Non-endemic, Non-endemic_MaxAge
species and Endemic species with no close relatives on the island. For
endemic clades with more than one species on the island (cladogenetic
clades/ radiations) these should be island age followed by the branching
times of the island clade including the stem age of the clade.  
`$stac` - the status of the colonist  
  
\* Non_endemic_MaxAge: 1  
\* Endemic: 2  
\* Endemic&Non_Endemic: 3  
\* Non_endemic: 4  
\* Endemic_MaxAge: 5 or 6  
  
`$missing_species` - number of island species that were not sampled for
particular clade (only applicable for endemic clades)  
`$type1or2` - whether the colonist belongs to type 1 or type 2. In this
dataset all are equal to 1.  

## Format

A list with 40 elements, the first of which contains 2 elements and the
following 39 containing 5 components.

## Source

Valente L, Etienne RS, Garcia-R JC (2019) Deep Macroevolutionary Impact
of Humans on New Zealand’s Unique Avifauna. Current Biology, 29,
2563–2569.  

## See also

[`DAISIE_dataprep`](https://rsetienne.github.io/DAISIE/reference/DAISIE_dataprep.md),
[`DAISIE_ML`](https://rsetienne.github.io/DAISIE/reference/DAISIE_ML.md),
[`DAISIE_SR_ML`](https://rsetienne.github.io/DAISIE/reference/DAISIE_SR_ML.md)
