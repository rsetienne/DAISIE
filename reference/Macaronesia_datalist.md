# Colonization and branching times of terrestrial avifaunal clades from Azores, Canary Islands, Cape Verde and Madeira in list format, accepted by DAISIE_ML and DAISIE_loglik_all

A list containing the colonization and branching times of the
terrestrial avifauna in 4 archipelagos: Azores, Canary Islands, Cape
Verde and Madeira. It is an R list object with the 4 main elements
corresponding to each of the archipelagos (e.g.
Macaronesia_datalist\[\[1\]\] calls the Azores data). Each of the four
elements is then made of several elemants:  
  
The first element of the list for an archipelago has two components:  
  
`$island_age` - the island age  
`$not_present` - the number of mainland lineages that are not present on
the island  
  
The following elements of the list each contains information on a single
colonist lineage on the island and has 5 components:  
  
`$colonist_name` - the name of the species or clade that colonized the
island  
`$branching_times` - island age followed by stem age of the
population/species in the case of Non-endemic, Non-endemic_MaxAge
species and Endemic species with no close relatives on the island. For
endemic clades with more than one species on the island (cladogenetic
clades/ radiations) these should be island age followed by the branching
times of the island clade including the stem age of the clade.  
  
\* Non_endemic_MaxAge: 1  
\* Endemic: 2  
\* Endemic&Non_Endemic: 3  
\* Non_endemic: 4  
  
\* Endemic_MaxAge: 5  
  
`$missing_species` - number of island species that were not sampled for
particular clade (only applicable for endemic clades)  
`$type1or2` - whether the colonist belongs to type 1 or type 2. In this
dataset all are equal to 1.  

## Format

A list with 4 main elements for each archipelago. Each element has
several sub-elements.

## Source

Valente L., Illera J.C, Havenstein K., Pallien T., Etienne R.S.,
Tiedemann R. Equilibrium bird species diversity in Atlantic islands.
2017 Current Biology, 27, 1660-1666.

## See also

[`DAISIE_dataprep`](https://rsetienne.github.io/DAISIE/reference/DAISIE_dataprep.md),
[`DAISIE_ML`](https://rsetienne.github.io/DAISIE/reference/DAISIE_ML.md)
