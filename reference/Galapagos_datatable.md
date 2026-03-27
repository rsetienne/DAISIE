# Colonization and branching times of 8 terrestrial avifaunal Galápagos clades in table format.

A table containing the colonization and branching times of the
terrestrial avifauna in the Galápagos. Each row on the table represents
and independent colonisation event. The table has four columns.  
  
`$Clade_name` - name of independent colonization event  
`$Status` - One of the following categories:  
\* Non_endemic: for non-endemic island species when an approximate time
of colonisation is known  
\* Non_endemic_MaxAge: for non-endemic island species when colonisation
time is unknown  
\* Endemic: for endemic species when an approximate colonisation time is
known  
\* "Endemic_MaxAge": applies to endemic species or endemic clades for
cases where the colonisation time is unknown, or when the user wants to
specify an upper bound for colonisation. This could for example apply to
endemic species that have recently gone extinct because of anthropogenic
causes, and which are not included in the phylogeny ("NA" should be
given in the branching times column). It could also apply to insular
radiations with long stem branches, for which the time of the first
cladogenetic event is known, but the precise time of colonisation is
not.  
\* Endemic&Non_Endemic: when endemic clade and mainland ancestor has
re-colonized  

`$Missing_species` - Number of island species that were not sampled for
particular clade (only applicable for endemic clades)  
`$Branching_times` - Stem age of the population/species in the case of
"Non_endemic", "Non_endemic_MaxAge" and "Endemic" species with no extant
close relatives on the island. Set "NA" if colonisation time unknown and
no upper bound is known. For "Endemic" cladogenetic species these should
be branching times of the radiation, including the stem age of the
radiation (colonisation time estimate).  

## Format

A table with 8 rows and 4 columns.

## Source

Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium and
non-equilibrium dynamics simultaneously operate in the Galapagos
islands. Ecology Letters 18: 844-852.

## See also

[`DAISIE_dataprep`](https://rsetienne.github.io/DAISIE/reference/DAISIE_dataprep.md),
[`DAISIE_ML`](https://rsetienne.github.io/DAISIE/reference/DAISIE_ML.md)
