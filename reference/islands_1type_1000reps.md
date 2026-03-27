# 1000 islands in DAISIE format simulated with the ML parameters of the CR model for the Galapagos data

Each simulated dataset is an element of the list, which can be called
using e.g. islands_1type_1000reps\[\[1\]\] Each of the island replicates
is a list in itself. The first (e.g.
islands_1type_1000reps\[\[x\]\]\[\[1\]\]) element of that list has the
following components:  
`$island_age` - the island age  
`$not_present` - the number of mainland lineages that are not present on
the island  
`$stt_all` - STT table for all species on the island (nI - number of
non-endemic species; nA - number of anagenetic species, nC - number of
cladogenetic species, present - number of independent colonisations
present )  
The subsequent elements of the list each contain information on a single
colonist lineage on the island and has 3 components:  
`$branching_times` - island age and stem age of the population/species
in the case of Non-endemic, Non-endemic_MaxAge and Endemic anagenetic
species. For cladogenetic species these should be island age and
branching times of the radiation including the stem age of the
radiation.  
`$stac` - the status of the colonist  
\* Non_endemic_MaxAge: 1  
\* Endemic: 2  
\* Endemic&Non_Endemic: 3  
\* Non_endemic: 4  
`$missing_species` - number of island species that were not sampled for
particular clade (only applicable for endemic clades)  

## Format

A list with 1000 items.

## Source

Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium and
non-equilibrium dynamics simultaneously operate in the Galapagos
islands. Ecology Letters 18: 844-852.

## See also

[`DAISIE_sim_cr()`](https://rsetienne.github.io/DAISIE/reference/DAISIE_sim.md),
[`DAISIE_plot_sims`](https://rsetienne.github.io/DAISIE/reference/DAISIE_plot_sims.md)
