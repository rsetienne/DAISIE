\name{islands_1type_1000reps}
\alias{islands_1type_1000reps}
\title{1000 islands in DAISIE format simulated with the ML parameters of the CR model for the Galapagos data}
\description{
Each simulated dataset is an element of the list, which can be called using e.g.islands_1type_1000reps[[1]]
Each of the island replicates is a list in itself. The first (e.g. islands_1type_1000reps[[x]][[1]]) element of that list has the following components: \cr
  \code{$island_age} - the island age \cr
  \code{$not_present} - the number of mainland lineages that are not present on the island \cr  
\code{$stt_all} - STT table for all species on the island (nI - number of non-endemic species; nA - number of anagenetic species, nC - number of cladogenetic species, present - number of independent colonisations present )\cr

  The subsequent elements of the list each contain information on a single colonist lineage on the island and has 3 components:\cr
  \code{$branching_times} - island age and stem age of the population/species in the case of Non-endemic, Non-endemic_MaxAge and Endemic anagenetic species. For cladogenetic species these should be island age and branching times of the radiation including the stem age of the radiation.\cr                          
  \code{$stac} - the status of the colonist \cr
  * Non_endemic_MaxAge: 1 \cr
  * Endemic: 2 \cr
  * Endemic&Non_Endemic: 3 \cr
  * Non_endemic: 4 \cr
  \code{$missing_species} - number of island species that were not sampled for particular clade (only applicable for endemic clades) \cr
}
\usage{
data(islands_1type_1000reps)
}
\format{
A list with 1000 items.
}
\source{
Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium and non-equilibrium dynamics simultaneously operate in the Galapagos islands. Ecology Letters 18: 844-852.}
\seealso{
  \code{\link{DAISIE_sim}},
  \code{\link{DAISIE_plot_sims}}
}
\keyword{datasets}
