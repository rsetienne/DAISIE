\name{Galapagos_datalist_2types}
\docType{data}
\alias{Galapagos_datalist_2types}
\title{Colonization and branching times of 8 terrestrial avifaunal clades in list format, accepted by DAISIE_ML and DAISIE_loglik_all}
\description{
  A list containing the colonization and branching times of the terrestrial avifauna in the Galapagos. This list can be generated using the DAISIE_dataprep function, which converts a user-specified data table into a data object, but the object can of course also be entered directly. It is an R list object with the following elements.\cr \cr
  The first element of the list has three components: \cr \cr
  \code{$island_age} - the island age \cr
  \code{$not_present_type1} - the number of mainland lineages of type 1 that are not present on the island \cr  
  \code{$not_present_type2} - the number of mainland lineages of type 2 that are not present on the island \cr \cr 
  The following 8 elements of the list each contains information on a single colonist lineage on the island and has 5 components:\cr \cr
  \code{$colonist_name} - the name of the species or clade that colonized the island \cr  
  \code{$branching_times} - island age and stem age of the population/species in the case of Non-endemic, Non-endemic_MaxAge and Endemic anagenetic species. For cladogenetic species these should be island age and branching times of the radiation including the stem age of the radiation.\cr                          
  \code{$stac} - the status of the colonist \cr \cr
  * Non_endemic_MaxAge: 1 \cr
  * Endemic: 2 \cr
  * Endemic&Non_Endemic: 3 \cr
  * Non_endemic: 4 \cr \cr
  \code{$missing_species} - number of island species that were not sampled for particular clade (only applicable for endemic clades) \cr
  \code{$type1or2} - whether the colonist belongs to type 1 or type 2. In this dataset only the finches are type 2 \cr
}
\usage{data(Galapagos_datalist_2types)}
\format{
  A list with 9 elements the first of which contains 3 elements and the following 8 containing 5 components.
}
\source{
Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium and non-equilibrium dynamics simultaneously operate in the Galapagos islands. Ecology Letters 18: 844-852.}
\seealso{
  \code{\link{DAISIE_dataprep}},
  \code{\link{DAISIE_ML}}
}
\keyword{datasets}
