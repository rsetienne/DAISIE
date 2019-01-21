\name{Macaronesia_datalist}
\docType{data}
\alias{Macaronesia_datalist}
\title{Colonization and branching times of terrestrial avifaunal clades from Azores, Canary Islands, Cape Verde and Madeira in list format, accepted by DAISIE_ML and DAISIE_loglik_all}
\description{
  A list containing the colonization and branching times of the terrestrial avifauna in 4 archipelagos: Azores, Canary Islands, Cape Verde and Madeira. It is an R list object with the 4 main elements corresponding to each of the archipelagos (e.g. Macaronesia_datalist[[1]] calls the Azores data). Each of the four elements is then made of several elemants:\cr \cr
  The first element of the list for an archipelago has two components: \cr \cr
  \code{$island_age} - the island age \cr
  \code{$not_present} - the number of mainland lineages that are not present on the island \cr \cr 
  The following elements of the list each contains information on a single colonist lineage on the island and has 5 components:\cr \cr
  \code{$colonist_name} - the name of the species or clade that colonized the island \cr  
  \code{$branching_times} - island age and stem age of the population/species in the case of Non-endemic, Non-endemic_MaxAge, Endemic_MaxAge and Endemic anagenetic species. For cladogenetic species the island age and branching times of the radiation including the stem age of the radiation are shown.\cr                          
  \code{$stac} - the status of the colonist \cr \cr
  * Non_endemic_MaxAge: 1 \cr
  * Endemic: 2 \cr
  * Endemic&Non_Endemic: 3 \cr
  * Non_endemic: 4 \cr \cr
  * Endemic_MaxAge: 5 \cr \cr

  \code{$missing_species} - number of island species that were not sampled for particular clade (only applicable for endemic clades) \cr
  \code{$type1or2} - whether the colonist belongs to type 1 or type 2. In this dataset all are equal to 1. \cr
}
\usage{
data(Macaronesia_datalist)
}
\format{
  A list with 4 main elements for each archipelago. Each element has several sub-elements.
}
\source{

Valente L., Illera J.C, Havenstein K., Pallien T., Etienne R.S., Tiedemann R. Macroevolutionary dynamics in Atlantic island avifaunas: were MacArthur & Wilson right about equilibrium? Under review.
}
\seealso{
  \code{\link{DAISIE_dataprep}},
  \code{\link{DAISIE_ML}}
}
\keyword{datasets}
