\name{Galapagos_datatable}
\docType{data}
\alias{Galapagos_datatable}
\title{Colonization and branching times of 8 terrestrial avifaunal clades in table format}
\description{
  A table containing the colonization and branching times of the terrestrial avifauna in the Galapagos.
  Each row on the table represents and independent colonisation event. The table has four columns. \cr \cr
  \code{$Clade_name} - name of independent colonization event \cr
  \code{$Status} - One of the following categories: \cr
  * Non_endemic: for cases where both island and non-island populations of the species have been sampled) \cr
  * Non_endemic_MaxAge: for cases where island population of the species has not been sampled and only the age of the species is available) \cr
  * Endemic: applicable for both cladogenetic or anagenetic species \cr
  * Endemic&Non_Endemic: when endemic clade and mainland ancestor has re-colonized \cr
  \code{$Missing_species} - Number of island species that were not sampled for particular clade (only applicable for endemic clades) \cr
  \code{$Branching_times} - Stem age of the population/species in the case of Non-endemic, Non-endemic_MaxAge and Endemic anagenetic species. For cladogenetic species these should be branching times of the radiation including the stem age of the radiation.\cr
}
\usage{
data(Galapagos_datatable)
}
\format{
  A table with 8 rows and 4 columns.
}
\source{
Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium and non-equilibrium dynamics simultaneously operate in the Galapagos islands. Ecology Letters 18: 844-852.}
\seealso{
  \code{\link{DAISIE_dataprep}},
  \code{\link{DAISIE_ML}}
}
\keyword{datasets}
