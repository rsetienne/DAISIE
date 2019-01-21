\name{DAISIE_plot_island}
\alias{DAISIE_plot_island}
\title{ Plot colonisation and branching time of species found on an island dataset.}
\description{ Produces an image with the times of colonisation, branching times, clade name and clade status for all clades found in a given dataset. }
\usage{
DAISIE_plot_island(
island,
island_age=NA)
}

\arguments{
\item{island}{Island data object. Can be in DAISIE list format (see Galapagos_datalist and DAISIE_data_prep for examples) or in table format (see Galapagos_datatable for an example)}
\item{island_age}{If island input is in table format, the age of the island must be specified. If island input is in DAISIE list format, this option will override the island age specified in the island list. }
}
\details{R plot showing for each clade in the island object: time of colonisation, branching times, species status (endemic or non-endemic), total number of species in clade (n), number of species not sampled (m). }
\value{
R plot.
}
\references{
Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium and non-equilibrium dynamics simultaneously operate in the Galapagos islands. Ecology Letters 18: 844-852.}
\author{
Luis Valente
}
\seealso{
 \code{\link{DAISIE_sim}}
\code{\link{DAISIE_plot_age_diversity}}
}
\examples{ 

### Plot Galapagos islands dataset
data(Galapagos_datalist)
DAISIE_plot_island(Galapagos_datalist)

}
\keyword{ models }