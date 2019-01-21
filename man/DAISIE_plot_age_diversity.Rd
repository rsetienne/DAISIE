\name{DAISIE_plot_age_diversity}
\alias{DAISIE_plot_age_diversity}
\title{ Plot clade age against clade diversity.}
\description{ Plots clade age against clade diversity for all clades for which colonisation time is known. }
\usage{
DAISIE_plot_age_diversity(
  island,
  title = 'Clade age vs clade diversity',
  island_age = NA) 
}
\arguments{
\item{island}{Island data object. Can be in DAISIE list format (see Galapagos_datalist and DAISIE_data_prep for examples) or in table format (see Galapagos_datatable for an example)}
\item{title}{Title of the plot}
\item{island_age}{If island input is in table format, the age of the island must be specified. If island input is in DAISIE list format, this option will override the island age specified in the island list. }
}
\details{R plot showing for each clade in the island object the time of colonisation and the diversity of the clade. Only clades for which colonisation time is known are plotted. Blue - endemic; black - non-endemic. C = number of independent colonisations on island, N = number of species on the island }
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
\code{\link{DAISIE_plot_island}}
}
\examples{ 

### Plot Galapagos age-diversity for Galapagos dataset
data(Galapagos_datalist)
DAISIE_plot_age_diversity(Galapagos_datalist)

}
\keyword{ models }