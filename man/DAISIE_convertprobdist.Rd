\name{DAISIE_convertprobdist}                               
\alias{DAISIE_convertprobdist}
\title{Converts the joint distribution of endemics and non-endemics under the DAISIE model to list format}
\description{ This function converts the joint distribution of the number of endemics and non-endemics from the matrix format of DAISIE_probdist to a list format}
\usage{
DAISIE_convertprobdist(
   pb
   )
}

\arguments{
  \item{pb}{ Probability distribution in matrix format as output by DAISIE_probdist}
}
\value{
 A list of length nrow(pb) containing matrices of square dimensions of size sqrt(ncol - 1) containing the joint probabilities with endemics in the rows and non-endemics in the columns. The last element of the list is a vector a times at which the joint probability distribution is evaluated. 
}
\references{
Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium and non-equilibrium dynamics simultaneously operate in the Galapagos islands. Ecology Letters 18: 844-852.}
\author{ Rampal S. Etienne }
\examples{
### Compute the probability distribution at t = 4 and t = 8, for a mainland pool
# size of 250 potential colonists and a vector of 5 parameters (cladogenesis, extinction,
# clade-level carrying capacity, immigration, anagenesis) starting from an empty
# island; output in list format

pb = DAISIE_probdist(
   pars1 = c(0.3,0.35,Inf,0.75,0.012),
   pars2 = c(100,250),
   tvec = c(4,8),
   initEI = c(0,0),
   initprobs = NULL
   )
DAISIE_convertprobdist(pb)
}

\keyword{models}