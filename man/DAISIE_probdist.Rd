\name{DAISIE_probdist}                               
\alias{DAISIE_probdist}
\title{The joint distribution of endemics and non-endemics under the DAISIE model}
\description{ This function calculates the joint distribution of the number of endemics and non-endemics for a given set of parameter values, a given mainland species pool size and a given set of times}
\usage{
DAISIE_probdist(
   pars1,
   pars2,
   tvec,
   initEI = c(0,0),
   initprobs = NULL
   )
}

\arguments{
  \item{pars1}{ Vector of model parameters: \cr \cr
    \code{pars1[1]} corresponds to lambda^c (cladogenesis rate) \cr
    \code{pars1[2]} corresponds to mu (extinction rate)  \cr
    \code{pars1[3]} corresponds to K (clade-level carrying capacity) \cr
    \code{pars1[4]} corresponds to gamma (immigration rate) \cr
    \code{pars1[5]} corresponds to lambda^a (anagenesis rate) 
  }
  \item{pars2}{ Vector of settings: \cr \cr
    \code{pars2[1]} corresponds to res, the maximum number of endemics or non-endemics for which the ODE system is solved; this must be much larger than the actual number for which the probability needs to be calculated.) \cr
    \code{pars2[2]} corresponds to M, size of the mainland pool, i.e the number of species that can potentially colonize the island. 
  }
  \item{tvec}{ The times at which the probabilities need to be computed.}
  \item{initEI}{ The initial values for the number of endemics and non-endemics; either this or initprobs must be NULL}
  \item{initprobs}{ The initial probability distribution for the number of endemics and non-endemics; either this or initEI must be NULL}

}
\details{
To obtain a matrix of probabilities with endemics in rows and non-endemics in columns for a certain time, one can run DAISIE_convertprobdist
}
\value{
  A matrix of dimensions 1 + length(tvec) and pars[2]^2 + 1] where the first column contains the times at which the probabilities are evaluated and the other columns contain the joint probabilities. 
}
\references{
Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium and non-equilibrium dynamics simultaneously operate in the Galapagos islands. Ecology Letters 18: 844-852.}
\author{ Rampal S. Etienne }
\examples{
### Compute the probability distribution at t = 4 and t = 8, for a mainland pool
# size of 250 potential colonists and a vector of 5 parameters (cladogenesis, 
# extinction, clade-level carrying capacity, immigration, anagenesis) starting
# from an empty island

DAISIE_probdist(
   pars1 = c(0.3,0.35,Inf,0.75,0.012),
   pars2 = c(100,250),
   tvec = c(4,8),
   initEI = c(0,0),
   initprobs = NULL
   )
}

\keyword{models}