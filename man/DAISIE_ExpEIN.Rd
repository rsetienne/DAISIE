\name{DAISIE_ExpEIN}                               
\alias{DAISIE_ExpEIN}
\title{The expected number of endemics and non-endemics under the DAISIE model}
\description{ This function calculates the expected number of endemics, non-endemics and the sum of these for a given set of parameter values, a given mainland species pool size and a given time}
\usage{
DAISIE_ExpEIN(
   t,
   pars,
   M,
   initEI = c(0,0)
   )
}
\arguments{
  \item{t}{ The time at which the expectations need to be computed }
  \item{pars}{ Vector of parameters: \cr \cr
    \code{pars[1]} corresponds to lambda^c (cladogenesis rate) \cr
    \code{pars[2]} corresponds to mu (extinction rate)  \cr
    \code{pars[3]} corresponds to K (clade-level carrying capacity) \cr
    \code{pars[4]} corresponds to gamma (immigration rate) \cr
    \code{pars[5]} corresponds to lambda^a (anagenesis rate) 
  }  
  \item{M}{ The size of the mainland pool, i.e the number of species that can potentially colonize the island }
  \item{initEI}{ The initial values of the number of endemics and non-endemics }

}

\value{
 \item{out}{The output is a list with three elements: \cr \cr
  \code{ExpE} The number of endemic species \cr
  \code{ExpI} The number of non-endemic species \cr
  \code{ExpN} The sum of the number of endemics and non-endemics
}
}
\references{
Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium and non-equilibrium dynamics simultaneously operate in the Galapagos islands. Ecology Letters 18: 844-852.}
\author{ Rampal S. Etienne }
\examples{
### Compute the expected values at t = 4, for a mainland pool size of 1000 potential
# colonists and a vector of 5 parameters (cladogenesis, extinction, clade-level carrying
# capacity, immigration, anagenesis)
 
DAISIE_ExpEIN(
   t = 4,
   pars = c(0.5,0.1,Inf,0.01,0.4),
   M = 1000
   )
}
\keyword{models}