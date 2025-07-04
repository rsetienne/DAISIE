% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DAISIE_ExpEIN.R
\name{DAISIE_ExpEIN2}
\alias{DAISIE_ExpEIN2}
\title{The expected number of endemics and non-endemics under the DAISIE model}
\usage{
DAISIE_ExpEIN2(
  tvec,
  pars,
  M,
  initEI = NULL,
  res = 1000,
  ddmodel = 11,
  methode = "ode45",
  reltolint = 1e-16,
  abstolint = 1e-16
)
}
\arguments{
\item{tvec}{The times at which the probabilities need to be computed.}

\item{pars}{A numeric vector containing the model parameters:
\itemize{
  \item{\code{pars[1]}: lambda^c (cladogenesis rate)}
  \item{\code{pars[2]}: mu (extinction rate)}
  \item{\code{pars[3]}: K (carrying capacity), set K=Inf for diversity
  independence.}
  \item{\code{pars[4]}: gamma (immigration rate)}
  \item{\code{pars[5]}: lambda^a (anagenesis rate)}
  \item{\code{pars[6]}: lambda^c (cladogenesis rate) for either type 2 species
  or rate set 2 in rate shift model}
  \item{\code{pars[7]}: mu (extinction rate) for either type 2 species or rate
  set 2 in rate shift model}
  \item{\code{pars[8]}: K (carrying capacity) for either type 2 species or rate
  set 2 in rate shift model, set K=Inf for diversity independence.}
  \item{\code{pars[9]}: gamma (immigration rate) for either type 2 species
  or rate set 2 in rate shift model}
  \item{\code{pars[10]}: lambda^a (anagenesis rate) for either type 2
  species or rate set 2 in rate shift model}
}
Elements 6:10 are required only when type 2 species are included
or in the rate shift model. For \code{\link{DAISIE_sim_relaxed_rate}()}
\code{pars[6]} is the standard deviation of the gamma distribution for the
relaxed parameter and the parameter chosen by the \code{relaxed_par}
argument is the mean of the gamma distribution for the relaxed parameter.}

\item{M}{Numeric defining the size of mainland pool, i.e. the number of
species that can potentially colonize the island.}

\item{initEI}{The initial values for the number of endemics and
  non-endemics. In \code{\link{DAISIE_probdist}()} or
\code{\link{DAISIE_margprobdist}()} either this or initprobs must be NULL. In
\code{\link{DAISIE_numcol}()} when it is NULL, it is assumed that the island
  is empty.}

\item{res}{Sets the maximum number of species for which a probability must
be computed, must be larger than the size of the largest clade.}

\item{ddmodel}{Sets the model of diversity-dependence: \cr \cr
\itemize{
  \item{ddmodel = 0 : no diversity dependence}
  \item{ddmodel = 1 : linear dependence in speciation rate}
  \item{ddmodel = 11: linear dependence in speciation rate and in immigration rate}
  \item{ddmodel = 2 : exponential dependence in speciation rate}
  \item{ddmodel = 21: exponential dependence in speciation rate and in immigration
 rate}}}

\item{methode}{Method of the ODE-solver. Supported Boost \code{ODEINT}
solvers (steppers) are:
\code{"odeint::runge_kutta_cash_karp54"}
\code{"odeint::runge_kutta_fehlberg78"}
\code{"odeint::runge_kutta_dopri5"}
\code{"odeint::bulirsch_stoer"}
without \code{odeint::}-prefix, \code{\link[deSolve]{ode}} method is
assumed. The default method overall is
\code{"lsodes"} for \code{\link{DAISIE_ML_CS}()}
and \code{"ode45"} from \code{\link[deSolve]{ode}()} for
\code{\link{DAISIE_ML_IW}()}.}

\item{reltolint}{Numeric relative tolerance of the integration}

\item{abstolint}{Numeric absolute tolerance of the integration}
}
\value{
\item{tot_expEIN}{The output is a list with three elements: \cr \cr
\code{ExpE} The number of endemic species at the times in tvec\cr
\code{ExpI} The number of non-endemic species at the times in tvec\cr
\code{ExpN} The sum of the number of endemics and non-endemics at the times
in tvec}
}
\description{
This function calculates the expected number of endemics, non-endemics and
the sum of these for a given set of parameter values, a given mainland
species pool size and a given time, where there can be diversity-dependence
}
\examples{
DAISIE_ExpEIN2(tvec = c(0.000001,0.5,0.75,1),
                         pars = c(0.3,0.1,10,1,0.1),
                         M = 1000,
                         initEI = rbind(c(1,0),c(2,0),c(0,1)))
}
\author{
Rampal S. Etienne
}
