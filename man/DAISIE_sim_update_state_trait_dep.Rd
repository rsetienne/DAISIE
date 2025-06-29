% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DAISIE_sim_update_state_trait_dep.R
\name{DAISIE_sim_update_state_trait_dep}
\alias{DAISIE_sim_update_state_trait_dep}
\title{Updates state of island given sampled event with two trait states.}
\usage{
DAISIE_sim_update_state_trait_dep(
  timeval,
  total_time,
  possible_event,
  maxspecID,
  mainland_spec,
  island_spec,
  stt_table,
  trait_pars
)
}
\arguments{
\item{timeval}{Numeric defining current time of simulation.}

\item{total_time}{Numeric defining the length of the simulation in time
units.}

\item{possible_event}{Numeric defining what event will happen.}

\item{maxspecID}{Current species IDs.}

\item{mainland_spec}{Number of mainland species.}

\item{island_spec}{Matrix with current state of simulation containing number
of species.}

\item{stt_table}{Matrix with number of species at each time step.}

\item{trait_pars}{A named list containing diversification rates considering
  two trait states created by \code{\link{create_trait_pars}}:
\itemize{
  \item{[1]:A numeric with the per capita transition rate with state 1}
  \item{[2]:A numeric with the per capita immigration rate with state 2}
  \item{[3]:A numeric with the per capita extinction rate with state 2}
  \item{[4]:A numeric with the per capita anagenesis rate with state 2}
  \item{[5]:A numeric with the per capita cladogenesis rate with state 2}
  \item{[6]:A numeric with the per capita transition rate with state 2}
  \item{[7]:A numeric with the number of species with trait state 2 on
   mainland}
}}
}
\value{
The updated state of the system, which is a list with the
\code{island_spec} matrix, an integer \code{maxspecID} with the most recent
ID of species and the \code{stt_table}, a matrix with the current species
through time table.
}
\description{
Makes the event happen by updating island species matrix and species IDs.
What event happens is determined by the sampling in the algorithm.
}
\seealso{
\link{DAISIE_sim_core_trait_dep}
}
\keyword{internal}
