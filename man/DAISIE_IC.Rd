% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DAISIE_IC.R
\name{DAISIE_IC}
\alias{DAISIE_IC}
\title{Calculates information criterion from DAISIE ML estimates?}
\usage{
DAISIE_IC(
  datalist,
  initparsopt,
  idparsopt,
  parsfix,
  idparsfix,
  endmc = 1000,
  res = 100,
  cond = 0,
  ddmodel = 0
)
}
\arguments{
\item{datalist}{Data object containing information on colonisation and
  branching times. This object can be generated using the DAISIE_dataprep
  function, which converts a user-specified data table into a data object,
  but the object can of course also be entered directly.
  It is an R list object with the following elements.\cr
  The first element of the list has two or three components: \cr \cr
  \code{$island_age} - the island age \cr
  Then, depending on whether a distinction between types is
  made, we have:\cr
 \code{$not_present} - the number of mainland lineages
  that are not present on the island \cr
  or:\cr
 \code{$not_present_type1} -
  the number of mainland lineages of type 1 that are not present on the
  island \cr
 \code{$not_present_type2} - the number of mainland lineages of
  type 2 that are not present on the island \cr \cr
 The remaining elements of
  the list each contains information on a single colonist lineage on the
  island and has 5 components:\cr \cr
 \code{$colonist_name} - the name of the
  species or clade that colonized the island \cr
 \code{$branching_times} -
  island age followed by stem age of the population/species in the case of
  Non-endemic, Non-endemic_MaxAge species and Endemic species with no close relatives
  on the island. For endemic clades with more than one species on the island
  (cladogenetic clades/ radiations) these should be island age followed by the
  branching times of the island clade including the stem age of the clade\cr
 \code{$stac} - the status of the colonist \cr \cr
- Non_endemic_MaxAge: 1 \cr
- Endemic: 2 \cr
- Endemic&Non_Endemic: 3 \cr
- Non_Endemic: 4 \cr
- Endemic_Singleton_MaxAge: 5 \cr
- Endemic_Clade_MaxAge: 6 \cr
- Endemic&Non_Endemic_Clade_MaxAge: 7 \cr
- Non_endemic_MaxAge_MinAge: 8 \cr
- Endemic_Singleton_MaxAge_MinAge: 9 \cr \cr
 \code{$missing_species} - number of island species that were not
  sampled for particular clade (only applicable for endemic clades) \cr
 \code{$type1or2} - whether the colonist belongs to type 1 or type 2 \cr}

\item{initparsopt}{The initial values of the parameters that must be
optimized, they are all positive.}

\item{idparsopt}{The ids of the parameters that must be optimized. The ids
are defined as follows: \cr \cr id = 1 corresponds to lambda^c
(cladogenesis rate) \cr id = 2 corresponds to mu (extinction rate) \cr
id = 3 corresponds to K (clade-level carrying capacity) \cr id = 4
corresponds to gamma (immigration rate) \cr id = 5 corresponds to lambda^a
(anagenesis rate) \cr id = 6 corresponds to lambda^c (cladogenesis rate)
for an optional subset of the species \cr id = 7 corresponds to mu
(extinction rate) for an optional subset of the species\cr id = 8
corresponds to K (clade-level carrying capacity) for an optional subset of
the species\cr id = 9 corresponds to gamma (immigration rate) for an
optional subset of the species\cr id = 10 corresponds to lambda^a
(anagenesis rate) for an optional subset of the species\cr id = 11
corresponds to p_f (fraction of mainland species that belongs to the second
subset of species.}

\item{parsfix}{The values of the parameters that should not be optimized.}

\item{idparsfix}{The ids of the parameters that should not be optimized,
e.g. c(1,3) if lambda^c and K should not be optimized.}

\item{endmc}{Numeric for how many simulations should run.}

\item{res}{Sets the maximum number of species for which a probability must
be computed, must be larger than the size of the largest clade.}

\item{cond}{cond = 0 : conditioning on island age \cr cond = 1 :
conditioning on island age and non-extinction of the island biota \cr.
cond > 1 : conditioning on island age and having at least cond colonizations
on the island. This last option is not yet available for the IW model \cr}

\item{ddmodel}{Sets the model of diversity-dependence: \cr \cr
\itemize{
  \item{ddmodel = 0 : no diversity dependence}
  \item{ddmodel = 1 : linear dependence in speciation rate}
  \item{ddmodel = 11: linear dependence in speciation rate and in immigration rate}
  \item{ddmodel = 2 : exponential dependence in speciation rate}
  \item{ddmodel = 21: exponential dependence in speciation rate and in immigration
 rate}}}
}
\value{
List of two numerics with WIC and AICb
}
\description{
Calculates information criterion from DAISIE ML estimates?
}
