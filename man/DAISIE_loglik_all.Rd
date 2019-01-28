\name{DAISIE_loglik_all}
\alias{DAISIE_loglik_all}
\alias{DAISIE_loglik_CS}
\title{
Computes the loglikelihood of the DAISIE model given data and a set of model parameters
}
\description{
Computes the loglikelihood of the DAISIE model given colonization and branching times for lineages on an island, and a set of model parameters
}
\usage{
DAISIE_loglik_all(
   pars1,
   pars2,
   datalist,
   methode = "lsodes"
   )
}
\arguments{
  \item{pars1}{ Contains the model parameters: \cr \cr
   \code{pars1[1]} corresponds to lambda^c (cladogenesis rate) \cr
   \code{pars1[2]} corresponds to mu (extinction rate) \cr
   \code{pars1[3]} corresponds to K (clade-level carrying capacity) \cr
   \code{pars1[4]} corresponds to gamma (immigration rate) \cr
   \code{pars1[5]} corresponds to lambda^a (anagenesis rate) \cr   
   \code{pars1[6]} corresponds to lambda^c (cladogenesis rate) for an optional subset of the species \cr
   \code{pars1[7]} corresponds to mu (extinction rate) for an optional subset of the species\cr
   \code{pars1[8]} corresponds to K (clade-level carrying capacity) for an optional subset of the species\cr
   \code{pars1[9]} corresponds to gamma (immigration rate) for an optional subset of the species\cr
   \code{pars1[10]} corresponds to lambda^a (anagenesis rate) for an optional subset of the species\cr   
   \code{pars1[11]} corresponds to p_f (fraction of mainland species that belongs to the second subset of species\cr
   The elements 6:10 and 11 are optional, that is,pars1 should either contain 5, 10 or 11 elements. If 10, then the fraction of potential colonists of type 2 is computed from the data. If 11, then pars1[11] is used, overruling any information in the data. 
   }
  \item{pars2}{ Contains the model settings \cr \cr
   \code{pars2[1]} corresponds to lx = length of ODE variable x \cr
   \code{pars2[2]} corresponds to ddmodel = diversity-dependent model, model of diversity-dependence, which can be one of\cr \cr
    ddmodel = 0 : no diversity dependence \cr
    ddmodel = 1 : linear dependence in speciation rate \cr
    ddmodel = 11: linear dependence in speciation rate and in immigration rate \cr
    ddmodel = 2 : exponential dependence in speciation rate\cr
    ddmodel = 21: exponential dependence in speciation rate and in immigration rate\cr \cr
   \code{pars2[3]} corresponds to cond = setting of conditioning\cr \cr
    cond = 0 : conditioning on island age \cr
    cond = 1 : conditioning on island age and non-extinction of the island biota \cr \cr
   \code{pars2[4]} sets whether parameters and likelihood should be printed (1) or not (0)
   }
  \item{datalist}{Data object containing information on colonisation and branching times. This object can be generated using the DAISIE_dataprep function, which converts a user-specified data table into a data object, but the object can of course also be entered directly. It is an R list object with the following elements.\cr
  The first element of the list has two or three components: \cr \cr
   \code{$island_age} - the island age \cr
  Then, depending on whether a distinction between types is made, we have:\cr
   \code{$not_present} - the number of mainland lineages that are not present on the island \cr  
  or:\cr
   \code{$not_present_type1} - the number of mainland lineages of type 1 that are not present on the island \cr  
   \code{$not_present_type2} - the number of mainland lineages of type 2 that are not present on the island \cr \cr  
  The remaining elements of the list each contains information on a single colonist lineage on the island and has 5 components:\cr \cr
   \code{$colonist_name} - the name of the species or clade that colonized the island \cr  
   \code{$branching_times} - island age and stem age of the population/species in the case of Non-endemic, Non-endemic_MaxAge and Endemic anagenetic species. For cladogenetic species these should be island age and branching times of the radiation including the stem age of the radiation.\cr                          
   \code{$stac} - the status of the colonist \cr \cr
  * Non_endemic_MaxAge: 1 \cr
  * Endemic: 2 \cr
  * Endemic&Non_Endemic: 3 \cr
  * Non_endemic: 4 \cr
  * Endemic_MaxAge: 5 \cr \cr
   \code{$missing_species} - number of island species that were not sampled for particular clade (only applicable for endemic clades) \cr
   \code{$type1or2} - whether the colonist belongs to type 1 or type 2 \cr
  }
  \item{methode}{Method of the ODE-solver. See package deSolve for details. Default is "lsodes"}

}
\details{
The output is a loglikelihood value
}
\value{The loglikelihood}
\references{
Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium and non-equilibrium dynamics simultaneously operate in the Galapagos islands. Ecology Letters 18: 844-852.}
\author{
Rampal S. Etienne & Bart Haegeman
}
\seealso{
 \code{\link{DAISIE_ML}},
 \code{\link{DAISIE_sim}}
}
\examples{
data(Galapagos_datalist_2types)
pars1 = c(0.195442017,0.087959583,Inf,0.002247364,0.873605049,
          3755.202241,8.909285094,14.99999923,0.002247364,0.873605049,0.163)
pars2 = c(100,11,0,1)
DAISIE_loglik_all(pars1,pars2,Galapagos_datalist_2types)
}
\keyword{models}