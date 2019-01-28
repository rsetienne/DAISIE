\name{DAISIE_ML}
\alias{DAISIE_ML}
\title{Maximization of the loglikelihood under the DAISIE model}
\description{
This function computes the maximum likelihood estimates of the parameters of the
DAISIE model for data from lineages colonizing an island. It also outputs the
corresponding loglikelihood that can be used in model comparisons.
}
\usage{
DAISIE_ML(  
   datalist,
   datatype = 'single',
   initparsopt,
   idparsopt,
   parsfix,
   idparsfix,
   idparsnoshift = 6:10,
   idparsmat = NULL,
   res = 100,
   ddmodel = 0,
   cond = 0,
   eqmodel = 0,
   x_E = 0.95,
   x_I = 0.98,
   tol = c(1e-04, 1e-05, 1e-07),
   maxiter = 1000 * round((1.25)^length(idparsopt)),
   methode = 'lsodes',
   optimmethod = 'subplex'
   )
}
\arguments{
  \item{datalist}{Data object containing information on colonisation and branching times. This object can be generated using the DAISIE_dataprep function, which converts a user-specified data table into a data object, but the object can of course also be entered directly. It is an R list object with the following elements.\cr
  The first element of the list has two three components: \cr \cr
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
  * Non_endemic: 4 \cr \cr
  \code{$missing_species} - number of island species that were not sampled for particular clade (only applicable for endemic clades) \cr
  \code{$type1or2} - whether the colonist belongs to type 1 or type 2 \cr
  }
  \item{datatype}{Sets the type of data: 'single' for a single island or archipelago treated as one, and 'multiple' for multiple archipelagoes potentially sharing the same parameters} 
  \item{initparsopt}{The initial values of the parameters that must be optimized}
  \item{idparsopt}{
   The ids of the parameters that must be optimized. The ids are defined as follows: \cr \cr
   id = 1 corresponds to lambda^c (cladogenesis rate) \cr
   id = 2 corresponds to mu (extinction rate) \cr
   id = 3 corresponds to K (clade-level carrying capacity) \cr
   id = 4 corresponds to gamma (immigration rate) \cr
   id = 5 corresponds to lambda^a (anagenesis rate) \cr   
   id = 6 corresponds to lambda^c (cladogenesis rate) for an optional subset of the species \cr
   id = 7 corresponds to mu (extinction rate) for an optional subset of the species\cr
   id = 8 corresponds to K (clade-level carrying capacity) for an optional subset of the species\cr
   id = 9 corresponds to gamma (immigration rate) for an optional subset of the species\cr
   id = 10 corresponds to lambda^a (anagenesis rate) for an optional subset of the species\cr   
   id = 11 corresponds to p_f (fraction of mainland species that belongs to the second subset of species
   }
  \item{idparsfix}{The ids of the parameters that should not be optimized, e.g. c(1,3) if lambda^c and K should not be optimized.}
  \item{parsfix}{The values of the parameters that should not be optimized }
  \item{idparsnoshift}{For datatype = 'single' only: The ids of the parameters that should not be different between two groups of species; 
This can only apply to ids 6:10, e.g. idparsnoshift = c(6,7) means that lambda^c and mu have the same values for both groups}\cr
  \item{idparsmat}{For datatype = 'multiple' only: Matrix containing the ids of the parameters, linking them to initparsopt and parsfix. Per island system we use the following order: \cr \cr
  * lac = (initial) cladogenesis rate \cr
  * mu = extinction rate \cr
  * K = maximum number of species possible in the clade \cr
  * gam = (initial) immigration rate \cr
  * laa = (initial) anagenesis rate \cr
  Example: idparsmat = rbind(c(1,2,3,4,5),c(1,2,3,6,7)) has different rates of immigration and anagenesis for the two islands.} 
  \item{res}{Sets the maximum number of species for which a probability must be computed, must be larger than the size of the largest clade }
  \item{ddmodel}{Sets the model of diversity-dependence: \cr \cr
    ddmodel = 0 : no diversity dependence \cr
    ddmodel = 1 : linear dependence in speciation rate \cr
    ddmodel = 11: linear dependence in speciation rate and in immigration rate \cr
    ddmodel = 2 : exponential dependence in speciation rate\cr
    ddmodel = 21: exponential dependence in speciation rate and in immigration rate\cr
  }
  \item{cond}{
    cond = 0 : conditioning on island age \cr
    cond = 1 : conditioning on island age and non-extinction of the island biota \cr
  }
  \item{eqmodel}{Sets the equilibrium constraint that can be used during the likelihood optimization. Only available for datatype = 'single'.\cr\cr
    eqmodel = 0 : no equilibrium is assumed \cr
    eqmodel = 13 : near-equilibrium is assumed on endemics using deterministic equation for endemics and immigrants. Endemics must be within x_E of the equilibrium value\cr
    eqmodel = 15 : near-equilibrium is assumed on endemics and immigrants using deterministic equation for endemics and immigrants. Endemics must be within x_E of the equilibrium value, while non-endemics must be within x_I of the equilibrium value
  }
  \item{x_E}{ Sets the fraction of the equlibrium endemic diversity above which the endemics are assumed to be in equilibrium; only active for eqmodel = 13 or 15 }
  \item{x_I}{ Sets the fraction of the equlibrium non-endemic diversity above which the system is assumed to be in equilibrium; only active for eqmodel = 15 }
  \item{tol}{Sets the tolerances in the optimization. Consists of: \cr
   reltolx = relative tolerance of parameter values in optimization \cr
   reltolf = relative tolerance of function value in optimization \cr
   abstolx = absolute tolerance of parameter values in optimization
  }
  \item{maxiter}{Sets the maximum number of iterations in the optimization}
  \item{methode}{Method of the ODE-solver. See package deSolve for details. Default is "lsodes"}
  \item{optimmethod}{Method used in likelihood optimization. Default is "subplex" (see subplex package). Alternative is 'simplex' which was the method in previous versions.}  
}
\details{The result of sort(c(idparsopt, idparsfix, idparsnoshift)) should be identical to c(1:10). If not, an error is reported that the input is incoherent. The same happens when the length of initparsopt is different from the length of idparsopt, and the length of parsfix is different from the length of idparsfix.\cr
Including the 11th parameter (p_f) in either idparsopt or idparsfix (and therefore initparsopt or parsfix) is optional. If this parameter is not specified, then the information in the data is used, otherwise the information in the data is overruled.
}
\value{The output is a dataframe containing estimated parameters and maximum loglikelihood.
  \item{lambda_c}{ gives the maximum likelihood estimate of lambda^c, the rate of cladogenesis}
  \item{mu}{ gives the maximum likelihood estimate of mu, the extinction rate}
  \item{K}{ gives the maximum likelihood estimate of K, the carrying-capacity}
  \item{gamma}{ gives the maximum likelihood estimate of gamma, the immigration rate }
  \item{lambda_a}{ gives the maximum likelihood estimate of lambda^a, the rate of anagenesis}
  \item{lambda_c2}{ gives the maximum likelihood estimate of lambda^c2, the rate of cladogenesis for the optional second group of species}
  \item{mu2}{ gives the maximum likelihood estimate of mu2, the extinction rate for the optional second group of species}
  \item{K2}{ gives the maximum likelihood estimate of K2, the carrying-capacity for the optional second group of species}
  \item{gamma2}{ gives the maximum likelihood estimate of gamma2, the immigration rate for the optional second group of species}
  \item{lambda_a2}{ gives the maximum likelihood estimate of lambda^a2, the rate of anagenesis for the optional second group of species}
  \item{loglik}{ gives the maximum loglikelihood}
  \item{df}{ gives the number of estimated parameters, i.e. degrees of feedom}
  \item{conv}{ gives a message on convergence of optimization; conv = 0 means convergence}
}
\references{
Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium and non-equilibrium dynamics simultaneously operate in the Galapagos islands. Ecology Letters 18: 844-852. <DOI:10.1111/ele.12461>.}
\author{
Rampal S. Etienne
}
\seealso{
 \code{\link{DAISIE_loglik_all}},
 \code{\link{DAISIE_sim}}
}
\examples{
cat("
### When all species have the same rates, and we want to optimize all 5 parameters,
# we use:

data(Galapagos_datalist)
DAISIE_ML(
   datalist = Galapagos_datalist,
   initparsopt = c(2.5,2.7,20,0.009,1.01),
   ddmodel = 11,
   idparsopt = 1:5,
   parsfix = NULL,
   idparsfix = NULL
)

### When all species have the same rates, and we want to optimize all parameters 
# except K (which we set equal to Inf), we use:

data(Galapagos_datalist)
DAISIE_ML(
   datalist = Galapagos_datalist,
   initparsopt = c(2.5,2.7,0.009,1.01),
   idparsopt = c(1,2,4,5),
   parsfix = Inf,
   idparsfix = 3
   )

### When all species have the same rates except that the finches have a different
# rate of cladogenesis, and we want to optimize all parameters except K (which we
# set equal to Inf), fixing the proportion of finch-type species at 0.163, we use:

data(Galapagos_datalist_2types)
DAISIE_ML(
   datalist = Galapagos_datalist_2types,
   initparsopt = c(0.38,0.55,0.004,1.1,2.28),
   idparsopt = c(1,2,4,5,6),
   parsfix = c(Inf,Inf,0.163),
   idparsfix = c(3,8,11),
   idparsnoshift = c(7,9,10)
   )

### When all species have the same rates except that the finches have a different
# rate of cladogenesis, extinction and a different K, and we want to optimize all
# parameters, fixing the proportion of finch-type species at 0.163, we use:

data(Galapagos_datalist_2types)
DAISIE_ML(
   datalist = Galapagos_datalist_2types,
   ddmodel = 11,   
   initparsopt = c(0.19,0.09,0.002,0.87,20,8.9,15),
   idparsopt = c(1,2,4,5,6,7,8),
   parsfix = c(Inf,0.163),
   idparsfix = c(3,11),
   idparsnoshift = c(9,10)
   )


### When all species have the same rates except that the finches have a different
# rate of extinction, and we want to optimize all parameters except K (which we 
# set equal to Inf), and we also# want to estimate the fraction of finch species
# in the mainland pool. we use:

data(Galapagos_datalist_2types)
DAISIE_ML(
   datalist = Galapagos_datalist_2types,
   initparsopt = c(2.48,2.7,0.009,1.01,2.25,0.163),
   idparsopt = c(1,2,4,5,7,11),
   parsfix = c(Inf,Inf),
   idparsfix = c(3,8),
   idparsnoshift = c(6,9,10)
   )

### When we have two islands with the same rates except for immigration and anagenesis rate,
# and we want to optimize all parameters, we use:

data(Galapagos_datalist)
DAISIE_ML(
   datalist = list(Galapagos_datalist,Galapagos_datalist),
   datatype = 'multiple',
   initparsopt = c(2.5,2.7,20,0.009,1.01,0.009,1.01),
   idparsmat = rbind(1:5,c(1:3,6,7)),
   idparsopt = 1:7,
   parsfix = NULL,
   idparsfix = NULL
)

### When we consider the four Macaronesia archipelagoes and set all parameters the same
# except for rates of cladogenesis, extinction and immigration for Canary Islands,
# rate of cladogenesis is fixed to 0 for the other archipelagoes,
# diversity-dependence is assumed to be absent
# and we want to optimize all parameters, we use:

data(Macaronesia_datalist)
DAISIE_ML(
   datalist = Macaronesia_datalist,
   datatype = 'multiple',
   initparsopt = c(1.053151832,0.052148979,0.512939011,0.133766934,0.152763179),
   idparsmat = rbind(1:5,c(6,2,3,7,5),1:5,1:5),
   idparsopt = c(2,4,5,6,7),
   parsfix = c(0,Inf),
   idparsfix = c(1,3)
)
   
")
}
\keyword{models}