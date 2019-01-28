\name{DAISIE-package}
\alias{DAISIE-package}
\alias{DAISIE}
\docType{package}
\title{
Dynamical Assembly of Islands by Speciation, Immigration and Extinction
}
\description{
DAISIE is an island biogeography model that assumes island biota assembly is governed by immigration, extinction and speciation through cladogenesis and anagenesis. This likelihood-based statistical package can simulate islands and estimate parameters of the DAISIE model based on phylogenetic/phylogeographic data. Cladogenesis and immigration rates can be dependent on diversity.
}
\details{
\tabular{ll}{
Package: \tab DAISIE\cr
Type: \tab Package\cr
Version: \tab 1.4\cr
Date: \tab 2017-04-12\cr
License: \tab GPL 2.0\cr
}

DAISIE_ML computes the parameters that maximizes the loglikelihood computed by DAISIE_loglik_all for one island/archipelago

DAISIE_ML2 computes the parameters that maximizes the loglikelihood for a set of islands / archipelagoes

DAISIE_sim simulates islands with given parameters

DAISIE_tutorial opens a step-by-step tutorial with example from Galapagos islands

DAISIE_loglik_all computes the loglikelihood of a given parameterization of the DAISIE model for a data set of colonization and branching times

DAISIE_plot_sims produces species-through-time plots based on simulated islands

DAISIE_dataprep converts data in table format to the list format accepted by DAISIE_loglik_all and DAISIE_ML

DAISIE_ExpEIN computes the expected number of endemics, non-endemics and their sum for a given time and set of parameters

DAISIE_probdist calculates the joint distribution of the number of endemics and non-endemics for a given set of parameter values, a given mainland species pool size and a given set of times

DAISIE_margprobdist calculates the marginal distribution of the number of endemics and non-endemics and their sum for a given set of parameter values, a given mainland species pool size and a given set of times

DAISIE_numcol calculates expectation and marginal distribution of the number of colonizations (lineages) for a given set of parameter values, a given mainland species pool size and a given set of times
}
\author{
Rampal S. Etienne, Luis M. Valente, Albert B. Phillimore & Bart Haegeman\cr
Maintainer: Rampal S. Etienne <r.s.etienne@rug.nl>rep(0:(D1 - 1), each = D2)
}
\references{
Valente, L., A.B. Phillimore & R.S. Etienne (2015). Equilibrium and non-equilibrium dynamics simultaneously operate in the Galapagos islands. Ecology Letters 18: 844-852, <DOI:10.1111/ele.12461>.

Valente, L., R.S. Etienne & L. Davalos (2017). Recent extinctions disturb path to equilibrium diversity in Caribbean bats. Nature Ecology and Evolution 1: 0026. <DOI:10.1038/s41559-016-0026>.
}
\keyword{ models }
\seealso{
\code{DDD}
}

