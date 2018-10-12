\name{DAISIE_dataprep}
\alias{DAISIE_dataprep}
\title{ Prepare colonisation and branching time data to run in DAISIE.}
\description{
   This function produces a data object that can be run in DAISIE likelihood 
computation/optimization functions. The function converts a user-specified table
to a DAISIE-compatible format. See Galapagos_datatable.Rdata for a template of an input table.)
}
\usage{
DAISIE_dataprep(
	datatable,
	island_age,
	M,
	number_clade_types = 1,
	list_type2_clades = NA,
	prop_type2_pool = "proportional",
	epss = 1E-5
	)
}

\arguments{
	\item{datatable}{Data frame (table) with user-specified data. See file Galapagos_datatable.Rdata 
for a template of an input table. Each row on the table represents and independent colonisation event. Table has the following four columns. \cr \cr
  \code{$Clade_name} - name of independent colonization event \cr
  \code{$Status} - One of the following categories: \cr
  * "Non_endemic": applies to non-endemic species for cases where both island and non-island populations of the species have been sampled) \cr
  * "Non_endemic_MaxAge": applies to non-endemic species for cases where island individuals of the species have not been sampled and only the age of the species is available) \cr
  * "Endemic": applies to endemic species and is applicable for both cladogenetic and anagenetic species \cr
  * "Endemic_MaxAge": applies to endemic species for cases where island individuals of the species have not been sampled and only the age of the species is available. This could apply to endemic species that have recently gone extinct because of antropogenic causes that are (evidently) not modelled, and for which no DNA data is available.\cr
  * "Endemic&Non_Endemic": when endemic clade is present and its mainland ancestor has re-colonized \cr
  \code{$Missing_species} - Number of island species that were not sampled for particular clade (only applicable for "Endemic" clades) \cr
  \code{$Branching_times} - Stem age of the population/species in the case of "Non-endemic", "Non-endemic_MaxAge" and "Endemic" anagenetic species. For "Endemic" cladogenetic 
species these should be branching times of the radiation including the stem age of the radiation.\cr
}
	\item{island_age}{Age of island in appropriate units}
	\item{M}{The size of the mainland pool, i.e the number of species that can potentially colonize the island}
	\item{number_clade_types}{Number of clade types. Default: number_clade_types = 1 all species 
are considered to belong to same macroevolutionary process. 
If number_clade_types = 2, there are two types of clades with distinct macroevolutionary processes.}
	\item{list_type2_clades}{If number_clade_types = 2, list_type2_clades specifies the 
names of the clades that have a distinct macroevolutionary process. The names must match 
those in the $Clade_name column of the source data table (e.g. list_type2_clades = "Finches"). 
If number_clade_types = 1, then list_type2_clades = NA should be specified (default)}
	\item{prop_type2_pool}{Specifies the fraction of potential mainland colonists that 
have a distinct macroevolutionary process. Applies only if number_clade_types = 2. 
Default "proportional" sets the fraction to be proportional to the number of clades of distinct
 macroevolutionary process that have colonised the island. Alternatively, the user
 can specify a value between 0 and 1 (e.g. if mainland pool size is 1000 and prop_type2_pool = 0.02 then number of type2 species is 20). }
	\item{epss}{ Default= 1E-5 should be appropriate in most cases. This value is used to set the maximum age of colonisation of "Non_endemic_MaxAge" and "Endemic_MaxAge" species to an age that is slightly younger than the island for cases when the age provided for that species is older than the island. The new maximum age is then used as an upper bound to integrate over all possible colonisation times. }
}
\details{The output is an R list containing the data formatted to be run
 on other DAISIE functions.}
\value{
  \item{datalist}{
  R list object containing data:\cr
  The first element of the list has two or three components: \cr
  \code{$island_age} - the island age \cr
  Then, depending on whether a distinction between species types is made, we have:\cr
  \code{$not_present} - the number of mainland lineages that are not present on the island \cr  
  or:\cr
  \code{$not_present_type1} - the number of mainland lineages of type 1 that are not present on the island \cr  
  \code{$not_present_type2} - the number of mainland lineages of type 2 that are not present on the island \cr  
  The following elements of the list each contains information on a single colonist lineage on the island and has 5 components:\cr
  \code{$colonist_name} - the name of the species or clade that colonized the island \cr  
  \code{$branching_times} - island age and stem age of the population/species in the case of "Non-endemic", "Non-endemic_MaxAge" and "Endemic" anagenetic species. For "Endemic" cladogenetic species these are island age and branching times of the radiation including the stem age of the radiation.\cr                          
  \code{$stac} - the status of the colonist \cr
  * Non_endemic_MaxAge: 1 \cr
  * Endemic: 2 \cr
  * Endemic&Non_Endemic: 3 \cr
  * Non_endemic: 4 \cr
  * Endemic_MaxAge: 5 \cr
  \code{$missing_species} - number of island species that were not sampled for particular clade (only applicable for endemic clades) \cr
  \code{$type_1or2} - whether the colonist belongs to type 1 or type 2 \cr
  }
}
\references{
Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium and non-equilibrium dynamics simultaneously operate in the Galapagos islands. Ecology Letters 18: 844-852.}
\author{
Luis M Valente
}
\examples{ 



### Create Galapagos data object where all taxa have the same macroevolutionary process

data(Galapagos_datatable)
DAISIE_dataprep(
   datatable = Galapagos_datatable,
   island_age = 4,
   M = 1000
   )

### Create Galapagos data object with a distinct macroevolutionary processes
# for the Darwin's finches. One process applies to type 1 species (all species 
# except for Darwin's finches) and the other applies only to type 2 species 
# (Darwin's finches). Set fraction of potential colonists of type 2 to be
# proportional to the number of type2 clades present on the island.

data(Galapagos_datatable)
DAISIE_dataprep(
   datatable = Galapagos_datatable,
   island_age = 4,
   M = 1000,
   number_clade_types = 2,
   list_type2_clades = "Finches"
   )

### Create Galapagos data object with a distinct macroevolutionary processes 
# for the Darwin's finches. One process applies to type 1 species (all species
# except for Darwin's finches) and the other applies only to type 2 species 
# (Darwin's finches). Set fraction of potential colonists of type 2 to be 0.163.

data(Galapagos_datatable)
DAISIE_dataprep(
   datatable = Galapagos_datatable,
   island_age = 4,
   M = 1000,
   number_clade_types = 2,
   list_type2_clades = "Finches",
   prop_type2_pool = 0.163
   )
}
\keyword{ models }