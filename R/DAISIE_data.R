#' @name Galapagos_datalist
#' @title Colonization and branching times of 8
#' terrestrial avifaunal clades in list format,
#' accepted by DAISIE_ML and DAISIE_loglik_all
#' @docType data
#' @format A list with 9 elements the first of which contains 2 elements and
#' the following 8 containing 5 components.
#' @description A list containing the colonization and branching times of the
#' terrestrial avifauna in the Galapagos where no distinction is made between
#' types of colonists. This list can be generated using the DAISIE_dataprep
#' function, which converts a user-specified data table into a data object,
#' but the object can of course also be entered directly. It is an R list
#' object with the following elements.\cr \cr
#' The first element of the list has two
#' components: \cr \cr
#' \code{$island_age} - the island age \cr
#' \code{$not_present} - the number of mainland lineages that are not present
#' on the island \cr \cr
#' The following 8 elements of the list each contains
#' information on a single colonist lineage on the island and has 5
#' components:\cr \cr
#' \code{$colonist_name} - the name of the species or clade
#' that colonized the island \cr
#' \code{$branching_times} - island age and stem
#' age of the population/species in the case of Non-endemic, Non-endemic_MaxAge
#' and Endemic anagenetic species. For cladogenetic species these should be
#' island age and branching times of the radiation including the stem age of
#' the radiation.\cr
#' \code{$stac} - the status of the colonist \cr \cr
#' * Non_endemic_MaxAge: 1 \cr
#' * Endemic: 2 \cr
#' * Endemic&Non_Endemic: 3 \cr
#' * Non_endemic: 4 \cr \cr
#' \code{$missing_species} - number of island species
#' that were not sampled for particular clade (only applicable for endemic
#' clades) \cr
#' \code{$type1or2} - whether the colonist belongs to type 1 or
#' type 2. In this dataset all are equal to 1. \cr
#' @seealso \code{\link{DAISIE_dataprep}}, \code{\link{DAISIE_ML}}
#' @source Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium
#' and non-equilibrium dynamics simultaneously operate in the Galapagos
#' islands. Ecology Letters 18: 844-852.
#' @keywords datasets
NULL


#' @name Galapagos_datalist_2types
#' @title Colonization and branching times of 8 terrestrial avifaunal clades in
#' list format, accepted by DAISIE_ML and DAISIE_loglik_all
#' @docType data
#' @format A list with 9 elements the first of which contains 3 elements and
#' the following 8 containing 5 components.
#' @description A list containing the colonization and branching times of the
#' terrestrial avifauna in the Galapagos. This list can be generated using the
#' DAISIE_dataprep function, which converts a user-specified data table into a
#' data object, but the object can of course also be entered directly. It is an
#' R list object with the following elements.\cr \cr
#' The first element of the
#' list has three components: \cr \cr
#' \code{$island_age} - the island age \cr
#' \code{$not_present_type1} - the number of mainland lineages of type 1 that
#' are not present on the island \cr
#' \code{$not_present_type2} - the number of
#' mainland lineages of type 2 that are not present on the island \cr \cr The
#' following 8 elements of the list each contains information on a single
#' colonist lineage on the island and has 5 components:\cr \cr
#' \code{$colonist_name} - the name of the species or clade that colonized the
#' island \cr
#' \code{$branching_times} - island age and stem age of the
#' population/species in the case of Non-endemic, Non-endemic_MaxAge and
#' Endemic anagenetic species. For cladogenetic species these should be island
#' age and branching times of the radiation including the stem age of the
#' radiation.\cr
#' \code{$stac} - the status of the colonist \cr \cr
#' * Non_endemic_MaxAge: 1 \cr
#' * Endemic: 2 \cr
#' * Endemic&Non_Endemic: 3 \cr
#' * Non_endemic: 4 \cr \cr
#' \code{$missing_species} - number of island species
#' that were not sampled for particular clade (only applicable for endemic
#' clades) \cr
#' \code{$type1or2} - whether the colonist belongs to type 1 or
#' type 2. In this dataset only the finches are type 2 \cr
#' @seealso \code{\link{DAISIE_dataprep}}, \code{\link{DAISIE_ML}}
#' @source Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium
#' and non-equilibrium dynamics simultaneously operate in the Galapagos
#' islands. Ecology Letters 18: 844-852.
#' @keywords datasets
NULL


#' @name Galapagos_datatable
#' @title Colonization and branching times of 8 terrestrial avifaunal clades in
#' table format.
#' @docType data
#' @format A table with 8 rows and 4 columns.
#' @description A table containing the colonization and branching times of the
#' terrestrial avifauna in the Galapagos.  Each row on the table represents and
#' independent colonisation event. The table has four columns. \cr \cr
#' \code{$Clade_name} -
#' name of independent colonization event \cr
#' \code{$Status} - One of the
#' following categories: \cr
#' * Non_endemic: for cases where both island and non-island populations of the species have been sampled) \cr
#' * Non_endemic_MaxAge: for cases where island population of the species has not
#' been sampled and only the age of the species is available) \cr
#' * Endemic: applicable for both cladogenetic or anagenetic species \cr
#' * Endemic&Non_Endemic: when endemic clade and mainland ancestor has
#' re-colonized \cr
#' \code{$Missing_species} - Number of island species that
#' were not sampled for particular clade (only applicable for endemic clades)\cr
#' \code{$Branching_times} - Stem age of the population/species in the case
#' of Non-endemic, Non-endemic_MaxAge and Endemic anagenetic species. For
#' cladogenetic species these should be branching times of the radiation
#' including the stem age of the radiation.\cr
#' @seealso \code{\link{DAISIE_dataprep}}, \code{\link{DAISIE_ML}}
#' @source Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium
#' and non-equilibrium dynamics simultaneously operate in the Galapagos
#' islands. Ecology Letters 18: 844-852.
#' @keywords datasets
NULL


#' @name islands_10reps_RAW
#' @title 1000 islands in RAW format simulated with the ML parameters of the CR model
#' for the Galapagos data.
#' @format A list with 10 items.
#' @description Each simulated dataset is an element of the list, which can be called using
#' e.g. islands_10reps_RAW[[1]] Each of the island replicates is a list in
#' itself. The first (e.g. islands_10reps_RAW[[x]][[1]]) element of that list
#' has the following components: \cr The following elements of the RAW list
#' each contain information on a single colonist lineage on the island and has
#' 5 components:\cr
#' \code{$branching_times} - island age and stem age of the
#' population/species in the case of Non-endemic, Non-endemic_MaxAge and
#' Endemic anagenetic species. For cladogenetic species these should be island
#' age and branching times of the radiation including the stem age of the
#' radiation.\cr
#' \code{$stac} - the status of the colonist \cr
#' * Not_present: 0 \cr
#' * Non_endemic_MaxAge: 1 \cr
#' * Endemic: 2 \cr
#' * Endemic&Non_Endemic: 3 \cr
#' * Non_endemic: 4 \cr
#' \code{$stt_table} - Species-through-time table for
#' the descendants of the mainland species (nI - number of non-endemic species;
#' nA - number of anagenetic species, nC - number of cladogenetic species)\cr
#' \code{$missing_species} - number of island species that were not sampled for
#' particular clade (only applicable for endemic clades) \cr
#' @seealso \code{\link{DAISIE_sim_constant_rate}},
#' \code{\link{DAISIE_plot_sims}}
#' @source Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium
#' and non-equilibrium dynamics simultaneously operate in the Galapagos
#' islands. Ecology Letters 18: 844-852.
#' @keywords datasets
NULL


#' @name islands_1type_1000reps
#' @title 1000 islands in DAISIE format simulated with the ML parameters of the CR
#' model for the Galapagos data
#' @format A list with 1000 items.
#' @description Each simulated dataset is an element of the list, which can be called using
#' e.g. islands_1type_1000reps[[1]] Each of the island replicates is a list in
#' itself. The first (e.g. islands_1type_1000reps[[x]][[1]]) element of that
#' list has the following components: \cr
#' \code{$island_age} - the island age\cr
#' \code{$not_present} - the number of mainland lineages that are not
#' present on the island \cr
#' \code{$stt_all} - STT table for all species on the
#' island (nI - number of non-endemic species; nA - number of anagenetic
#' species, nC - number of cladogenetic species, present - number of
#' independent colonisations present )\cr
#' The subsequent elements of the list each contain information on a single
#' colonist lineage on the island and has 3 components:\cr
#' \code{$branching_times} - island age and stem age of the population/species
#' in the case of Non-endemic, Non-endemic_MaxAge and Endemic anagenetic
#' species. For cladogenetic species these should be island age and branching
#' times of the radiation including the stem age of the radiation.\cr
#' \code{$stac} - the status of the colonist \cr
#' * Non_endemic_MaxAge: 1 \cr
#' * Endemic: 2 \cr
#' * Endemic&Non_Endemic: 3 \cr
#' * Non_endemic: 4 \cr
#' \code{$missing_species} - number of island species that were not sampled for
#' particular clade (only applicable for endemic clades) \cr
#' @seealso \code{\link{DAISIE_sim_constant_rate}},
#' \code{\link{DAISIE_plot_sims}}
#' @source Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium
#' and non-equilibrium dynamics simultaneously operate in the Galapagos
#' islands. Ecology Letters 18: 844-852.
#' @keywords datasets
NULL


#' @name islands_2types_1000reps
#' @title 1000 islands in DAISIE format simulated with the ML parameters of the
#' CR_lamc_mu_K model for the Galapagos data (2 types of species)
#' @format A list with 1000 items.
#' @description Each simulated dataset is an element of the list, which can be called using
#' e.g. islands_2types_1000reps[[1]] Each of the island replicates is a list in
#' itself. The first (e.g. islands_2types_1000reps[[x]][[1]]) element of that
#' list has the following components: \cr
#' \code{$island_age} - the island age\cr
#'\code{$not_present_type1} - the number of mainland lineages of type 1
#' that are not present on the island \cr
#' \code{$not_present_type2} - the
#' number of mainland lineages of type 2 that are not present on the island \cr
#' \code{$stt_all} - STT table for all species on the island (nI - number of
#' non-endemic species; nA - number of anagenetic species, nC - number of
#' cladogenetic species, present - number of independent colonisations present
#' )\cr
#' \code{$stt_stt_type1} - STT table for type 1 species on the island (nI
#' - number of non-endemic species; nA - number of anagenetic species, nC -
#' number of cladogenetic species, present - number of independent
#' colonisations present )\cr
#' \code{$stt_stt_type2} - STT table for type 2
#' species on the island (nI - number of non-endemic species; nA - number of
#' anagenetic species, nC - number of cladogenetic species, present - number of
#' independent colonisations present )\cr
#' The subsequent elements of the list each contain information on a single
#' colonist lineage on the island and has 4 components:\cr
#' \code{$branching_times} - island age and stem age of the population/species
#' in the case of Non-endemic, Non-endemic_MaxAge and Endemic anagenetic
#' species. For cladogenetic species these should be island age and branching
#' times of the radiation including the stem age of the radiation.\cr
#' \code{$stac} - the status of the colonist \cr
#' * Non_endemic_MaxAge: 1 \cr
#' * Endemic: 2 \cr
#' * Endemic&Non_Endemic: 3 \cr
#' * Non_endemic: 4 \cr
#' \code{$missing_species} - number of island species that were not sampled for
#' particular clade (only applicable for endemic clades) \cr
#' \code{$type_1or2} - whether the colonist belongs to type 1 or type 2 \cr
#' @seealso \code{\link{DAISIE_sim_constant_rate}},
#' \code{\link{DAISIE_plot_sims}}
#' @source Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium
#' and non-equilibrium dynamics simultaneously operate in the Galapagos
#' islands. Ecology Letters 18: 844-852.
#' @keywords datasets
NULL


#' @name Macaronesia_datalist
#' @title Colonization and branching times of terrestrial avifaunal clades from
#' Azores, Canary Islands, Cape Verde and Madeira in list format, accepted by
#' DAISIE_ML and DAISIE_loglik_all
#' @docType data
#' @format A list with 4 main elements for each archipelago. Each element has
#' several sub-elements.
#' @description A list containing the colonization and branching times of the terrestrial
#' avifauna in 4 archipelagos: Azores, Canary Islands, Cape Verde and Madeira.
#' It is an R list object with the 4 main elements corresponding to each of the
#' archipelagos (e.g. Macaronesia_datalist[[1]] calls the Azores data). Each of
#' the four elements is then made of several elemants:\cr \cr
#' The first element of the list for an archipelago has two components: \cr \cr
#' \code{$island_age} - the island age \cr
#' \code{$not_present} - the number of
#' mainland lineages that are not present on the island \cr \cr
#' The following elements of the list each contains information on a single colonist lineage
#' on the island and has 5 components:\cr \cr
#' \code{$colonist_name} - the name
#' of the species or clade that colonized the island \cr
#' \code{$branching_times} - island age and stem age of the population/species
#' in the case of Non-endemic, Non-endemic_MaxAge, Endemic_MaxAge and Endemic
#' anagenetic species. For cladogenetic species the island age and branching
#' times of the radiation including the stem age of the radiation are shown.\cr
#' \code{$stac} - the status of the colonist \cr \cr
#' * Non_endemic_MaxAge: 1\cr
#' * Endemic: 2 \cr
#' * Endemic&Non_Endemic: 3 \cr
#' * Non_endemic: 4 \cr \cr
#' * Endemic_MaxAge: 5 \cr \cr
#' #' \code{$missing_species} - number of island species that were not sampled for
#' particular clade (only applicable for endemic clades) \cr
#' \code{$type1or2} -
#' whether the colonist belongs to type 1 or type 2. In this dataset all are
#' equal to 1. \cr
#' @seealso \code{\link{DAISIE_dataprep}}, \code{\link{DAISIE_ML}}
#' @source
#' Valente L., Illera J.C, Havenstein K., Pallien T., Etienne R.S., Tiedemann
#' R. Macroevolutionary dynamics in Atlantic island avifaunas: were MacArthur &
#' Wilson right about equilibrium? Under review.
#' @keywords datasets
NULL


#' @name Biwa_datalist
#' @title Colonization and branching times of 68 fish clades in list
#' format, accepted by DAISIE_ML, DAISIE_SR_ML, DAISIE_loglik_all and DAISIE_SR_loglik_all
#' @docType data
#' @format A list with 69 elements the first of which contains 2 elements and
#' the following 68 containing 5 components.
#' #' @description A list containing the colonization and branching times of the
#' Lake Biwa (Japan) fishes where no distinction is made between types of
#' colonists. This list can be generated using the DAISIE_dataprep function,
#' which converts a user-specified data table into a data object, but the
#' object can of course also be entered directly. It is an R list object with
#' the following elements.\cr \cr
#' The first element of the list has two
#' components: \cr \cr
#' \code{$island_age} - the island age \cr
#' \code{$not_present} - the number of mainland lineages that are not present
#' on the island \cr \cr
#' The following 8 elements of the list each contains
#' information on a single colonist lineage on the island and has 5
#' components:\cr \cr
#' \code{$colonist_name} - the name of the species or clade
#' that colonized the island \cr
#' \code{$branching_times} - island age and stem
#' age of the population/species in the case of Non-endemic, Non-endemic_MaxAge
#' and Endemic anagenetic species. For cladogenetic species these should be
#' island age and branching times of the radiation including the stem age of
#' the radiation.\cr
#' \code{$stac} - the status of the colonist \cr \cr
#' * Non_endemic_MaxAge: 1 \cr
#' * Endemic: 2 \cr
#' * Endemic&Non_Endemic: 3 \cr
#' * Non_endemic: 4 \cr \cr
#' \code{$missing_species} - number of island species
#' that were not sampled for particular clade (only applicable for endemic
#' clades) \cr
#' \code{$type1or2} - whether the colonist belongs to type 1 or
#' type 2. In this dataset all are equal to 1. \cr
#' @seealso \code{\link{DAISIE_dataprep}}, \code{\link{DAISIE_ML}}, \code{\link{DAISIE_SR_ML}}
#' @source Hauffe, T., D. Delicado, R.S. Etienne and L. Valente (submitted).
#' Lake expansion increases equilibrium diversity via the target effect of
#' island biogeography
#' @keywords datasets
NULL
