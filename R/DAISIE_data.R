#' @name Galapagos_datalist
#' @title Colonization and branching times of 8 terrestrial avifaunal clades in list
#' format, accepted by DAISIE_ML and DAISIE_loglik_all
#' @docType data
#' @format A list with 9 elements the first of which contains 2 elements and
#' the following 8 containing 5 components.
#' @description A list containing the colonization and branching times of the terrestrial
#' avifauna in the Galapagos where no distinction is made between types of
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
#' \code{$branching_times} - island age followed by stem age of the population/species
#' in the case of Non-endemic, Non-endemic_MaxAge species and Endemic species with no close
#' relatives on the island. For endemic clades with more than one species on the island
#' (cladogenetic clades/ radiations) these should be island age followed by the
#' branching times of the island cladeincluding the stem age of the clade.\cr
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
#' @title Colonization and branching times of 8 terrestrial avifaunal clades in list
#' format, accepted by DAISIE_ML and DAISIE_loglik_all
#' @docType data
#' @format A list with 9 elements the first of which contains 3 elements and
#' the following 8 containing 5 components.
#' @description A list containing the colonization and branching times of the terrestrial
#' avifauna in the Galapagos. This list can be generated using the
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
#' \code{$branching_times} - island age followed by stem age of the population/species
#' in the case of Non-endemic, Non-endemic_MaxAge species and Endemic species with no close relatives
#' on the island. For endemic clades with more than one species on the island
#' (cladogenetic clades/ radiations) these should be island age followed by the
#' branching times of the island cladeincluding the stem age of the clade.\cr
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
#' @title Colonization and branching times of 8 terrestrial avifaunal Galápagos clades in table
#' format.
#' @docType data
#' @format A table with 8 rows and 4 columns.
#' @description A table containing the colonization and branching times of the terrestrial
#' avifauna in the Galápagos. Each row on the table represents and independent
#' colonisation event. The table has four columns. \cr \cr
#' \code{$Clade_name} -
#' name of independent colonization event \cr
#' \code{$Status} - One of the
#' following categories: \cr
#' * Non_endemic: for non-endemic island species when an approximate time of colonisation is known \cr
#' * Non_endemic_MaxAge: for non-endemic island species when colonisation time is unknown \cr
#' * Endemic: for endemic species when an approximate colonisation time is known \cr
#' * "Endemic_MaxAge": applies to endemic species or endemic clades for cases where the
#' colonisation time is unknown, or when
#' the user wants to specify an upper bound for colonisation.
#' This could for example apply to endemic species that have recently gone extinct because
#' of anthropogenic causes, and which are not included
#' in the phylogeny ("NA" should be given in the branching times column). It
#' could also apply to insular radiations with long stem branches, for which the
#' time of the first cladogenetic event is known, but the precise time of colonisation
#' is not.\cr
#' * Endemic&Non_Endemic: when endemic clade and mainland ancestor has
#' re-colonized \cr
#'
#' \code{$Missing_species} - Number of island species that
#' were not sampled for particular clade (only applicable for endemic clades)\cr
#' \code{$Branching_times} - Stem age of the population/species in the case of "Non_endemic",
#'  "Non_endemic_MaxAge" and "Endemic" species with no extant close relatives on the island.
#'  Set "NA" if colonisation time unknown and no upper bound is known.
#' For "Endemic" cladogenetic species these should be branching times of the
#' radiation, including the stem age of the radiation (colonisation time estimate).\cr
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
#' @seealso \code{\link{DAISIE_sim_cr}()}, \code{\link{DAISIE_plot_sims}}
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
#' @seealso \code{\link{DAISIE_sim_cr}()}, \code{\link{DAISIE_plot_sims}}
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
#' @seealso \code{\link{DAISIE_sim_cr}()}, \code{\link{DAISIE_plot_sims}}
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
#' \code{$branching_times} - island age followed by stem age of the population/species
#' in the case of Non-endemic, Non-endemic_MaxAge species and Endemic species with no close
#' relatives on the island. For endemic clades with more than one species on the island
#' (cladogenetic clades/ radiations) these should be island age followed by the
#' branching times of the island clade including the stem age of the clade. \cr \cr
#' * Non_endemic_MaxAge: 1\cr
#' * Endemic: 2 \cr
#' * Endemic&Non_Endemic: 3 \cr
#' * Non_endemic: 4 \cr \cr
#' * Endemic_MaxAge: 5 \cr \cr
#' \code{$missing_species} - number of island species that were not sampled for
#' particular clade (only applicable for endemic clades) \cr
#' \code{$type1or2} -
#' whether the colonist belongs to type 1 or type 2. In this dataset all are
#' equal to 1. \cr
#' @seealso \code{\link{DAISIE_dataprep}}, \code{\link{DAISIE_ML}}
#' @source
#' Valente L., Illera J.C, Havenstein K., Pallien T., Etienne R.S., Tiedemann
#' R. Equilibrium bird species diversity in Atlantic islands. 2017 Current Biology, 27, 1660-1666.
#' @keywords datasets
NULL



#' @name NewZealand_birds_datalist
#' @title Colonization and branching times of New Zealand birds.
#' @format A list with 40 elements, the first of which contains 2 elements and
#' the following 39 containing 5 components.
#' @description A list containing the colonization and branching times of the birds of
#' New Zealand. Main dataset used in
#' Valente, Etienne, Garcia-R (2019) Current Biology. Island age 52 Myr and mainland
#' pool size of 1000 species. \cr
#' The first element of the list has two
#' components: \cr \cr
#' \code{$island_age} - the island age \cr
#' \code{$not_present} - the number of mainland lineages that are not present
#' on the island \cr \cr
#' The following elements of the list each contain
#' information on a single colonist lineage on the island and has 5
#' components:\cr \cr
#' \code{$colonist_name} - the name of the species or clade
#' that colonized the island \cr
#' \code{$branching_times} - island age followed by stem age of the population/species
#' in the case of Non-endemic, Non-endemic_MaxAge species and Endemic species with no close
#' relatives on the island. For endemic clades with more than one species on the island
#' (cladogenetic clades/ radiations) these should be island age followed by the
#' branching times of the island clade including the stem age of the clade.\cr
#' \code{$stac} - the status of the colonist \cr \cr
#' * Non_endemic_MaxAge: 1 \cr
#' * Endemic: 2 \cr
#' * Endemic&Non_Endemic: 3 \cr
#' * Non_endemic: 4 \cr
#' * Endemic_MaxAge: 5 or 6 \cr \cr
#' \code{$missing_species} - number of island species
#' that were not sampled for particular clade (only applicable for endemic
#' clades) \cr
#' \code{$type1or2} - whether the colonist belongs to type 1 or
#' type 2. In this dataset all are equal to 1. \cr
#' @seealso \code{\link{DAISIE_dataprep}}, \code{\link{DAISIE_ML}}, \code{\link{DAISIE_SR_ML}}
#' @source Valente L, Etienne RS, Garcia-R JC (2019) Deep Macroevolutionary Impact of
#' Humans on New Zealand’s
#' Unique Avifauna. Current Biology, 29, 2563–2569.\cr
#' @keywords datasets
NULL



#' @name Bats_GreaterAntilles
#' @title Colonization and branching times of Noctilionoid bats from the Greater Antilles.
#' @format A list with 17 elements, the first of which contains 2 elements and
#' the following 16 containing 5 components.
#' @description A list containing the colonization and branching times of the noctilionoid
#' bats of the Greater Antilles. Main dataset used in
#' Valente, Etienne and Dávalos (2017) Nature Ecology and Evolution. Island age 20 Myr and
#' mainland pool size of 100 species.\cr
#' The first element of the list has two
#' components: \cr \cr
#' \code{$island_age} - the island age \cr
#' \code{$not_present} - the number of mainland lineages that are not present
#' on the island \cr \cr
#' The following elements of the list each contain
#' information on a single colonist lineage on the island and has 5
#' components:\cr \cr
#' \code{$colonist_name} - the name of the species or clade
#' that colonized the island \cr
#' \code{$branching_times} - island age followed by stem age of the population/species
#' in the case of Non-endemic, Non-endemic_MaxAge species and Endemic species with no close
#' relatives on the island. For endemic clades with more than one species on the island
#' (cladogenetic clades/ radiations) these should be island age followed by the
#' branching times of the island clade including the stem age of the clade.\cr
#' \code{$stac} - the status of the colonist \cr \cr
#' * Non_endemic_MaxAge: 1 \cr
#' * Endemic: 2 \cr
#' * Endemic&Non_Endemic: 3 \cr
#' * Non_endemic: 4 \cr
#' * Endemic_MaxAge: 5 or 6 \cr \cr
#' \code{$missing_species} - number of island species
#' that were not sampled for particular clade (only applicable for endemic
#' clades) \cr
#' \code{$type1or2} - whether the colonist belongs to type 1 or
#' type 2. In this dataset all are equal to 1. \cr
#' @seealso \code{\link{DAISIE_dataprep}}, \code{\link{DAISIE_ML}}, \code{\link{DAISIE_SR_ML}}
#' @source Valente L, Etienne RS, Dávalos (2017) Recent extinctions disturb path to
#' equilibrium diversity in Caribbean bats. Nature Ecology and Evolution, 1, 26.\cr
#' @keywords datasets
NULL




#' @name Biwa_datalist
#' @title Colonization and branching times of 68 fish clades.
#' @docType data
#' @format A list with 69 elements, the first of which contains 2 elements and
#' the following 68 containing 5 components.
#' @description A list containing the colonization and branching times of the fishes of
#' Lake Biwa (Japan). Main dataset used in Hauffe et al (2020). This list can be generated
#' using the DAISIE_dataprep function,
#' which converts a user-specified data table into a data object, but the
#' object can of course also be entered directly. It is an R list object with
#' the following elements.\cr \cr
#' The first element of the list has two
#' components: \cr \cr
#' \code{$island_age} - the island age \cr
#' \code{$not_present} - the number of mainland lineages that are not present
#' on the island \cr \cr
#' The following elements of the list each contains
#' information on a single colonist lineage on the island and has 5
#' components:\cr \cr
#' \code{$colonist_name} - the name of the species or clade
#' that colonized the island \cr
#' \code{$branching_times} - island age followed by stem age of the population/species
#' in the case of Non-endemic, Non-endemic_MaxAge species and Endemic species with no close
#' relatives on the island. For endemic clades with more than one species on the island
#' (cladogenetic clades/ radiations) these should be island age followed by the
#' branching times of the island clade including the stem age of the clade.\cr
#' \code{$stac} - the status of the colonist \cr \cr
#' * Non_endemic_MaxAge: 1 \cr
#' * Endemic: 2 \cr
#' * Endemic&Non_Endemic: 3 \cr
#' * Non_endemic: 4 \cr
#' * Endemic_MaxAge: 5 or 6 \cr \cr
#' \code{$missing_species} - number of island species
#' that were not sampled for particular clade (only applicable for endemic
#' clades) \cr
#' \code{$type1or2} - whether the colonist belongs to type 1 or
#' type 2. In this dataset all are equal to 1. \cr
#' @seealso \code{\link{DAISIE_dataprep}}, \code{\link{DAISIE_ML}}, \code{\link{DAISIE_SR_ML}}
#' @source Hauffe, T., D. Delicado, R.S. Etienne and L. Valente. Lake expansion
#' elevates equilibrium diversity via increasing colonisation. (2020) Journal of Biogeography \cr
#' @keywords datasets
NULL

#' @name archipelago_data
#' @title Physical data on 41 archipelagos
#' @docType data
#' @format A dataframe containing information on archipelago name, area, age and distance from the mainland
#' @description A dataframe with in subsequent columns the name of the archipelago (Archipelago)
#' the area of the archipelago (Area), the age (Age) and the distance from the mainland (Distance)
#' @source Valente L, Phillimore AB, Melo M, Warren BH, Clegg SM, Havenstein K,
#'  Tiedemann R, Illera JC, Thébaud C, Aschenbach T, Etienne RS. A simple dynamic model
#'  explains island bird diversity worldwide (2020) Nature, 579, 92-96 \cr
#' @keywords datasets
NULL


#' @name archipelagos41
#' @title DAISIE datalist object including bird phylogenetic data and
#' physical data for 41 archipelagos
#' @docType data
#' @format A datalist containing data on the 41 archipelagos studied in Valente et
#' al 2020 (Main Dataset D1). Contains colonisation and branching times for bird species
#' in each of the archipelagos. It also contains information on archipelago name, area,
#' age and distance from the nearest mainland.
#' @description A datalist with 41 items representing the 41 archipelagos.
#' Each archipelago can be called separately using archipelagos41[[x]] with x being a
#' number between 1 and 41. Using archipelagos41[[x]][[1]] will show just the
#' top part of the archipelago item where the archipelago name and physical features
#' are displayed. The structure of each of the archipelagos is the same as regular DAISIE
#' datalist generated using DAISIE_dataprep.
#' @source Valente L, Phillimore AB, Melo M, Warren BH, Clegg SM, Havenstein K,
#'  Tiedemann R, Illera JC, Thébaud C, Aschenbach T, Etienne RS. A simple dynamic model
#'  explains island bird diversity worldwide (2020) Nature, 579, 92-96 \cr
#' @keywords datasets
NULL




#' @name frogs_datalist
#' @title Colonization and branching times of 5 Eleutherodactylus clades from Hispaniola island.
#' @docType data
#' @format A list with 6 elements, the first of which contains 2 elements and
#' the following 5 elements containing 5 components.
#' @description A list containing the colonization and branching times of the Eleutherodactylus frogs of
#' Hispaniola. Main dataset used in Etienne et al. This list can be generated
#' using the DAISIE_dataprep function, which converts a user-specified data table into a data object,
#' but the
#' object can of course also be entered directly. It is an R list object with
#' the following elements.\cr \cr
#' The first element of the list has two
#' components: \cr \cr
#' \code{$island_age} - the island age \cr
#' \code{$not_present} - the number of mainland lineages that are not present
#' on the island \cr \cr
#' The following elements of the list each contains
#' information on a single colonist lineage on the island and has 5
#' components:\cr \cr
#' \code{$colonist_name} - the name of the species or clade
#' that colonized the island \cr
#' \code{$branching_times} - island age followed by stem age of the population/species
#' in the case of Non-endemic, Non-endemic_MaxAge species and Endemic species with no close
#' relatives on the island. For endemic clades with more than one species on the island
#' (cladogenetic clades/ radiations) these should be island age followed by the
#' branching times of the island clade including the stem age of the clade.\cr
#' \code{$stac} - the status of the colonist \cr \cr
#' * Non_endemic_MaxAge: 1 \cr
#' * Endemic: 2 \cr
#' * Endemic&Non_Endemic: 3 \cr
#' * Non_endemic: 4 \cr
#' * Endemic_MaxAge: 5 or 6 \cr \cr
#' \code{$missing_species} - number of island species
#' that were not sampled for particular clade (only applicable for endemic
#' clades) \cr
#' \code{$type1or2} - whether the colonist belongs to type 1 or
#' type 2. In this dataset all are equal to 1. \cr
#' @seealso \code{\link{DAISIE_dataprep}}, \code{\link{DAISIE_ML}}, \code{\link{DAISIE_SR_ML}}
#' @source Etienne RS, Haegeman B, Dugo-Cota A, Vila C, Gonzalez-Voyer A & Valente L. The
#'  limits to ecological limits to diversification.\cr
#' @keywords datasets
NULL


#' @name frogs_datatable
#' @title Colonization and branching times of 5 Eleutherodactylus (frogs) clades from
#' the island of Hispaniola.
#' @docType data
#' @format A table with 5 rows and 4 columns.
#' @description A table containing the colonization and branching times of the Eleutherodacytlus
#'frogs of the island of Hispaniola (Greater Antilles). Each row on the table represents and independent
#' colonisation event. The table has four columns. \cr \cr
#' \code{$Clade_name} -
#' name of independent colonization event \cr
#' \code{$Status} - One of the
#' following categories: \cr
#' * Non_endemic: for non-endemic island species when an approximate time of colonisation is known \cr
#' * Non_endemic_MaxAge: for non-endemic island species when colonisation time is unknown \cr
#' * Endemic: for endemic species when an approximate colonisation time is known \cr
#' * "Endemic_MaxAge": applies to endemic species or endemic clades for cases where the
#' colonisation time is unknown, or when
#' the user wants to specify an upper bound for colonisation.
#' This could for example apply to endemic species that have recently gone extinct because
#' of anthropogenic causes, and which are not included
#' in the phylogeny ("NA" should be given in the branching times column). It
#' could also apply to insular radiations with long stem branches, for which the
#' time of the first cladogenetic event is known, but the precise time of colonisation
#' is not.\cr
#' * Endemic&Non_Endemic: when endemic clade and mainland ancestor has
#' re-colonized \cr
#' @source Etienne RS, Haegeman B, Dugo-Cota A, Vila C, Gonzalez-Voyer A & Valente L. The
#' limits to ecological limits to diversification.\cr
#' @keywords datasets
NULL

#' @name stac_table
#' @title Explanatory table on meaning and use of \code{stac} settings
#' @docType data
#' @format A table with 9 rows and 5 columns.
#' @description
#' A table containing the information regarding the meaning
#' of the \code{stac} codes utilised by DAISIE's ML functions. It is used to
#' render the "DAISIE \code{stac} values" vignette. \cr \cr
#' \code{stac} stands for "status of the clade" formed by the immigrant. It is
#' an important part of DAISIE objects that informs the likelihood functions
#' about the endemicity status and type of data available for each insular
#' clade. Each colonisation event that has extant species on the island needs to
#' have a \code{stac} value specified. This also has implications in what is
#' included in the DAISIE object \code{"branching_times"} vector, which also
#' described in this table.
#' The table is composed of the following columns:
#' \itemize{
#'   \item{$stac: A numeric with each stac code from 1 until 9.}
#'   \item{$Input for DAISIE_dataprep table "Clade_Name": A character
#'   with the possibilities for input in \link{DAISIE_dataprep}()
#'   of each clade, which are then translated to a numeric \code{stac} code as
#'   in \code{$stac}.}
#'   \item{$Type of species or clade": A character
#'   with the plain English explanation of the different possible kinds of
#'   island lineages that can be considered by DAISIE.}
#'   \item{$Input for DAISIE_dataprep table "Branching_times": A
#'   character with the plain English explanation of what the branching times
#'   vector in the DAISIE object should contain.}
#'   \item{$Colonisation times: A character with the plain English
#'   explanation of what the colonisation time in the branching times vector
#'   of the DAISIE object (first element of the vector) means.}
#' }
#' @keywords datasets
NULL
