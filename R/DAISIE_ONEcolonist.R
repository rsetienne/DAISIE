#' Title
#'
#' @inheritParams default_params_doc
#'
#' @return a list with these elements:
#' \itemize{
#'   item{[1]: stt_table, the same stt_table as put in}
#'   item{[2]: branching_times, branching times}
#'   item{[3]: stac, ?statuses}
#'   item{[4]: missing_species, ?the number of missing species}
#'   item{[5]: other_clades_same_ancestor, ?no idea}
#'   item{[6]: non-endemic species}
#'   item{[7]: endemic species}
#'   }
DAISIE_ONEcolonist <- function(time,
                              island_spec,
                              stt_table) {
  ### number of independent colonisations
  uniquecolonisation <- as.numeric(unique(
    island_spec[, "Colonisation time (BP)"]))
  number_colonisations <- length(uniquecolonisation)
  ### if there is only one independent colonisation - anagenetic and
  ### cladogenetic species are classed as stac=2; immigrant classed as stac=4:
  if (number_colonisations == 1) {
    if (island_spec[1, "Species type"] == "I") {
      descendants <- list(stt_table = stt_table,
                         branching_times = c(time, as.numeric(island_spec[1, "Colonisation time (BP)"])),
                         stac = 4,
                         missing_species = 0)
    }
    if (island_spec[1, "Species type"] == "A") {
      descendants <- list(stt_table = stt_table,
                         branching_times = c(time, as.numeric(island_spec[1, "Colonisation time (BP)"])),
                         stac = 2,
                         missing_species = 0)
    }
    if (island_spec[1, "Species type"] == "C") {
      descendants <- list(stt_table = stt_table,
                         branching_times = c(time, rev(sort(as.numeric(island_spec[, "branching time (BP)"])))),
                         stac = 2,
                         missing_species = 0)
    }
  }

  ### if there are two or more independent colonisations, all species are classed as stac=3 and put within same list item:
  if (number_colonisations > 1) {
    descendants <- list(stt_table = stt_table,
                       branching_times = NA,
                       stac = 3,
                       missing_species = 0,
                       other_clades_same_ancestor = list())

    btimes_all_clado_desc <- rev(sort(as.numeric(island_spec[, "branching time (BP)"])))

    if (length(btimes_all_clado_desc) != 0) { descendants$branching_times= c(time, btimes_all_clado_desc)}
    if (length(btimes_all_clado_desc)==0) { descendants$branching_times= c(time, max(as.numeric(island_spec[,"Colonisation time (BP)"])))}

    ### create table with information on other clades with same ancestor, but this information is not used in DAISIE_ML
    oldest <- which(as.numeric(island_spec[,"Colonisation time (BP)"]) == max(as.numeric(island_spec[,"Colonisation time (BP)"])))

    youngest_table = island_spec[-oldest,]
    if (is.character(youngest_table) && !is.matrix(youngest_table)) { ####### TODO: WILL FAIL IN DEVELOP ON rsetienne/DAISIE
      youngest_table = t(as.matrix(youngest_table))
    }

    uniquecol <- as.numeric(unique(youngest_table[,"Colonisation time (BP)"]))

    for (colonisation in 1:length(uniquecol)) {
      descendants$other_clades_same_ancestor[[colonisation]] = list(brts_miss = NA,species_type = NA)

      samecolonisation = which(as.numeric(youngest_table[,"Colonisation time (BP)"]) == uniquecol[colonisation])

      if (youngest_table[samecolonisation[1],"Species type"] == "I") {
        descendants$other_clades_same_ancestor[[colonisation]]$brts_miss = as.numeric(youngest_table[samecolonisation,"Colonisation time (BP)"])
        descendants$other_clades_same_ancestor[[colonisation]]$species_type = "I"
      }

      if (youngest_table[samecolonisation[1],"Species type"] == "A") {
        descendants$other_clades_same_ancestor[[colonisation]]$brts_miss =  as.numeric(youngest_table[samecolonisation,"Colonisation time (BP)"])
        descendants$other_clades_same_ancestor[[colonisation]]$species_type = "A"
      }

      if (youngest_table[samecolonisation[1],"Species type"] == "C") {
        descendants$other_clades_same_ancestor[[colonisation]]$brts_miss = rev(sort(as.numeric(youngest_table[samecolonisation,"branching time (BP)"])))
        descendants$other_clades_same_ancestor[[colonisation]]$species_type = "C"
      }
    }
  }
  return(descendants)
}
