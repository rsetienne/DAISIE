#' Converts simulation output into island output
#'
#' @inheritParams default_params_doc
#'
#' @return list with the island information, composed stt table,
#' branching times of extant species, status of species on
#' the island and number of missing species.
DAISIE_create_island <- function(stt_table,
                                 totaltime,
                                 island_spec,
                                 mainland_n,
                                 trait_pars = NULL) {
  
  if (!is.null(trait_pars)) {
    return(
      DAISIE_create_island_trait(
        stt_table = stt_table,
        totaltime = totaltime,
        island_spec = island_spec,
        mainland_n = mainland_n,
        trait_pars = trait_pars
      )
    )
  }
  ### if there are no species on the island branching_times = island_age,
  ### stac = 0, missing_species = 0
  if (length(island_spec[, 1]) == 0) {
    island <- list(stt_table = stt_table,
                   branching_times = totaltime,
                   stac = 0,
                   missing_species = 0)
  } else {
    cnames <- c("Species",
                "Mainland Ancestor",
                "Colonisation time (BP)",
                "Species type",
                "branch_code",
                "branching time (BP)",
                "Anagenetic_origin")
    colnames(island_spec) <- cnames
    ### set ages as counting backwards from present
    island_spec[, "branching time (BP)"] <- totaltime -
      as.numeric(island_spec[, "branching time (BP)"])
    island_spec[, "Colonisation time (BP)"] <- totaltime -
      as.numeric(island_spec[, "Colonisation time (BP)"])
    if (mainland_n == 1) {
      island <- DAISIE_ONEcolonist(totaltime,
                                   island_spec,
                                   stt_table)
    } else if (mainland_n > 1) {
      ### number of colonists present
      colonists_present <- sort(as.numeric(unique(
        island_spec[, "Mainland Ancestor"])))
      number_colonists_present <- length(colonists_present)
      island_clades_info <- list()
      for (i in 1:number_colonists_present) {
        subset_island <- island_spec[which(island_spec[, "Mainland Ancestor"] ==
                                             colonists_present[i]), ]
        if (!is.matrix(subset_island)) {
          subset_island <- rbind(subset_island[1:7])
          colnames(subset_island) <- cnames
        }
        island_clades_info[[i]] <- DAISIE_ONEcolonist(
          totaltime,
          island_spec = subset_island,
          stt_table = NULL)
        island_clades_info[[i]]$stt_table <- NULL
      }
        island <- list(stt_table = stt_table,
                       taxon_list = island_clades_info)
    }
  }
  return(island)
}

DAISIE_create_island_trait <- function(stt_table,
                                       totaltime,
                                       island_spec,
                                       mainland_n,
                                       trait_pars){
  
  ### if there are no species on the island branching_times = island_age, stac = 0, missing_species = 0
  if (length(island_spec[,1]) == 0) {
    island <- list(stt_table = stt_table,
                   branching_times = totaltime,
                   stac = 0,
                   missing_species = 0)
    
  } else {
    cnames <- c("Species",
                "Mainland Ancestor",
                "Colonisation time (BP)",
                "Species type",
                "branch_code",
                "branching time (BP)",
                "Anagenetic_origin",
                "trait_state")
    
    colnames(island_spec) <- cnames
    
    ### set ages as counting backwards from present
    island_spec[, "branching time (BP)"] <- totaltime - as.numeric(island_spec[, "branching time (BP)"])
    island_spec[, "Colonisation time (BP)"] <- totaltime - as.numeric(island_spec[, "Colonisation time (BP)"])
    
    mainland_ntotal = mainland_n + trait_pars$M2
    
    if (mainland_ntotal == 1) {
      island <- DAISIE_ONEcolonist(totaltime,
                                   island_spec,
                                   stt_table)
      
      
    } else if (mainland_ntotal > 1) {
      
      ### number of colonists present
      colonists_present <- sort(as.numeric(unique(island_spec[, 'Mainland Ancestor'])))
      number_colonists_present <- length(colonists_present) 
      
      island_clades_info <- list()  
      for (i in 1:number_colonists_present) {
        subset_island <- island_spec[which(island_spec[, "Mainland Ancestor"] ==
                                             colonists_present[i]), ]
        if (!is.matrix(subset_island)) {
          subset_island <- rbind(subset_island[1:8])
          colnames(subset_island) <- cnames
        }
        island_clades_info[[i]] <- DAISIE_ONEcolonist(
          totaltime,
          island_spec = subset_island,
          stt_table = NULL)
        island_clades_info[[i]]$stt_table <- NULL
      }
      island <- list(stt_table = stt_table,
                     taxon_list = island_clades_info)
    }
  }
  return(island)
}

