#' Converts simulation output into island output
#'
#' @param stt_table a species-through-time table
#' @param totaltime simulated amount of time
#' @param island_spec matrix with species on island (state of system at each time point)
#' @param mainland_n number of mainland species
#' @param keep_final_state logical indicating if final state of simulation 
#' should be returned. Default is \code{FALSE}
#'
#' @return list with the island information, composed stt table, branching times of extant
#' species, status os species on the island and number of missing species.
DAISIE_create_island <- function(stt_table, 
                                 totaltime, 
                                 island_spec, 
                                 mainland_n,
                                 keep_final_state = FALSE) {
  
  ### if there are no species on the island branching_times = island_age, stac = 0, missing_species = 0 
  if (length(island_spec[,1]) == 0) {
    
    
    
    
    
    
    if (keep_final_state == TRUE) {
      island <- list(stt_table = stt_table,
                     branching_times = totaltime,
                     stac = 0,
                     missing_species = 0, island_spec = island_spec)
    } else {
      island <- list(stt_table = stt_table,
                     branching_times = totaltime,
                     stac = 0,
                     missing_species = 0)
    }
    
    
    
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
    island_spec[, "branching time (BP)"] <- totaltime - as.numeric(island_spec[, "branching time (BP)"])
    island_spec[, "Colonisation time (BP)"] <- totaltime - as.numeric(island_spec[, "Colonisation time (BP)"])
    
    if (mainland_n == 1) {
      island <- DAISIE_ONEcolonist(totaltime,
                                   island_spec,
                                   stt_table,
                                   keep_final_state = keep_final_state)
    } else if (mainland_n > 1) {
      
      ### number of colonists present
      colonists_present <- sort(as.numeric(unique(island_spec[, 'Mainland Ancestor'])))
      number_colonists_present <- length(colonists_present) 
      
      island_clades_info <- list()  
      for (i in 1:number_colonists_present) {
        subset_island <- island_spec[which(island_spec[, 'Mainland Ancestor'] == colonists_present[i]),] 
        
        if (class(subset_island) != 'matrix') {
          subset_island <- rbind(subset_island[1:7])
          colnames(subset_island) <- cnames
        }
        
        island_clades_info[[i]] <- DAISIE_ONEcolonist(totaltime,
                                                      island_spec = subset_island,
                                                      stt_table = NULL,
                                                      keep_final_state = keep_final_state)
        island_clades_info[[i]]$stt_table <- NULL
      }
      if (keep_final_state == FALSE) {
        island <- list(stt_table = stt_table,
                       taxon_list = island_clades_info)
        
      } else {
        island <- list(stt_table = stt_table,
                       taxon_list = island_clades_info, island_spec = island_spec)
      }
    }
  }
  
  return(island)
}