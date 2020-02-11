#' Does somethings
#'
#' @inheritParams default_params_doc
#'
#' @return something
#' @export
DAISIE_make_archipelago <- function(archipelago,
                                    M,
                                    phylo_data,
                                    archipelago_data) {
  the_row <- which(archipelago_data[, "Archipelago"] == archipelago)
  island_age <- archipelago_data[the_row, "Age_oldest.Myr."]
  area <- archipelago_data[the_row, "Area.km2."]
  distance_continent <- archipelago_data[the_row, "Distance_continent"]
  distance_nearest_big <- archipelago_data[the_row, "Distance_nearest_big"]

  if (archipelago_data[the_row, "Total_species"] == 0) {
    archipelago_daisie <- list()
    archipelago_daisie[[1]] <- list(island_age = island_age, not_present = M)
    } else {
  archi_subset <- phylo_data[which(phylo_data[, "Archipelago"] == archipelago), ]
  da <- paste(archi_subset[, "Genus"], archi_subset[, "Species"], archi_subset[, "Subspecies"], sep = " ")
  db <- as.character(archi_subset[, "DAISIE_STATUS"])
  dc <- as.numeric(archi_subset[, "Missing_spec"])
  de <- as.character(archi_subset[, "BRANCHING_TIMES"])
  de[is.na(de)] <- island_age
  archi_daisie <- data.frame(da, db, dc, de)
  colnames(archi_daisie) <- c("Clade_name", "Status", "Missing_species", "Branching_times")
  archi_daisie <- archi_daisie[rev(order(archi_daisie[, "Branching_times"])), ]
  archipelago_daisie <- DAISIE_dataprep(archi_daisie, island_age, M = M)
    }
  archipelago_daisie[[1]]$area <- area
  archipelago_daisie[[1]]$distance_continent <- distance_continent
  archipelago_daisie[[1]]$distance_nearest_big <- distance_nearest_big
  archipelago_daisie[[1]]$name <- archipelago

  archipelago_daisie <- Add_brt_table(archipelago_daisie)
  archipelago_daisie[[1]]$brts_table <- NULL

  return(archipelago_daisie)
}

#' Does something
#'
#' @param archipelago_list  something
#' @param M  something
#' @param phylo_data  something
#' @param archipelago_data  something
#'
#' @return  something
#' @export
DAISIE_make_global <- function(archipelago_list, M, phylo_data, archipelago_data) {
  global_object <- list()
  for (i in 1:length(archipelago_list)) {
    the_archipelago <- as.character(archipelago_list[i])
    print(the_archipelago)
    global_object[[i]] <- DAISIE_make_archipelago(the_archipelago, M, phylo_data, archipelago_data)
  }
  return(global_object)
}
