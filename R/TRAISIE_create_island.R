#' @keywords internal
TRAISIE_create_island <- function(stt_table,
                                  total_time,
                                  island_spec,
                                  mainland,
                                  trait_pars,
                                  num_observed_states,
                                  num_hidden_states) {

  ### Check if island_spec is a matrix and convert it if not


  if (length((island_spec)) > 1) {
    island_spec <- island_spec[order(island_spec[, 5]), ]
    if (!is.matrix(island_spec)) {

      island_spec <- t(as.matrix(island_spec))
    }
  }

  ### if there are no species on the island branching_times = island_age, stac = 0, missing_species = 0
  if (length(island_spec[, 1]) == 0) {
    island <- list(stt_table = stt_table,
                   branching_times = total_time,
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
                "trait_state",
                "connection")

    colnames(island_spec) <- cnames

    ### Set ages as counting backwards from present
    island_spec[, "branching time (BP)"] <- total_time - as.numeric(island_spec[, "branching time (BP)"])
    island_spec[, "Colonisation time (BP)"] <- total_time - as.numeric(island_spec[, "Colonisation time (BP)"])

    mainland_total <- sum(unlist(mainland))



    ### Number of colonists present
    colonists_present <- sort(as.numeric(unique(island_spec[, "Mainland Ancestor"])))
    number_colonists_present <- length(colonists_present)

    island_clades_info <- list()

    for (i in seq_along(island_spec[, 1])) {
      mainland_spec <- island_spec[i, 2]

      all_spec <- island_spec[which(island_spec[, "Mainland Ancestor"] == mainland_spec), ]

      if (!is.matrix(all_spec)) {
        cnames <- names(all_spec)
        all_spec <- rbind(all_spec[cnames])
        colnames(all_spec) <- cnames
      }

      col_times <- unique(stats::na.omit(suppressWarnings(
        as.numeric(all_spec[, "Colonisation time (BP)"])
      )))


      if ((length(col_times) == 1) ||
          (length(col_times) > 1 && any(all_spec[, "Species type"] == "I", na.rm = TRUE))) {

        subset_island <- island_spec[which(island_spec[, "Mainland Ancestor"] == as.character(mainland_spec)), ]

      } else if (length(col_times) > 1 &&  any(all_spec[, "Species type"] != "I", na.rm = TRUE)) {

        subset_island <- all_spec[
          all_spec[, "Colonisation time (BP)"] == island_spec[i, ][["Colonisation time (BP)"]],
          , drop = FALSE
        ]
      }


      if (!is.matrix(subset_island)) {
        subset_island <- rbind(subset_island[1:9])
        colnames(subset_island) <- cnames
      }

      island_clades_info[[i]] <- DAISIE:::DAISIE_ONEcolonist(
        total_time,
        island_spec = subset_island,
        stt_table = stt_table)

      island_clades_info[[i]]$stt_table <- stt_table
    }


    # Extracting taxon_list and handling matching colonization times
    for (i in seq_along(island_clades_info)) {
      # Extract colonization times from island_spec (it is a vector)
      colonization_times <- as.numeric(island_spec[which(island_spec[, "Mainland Ancestor"] == colonists_present[i]), "Colonisation time (BP)"])

      # Loop through taxon_list to find matching colonization times
      matching_taxa_list <- list()

      for (j in seq_along(island_clades_info)) {
        # Check if any of the colonization times match the branching times
        if (any(colonization_times == island_clades_info[[j]]$branching_times[2])) {
          matching_taxa_list <- append(matching_taxa_list, j)
        }
      }

      # If matches are found, add them to the corresponding taxon
      if (length(matching_taxa_list) > 0) {
        for (match in matching_taxa_list) {
          # Prepare the subset of island_spec for the current match
          island_clades_info[[match]]$island_spec <- list(island_spec[which(island_spec[, "Colonisation time (BP)"] == island_clades_info[[match]]$branching_times[2]), ])[1]
          isla <- list(island_spec[which(island_spec[, "Colonisation time (BP)"] == island_clades_info[[match]]$branching_times[2]), ])[1]


          subset_island <- all_spec[
            all_spec[, "Colonisation time (BP)"] == island_spec[i, ][["Colonisation time (BP)"]],
            , drop = FALSE
          ]
          # Ensure isla[[1]] is a data frame before indexing
          if (!is.matrix(isla[[1]])) {
            isla[[1]] <- matrix(isla[[1]], nrow = 1)
            colnames(isla[[1]]) <- cnames

          }


          # Ensure the "Mainland Ancestor" value is numeric
          mainland_ancestor_value <- as.numeric(isla[[1]][, "Mainland Ancestor"][1])

          # Initialize root_state vector with zeros
          root_state <- rep(0, length(mainland))

          # Convert mainland list to numeric vector
          mainland_counts <- unlist(mainland)

          # Compute cumulative sums to identify trait group boundaries
          cumulative <- cumsum(mainland_counts)

          # Find which group the ancestor belongs to
          group_index <- which(mainland_ancestor_value <= cumulative)[1]

          # Set the corresponding position to 1 (one-hot encoding)
          root_state[group_index] <- 1

          island_clades_info[[match]]$root_state <- root_state

          if (ncol(isla[[1]]) >= 9) {
            island_clades_info[[match]]$traits <- as.numeric(isla[[1]][, 8])
          } else {
            island_clades_info[[match]]$traits <- NA
            warning("Eight column not found in isla[[1]], assigned NA")
          }

          isla <- isla[[1]]
          # Always assign sampling_fraction
          island_clades_info[[match]]$sampling_fraction <- rep(1, num_observed_states)

          # Convert colonisation time column to numeric
          colonisation_times <- as.numeric(isla[, "Colonisation time (BP)"])

          # Check if there are at least 2 unique colonisation times
          if (length(unique(colonisation_times)) > 1) {
            # Find the smallest colonisation time
            min_col_time <- min(colonisation_times, na.rm = TRUE)

            # Keep only the rows where colonisation time is NOT the smallest
            isla  <-  isla[colonisation_times != min_col_time, , drop = FALSE]
          }

          if (length(isla[, 9]) > 1) {



            phy <- build_phylo_tree_from_island_spec(island_spec = isla)

            island_clades_info[[match]]$phylogeny <- phy
          } else {
            island_clades_info[[match]]$phylogeny <- NA
          }

        }

      }



    }
    hashes <- vapply(island_clades_info, digest::digest, character(1), algo = "xxhash64")

    keep   <- !duplicated(hashes)

    island_clades_info <- island_clades_info[keep]

    island <- island_clades_info[[1]]
  }

  return(island)
}
