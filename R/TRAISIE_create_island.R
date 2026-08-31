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



            phy <- TRAISIE_build_phylo_tree_from_island_spec(island_spec = isla)

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

#' @keywords internal
TRAISIE_build_phylo_tree_from_island_spec <- function(island_spec) {


  # Ensure required package is loaded
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("The 'ape' package is required but not installed. Please install it with install.packages('ape').")
  }
  island_spec <- island_spec[order(island_spec[, 5]), ]
  # Convert colonisation time column to numeric
  colonisation_times <- as.numeric(island_spec [, "Colonisation time (BP)"])

  # Check if there are at least 2 unique colonisation times
  if (length(unique(colonisation_times)) > 1) {
    # Find the smallest colonisation time
    min_col_time <- min(colonisation_times, na.rm = TRUE)

    # Keep only the rows where colonisation time is NOT the smallest
    island_spec  <-  island_spec [colonisation_times != min_col_time, , drop = FALSE]
  }



  # Inputs
  history_vec  <- island_spec[, 9]  # e.g., character vector like c("NA 1.7 1.3", ...)
  bt_vec       <- island_spec[, 6]
  col_time     <- island_spec[, 3] # e.g., numeric or character vector
  branch_codes <- island_spec[, 5]
  traits_vec   <- island_spec[, 8]


  # Convert bt_vec to numeric for comparison
  bt_vec_num <- suppressWarnings(as.numeric(bt_vec))
  tolerance  <- 1e-8

  # Initialize cleaned vector
  cleaned_history <- character(length(history_vec))

  # Loop through each element of history_vec
  for (i in seq_along(history_vec)) {

    # Split the string into elements
    vals <- unlist(strsplit(history_vec[i], " "))

    # Initialize list to keep valid values
    valid_vals <- character(0)

    for (val in vals) {
      # Keep "NA" as-is
      if (val == "NA") {
        valid_vals <- c(valid_vals, "NA")
      } else {
        num_val <- suppressWarnings(as.numeric(val))
        # Keep if close to any value in bt_vec_num
        if (!is.na(num_val)) {
          if (any(abs(num_val - bt_vec_num) < tolerance, na.rm = TRUE)) {
            valid_vals <- c(valid_vals, val)
          }
        }
      }
    }

    # Join the valid values back into a string
    cleaned_history[i] <- paste(valid_vals, collapse = " ")
  }
  # Step 2: Create value → label map
  all_vals <- unique(unlist(strsplit(cleaned_history, " ")))
  all_vals <- all_vals[all_vals != "NA" & all_vals != ""]
  label_map <- stats::setNames(seq_along(all_vals), all_vals)

  # Step 3: Recode each row with unique terminal node
  used_numbers <- c()
  next_label <- max(label_map) + 1

  # Step 3: Recode each row with unique terminal node
  translated_result <- vector("character", length(cleaned_history))
  used_numbers <- c()
  next_label <- max(label_map) + 1

  for (i in seq_along(cleaned_history)) {
    row_str <- cleaned_history[i]
    vals <- unlist(strsplit(row_str, " "))
    vals <- vals[vals != "NA" & vals != ""]

    labels <- as.character(label_map[vals])
    current_numbers <- as.numeric(labels)

    while (next_label %in% current_numbers || next_label %in% used_numbers) {
      next_label <- next_label + 1
    }

    final_numbers <- c(current_numbers, next_label)
    used_numbers <- c(used_numbers, next_label)
    next_label <- next_label + 1

    translated_result[i] <- paste(final_numbers, collapse = " ")
  }

  # Final output
  translated_history <- translated_result

  # Step 4: Create node vectors
  vecs <- lapply(translated_history, function(str_row) {
    as.integer(unlist(strsplit(str_row, " ")))
  })

  # Step 5: Extract edges
  edges <- do.call(rbind, lapply(vecs, function(path) {
    if (length(path) < 2) return(NULL)
    cbind(parent = path[-length(path)], child = path[-1])
  }))

  # Step 6: Identify and remove duplicate edges
  duplicated_rows <- duplicated(edges)
  duplicate_indices <- which(duplicated_rows)
  edges <- edges[!duplicated_rows, ]

  # Step 7: Compute edge lengths
  time_vecs_with_zero <- lapply(cleaned_history, function(row_str) {
    vals <- unlist(strsplit(as.character(row_str), " "))
    vals <- vals[vals != "NA" & vals != ""]
    numeric_vals <- as.numeric(vals)
    c(numeric_vals, 0)
  })

  edge_lengths <- unlist(lapply(time_vecs_with_zero, function(times) {
    if (length(times) < 2) return(numeric(0))
    utils::head(times, -1) - utils::tail(times, -1)
  }))

  if (length(duplicate_indices) > 0) {
    edge_lengths <- edge_lengths[-duplicate_indices]
  }

  # Step 8: Remap edges to comply with ape format
  all_nodes <- unique(as.vector(edges))
  parent_nodes <- unique(edges[, 1])
  child_nodes  <- unique(edges[, 2])
  tip_nodes <- sort(setdiff(child_nodes, parent_nodes))
  ntip <- length(tip_nodes)
  internal_nodes <- sort(setdiff(all_nodes, tip_nodes))

  node_order <- c(tip_nodes, internal_nodes)
  node_map <- stats::setNames(seq_along(node_order), node_order)

  edge_remapped <- matrix(
    c(node_map[as.character(edges[, 1])], node_map[as.character(edges[, 2])]),
    ncol = 2
  )
  # Step 9: Build phylo object
  tree <- list()
  tree$edge <- matrix(as.integer(edge_remapped), ncol = 2)
  tree$edge.length <- as.numeric(edge_lengths)
  tree$Nnode <- length(internal_nodes)
  tree$tip.label <- branch_codes
  class(tree) <- "phylo"

  return(tree)
}


