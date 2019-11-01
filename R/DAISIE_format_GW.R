#' Formats guild-wide simulation output into standard DAISIE
#' list output
#'
#' @param island_replicates List output from DAISIE_sim_core.
#' @param time Numeric double with total time of simulation.
#' @param M Int stating number of mainland species.
#' @param sample_freq Int stating how often results are
#' sampled for plotting.
#' @param island_type type of island for simulation.
#' @param num_guilds number of guilds on the mainland.
#' @param start_midway Logical stating if simulation starts at t > 0.
#' @param verbose Logical controling if progress is printed to console.
#'
#' @return List with GW DAISIE simulation output
#' @export
#'
DAISIE_format_GW <- function(island_replicates,
                            time,
                            M,
                            sample_freq,
                            island_type,
                            num_guilds,
                            start_midway = FALSE,
                            verbose = TRUE) {
  totaltime <- time
  several_islands <- list()
  for (rep in 1:length(island_replicates)) {
    full_list <- island_replicates[[rep]]
    stac_list <- list()
    taxon_list_stac_vec <- c()
    for (i in 1:length(full_list)) {
      if (is.null(full_list[[i]]$taxon_list)) {
        stac_list[[i]] <- full_list[[i]]$stac
      } else {
        for (j in 1:length(full_list[[i]]$taxon_list)) {
         taxon_list_stac_vec[j] <- full_list[[i]]$taxon_list[[j]]$stac
         stac_list[[i]] <- taxon_list_stac_vec
        }
      }
    }
    stac_vec <- unlist(stac_list)
    number_not_present <- length(which(stac_vec == 0))
    present <- which(stac_vec != 0)
    number_present <- length(present)
    init_nonend_spec_list <- list()
    init_end_spec_list <- list()
    carrying_capacity_list <- list()
    for (i in 1:length(full_list)) {
      if (is.null(full_list[[i]]$taxon_list)) {
        init_nonend_spec_list[[i]] <- full_list[[i]]$init_nonend_spec
        init_end_spec_list[[i]] <- full_list[[i]]$init_end_spec
        carrying_capacity_list[[i]] <- full_list[[i]]$carrying_capacity
      } else {
        for (j in 1:length(full_list[[i]]$taxon_list)) {
          init_nonend_spec_list[[i]] <- full_list[[i]]$taxon_list[[1]]$init_nonend_spec
          init_end_spec_list[[i]] <- full_list[[i]]$taxon_list[[1]]$init_end_spec
          carrying_capacity_list[[i]] <- full_list[[i]]$taxon_list[[1]]$carrying_capacity
        }
      }
    }
    init_nonend_spec <- sum(unlist(init_nonend_spec_list))
    init_end_spec <- sum(unlist(init_end_spec_list))
    carrying_capacity_vec <- unlist(carrying_capacity_list)
    island_list <- list()
    stt_list <- list()
    for (i in 1:num_guilds) {
      stt_list[[i]] <- full_list[[i]]$stt_table
    }
    stt_all <- matrix(ncol = 5, nrow = sample_freq + 1)
    colnames(stt_all) <- c("Time", "nI", "nA", "nC", "present")
    stt_all[,"Time"] <- rev(seq(from = 0,
                               to = totaltime,
                               length.out = sample_freq + 1))
    if (start_midway == FALSE) { #where should this if statement stop?
      if (island_type == "oceanic") {
        stt_all[1, 2:5] <- c(0, 0, 0, 0)
      } else {
        stt_all[1, 2:5] <- c(init_nonend_spec, init_end_spec, 0, 0)
      }
    }
    if (start_midway == TRUE) {
        for (x in 1:M) {
          stt_all[1, 2:5] <- stt_list[[x]][max(which(stt_list[[x]][, "Time"] >=
                                                       totaltime)), 2:4]
        }
      }
      for (i in 2:nrow(stt_all)) {
        the_age <- stt_all[i, "Time"]
        store_richness_time_slice <- matrix(nrow = num_guilds, ncol = 3)
        colnames(store_richness_time_slice) <- c("I", "A", "C")
        for (x in 1:num_guilds) {
          row_index <- max(which(stt_list[[x]][, "Time"] >= the_age))
          store_richness_time_slice[x, ] <- stt_list[[x]][max(which(stt_list[[x]][, "Time"] >= the_age)), 2:4]
        }
        count_time_slice <- store_richness_time_slice[, 1] +
          store_richness_time_slice[, 2] +
          store_richness_time_slice[, 3]
        present_time_slice <- rep(0, num_guilds)
        present_time_slice[which(count_time_slice > 0)] <- 1
        store_richness_time_slice <- cbind(store_richness_time_slice,
                                           present_time_slice)
        stt_all[i, c(2:5)] <- apply(store_richness_time_slice, 2, sum)
      }
    island_list[[1]] <- list(island_age = totaltime,
                             not_present = number_not_present,
                             stt_all = stt_all)
      if (number_present > 0) {
          for (i in 1:num_guilds) {
            full_list[[i]]$stt_table <- NULL
              if (is.null(full_list[[i]]$taxon_list)) {
              island_list <- append(island_list, list(full_list[[i]]))
            } else {
              for (j in 1:length(full_list[[i]]$taxon_list)) {
                taxon_sublist_name <- paste0("taxon_sublist_", c(1:length(full_list[[i]]$taxon_list)))
                assign(taxon_sublist_name[j], full_list[[i]]$taxon_list[[j]])
                island_list <- append(island_list, list(get(taxon_sublist_name[[j]])))
                }
                }
              }
        island_list[[length(island_list) + 1]] <- list(
          init_nonend_spec = init_nonend_spec,
          init_end_spec = init_end_spec,
          all_carrying_capacities = carrying_capacity_vec)
      }
      if (number_present == 0) {
        island_list <- list()
        island_list[[1]] <- list(island_age = totaltime,
                                 not_present = M, #check what not_present should be
                                 stt_all = stt_all)
        island_list[[2]] <- list(branching_times = totaltime,
                                 stac = 0,
                                 missing_species = 0,
                                 init_nonend_spec = init_nonend_spec,
                                 init_end_spec = init_end_spec,
                                 carrying_capacity = "N/A",
                                 all_carrying_capacities = carrying_capacity_vec)
        }

      several_islands[[rep]] <- island_list
      if (verbose == TRUE) {
        print(paste("Island being formatted: ",
                    rep,
                    "/",
                    length(island_replicates),
                    sep = ""))
      }
    }
  return(several_islands)
}
