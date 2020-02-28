#' Formats guild-wide simulation output into standard DAISIE
#' list output
#'
#' @inheritParams default_params_doc
#'
#' @return List with GW DAISIE simulation output
#' @export
#'
DAISIE_format_GW <- function(island_replicates,
                            time,
                            M,
                            sample_freq,
                            num_guilds,
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
    island_list <- list()
    stt_list <- list()
    for (i in 1:num_guilds) {
      stt_list[[i]] <- full_list[[i]]$stt_table
    }
    stt_all <- matrix(ncol = 5, nrow = sample_freq + 1)
    colnames(stt_all) <- c("Time", "nI", "nA", "nC", "present")
    stt_all[, "Time"] <- rev(seq(from = 0,
                               to = totaltime,
                               length.out = sample_freq + 1))
    immig_spec <- c()
    ana_spec <- c()
    for (i in 1:num_guilds) {
      immig_spec[[i]] <- sum(full_list[[i]]$stt_table[1, 2])
      ana_spec[[i]] <- sum(full_list[[i]]$stt_table[1, 3])
      }
    immig_spec <- sum(immig_spec)
    ana_spec <- sum(ana_spec)
    stt_all[1, 2:5] <- c(immig_spec, ana_spec, 0, 0)
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
      }
      if (number_present == 0) {
        island_list <- list()
        island_list[[1]] <- list(island_age = totaltime,
                                 not_present = M, #check what not_present should be
                                 stt_all = stt_all)
        island_list[[2]] <- list(branching_times = totaltime,
                                 stac = 0,
                                 missing_species = 0)
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
