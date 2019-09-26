DAISIE_format_GW = function(island_replicates,
                            time,
                            M,
                            sample_freq,
                            island_type,
                            num_guild,
                            start_midway = FALSE,
                            verbose = TRUE) {
  totaltime <- time
  several_islands = list()
  for(rep in 1:length(island_replicates)) {
    full_list = island_replicates[[rep]]
    stac_vec = unlist(full_list)[which(names(unlist(full_list)) == "stac")]
    number_not_present = length(which(stac_vec == 0))
    present = which(stac_vec!=0)
    number_present = length(present)
    type_vec = unlist(full_list)[which(names(unlist(full_list)) == "type1or2")]
    prop_type2_pool = length(which(type_vec == 2)) / M
    number_type2_cols = length(which(match(which(stac_vec != 0),
                                           which(type_vec == 2)) > 0))
    number_type1_cols = number_present-number_type2_cols
    island_list = list()
    for(i in 1:(number_present + 1)) {
      island_list[[i]] = list()
    }
    ### all species
    stt_list = list()
    for(i in 1:(M / num_guild)) {
      stt_list[[i]] = full_list[[i]]$stt_table
    }
    stt_all = matrix(ncol = 5,nrow = sample_freq + 1)
    colnames(stt_all) = c("Time","nI","nA","nC","present")
    stt_all[,"Time"] = rev(seq(from = 0,
                               to = totaltime,
                               length.out = sample_freq + 1))
    if (start_midway == FALSE) {
      if (island_type == "oceanic") {
        stt_all[1, 2:5] <- c(0, 0, 0, 0)
      } else {
        immig_spec = c()
        ana_spec = c()
        for (i in 1:M) {
          immig_spec[[i]] <- sum(full_list[[i]]$stt_table[1,2])
          ana_spec[[i]] <- sum(full_list[[i]]$stt_table[1, 3])
        }
        immig_spec <- sum(immig_spec)
        ana_spec <- sum(ana_spec)
        stt_all[1, 2:5] = c(immig_spec, ana_spec, 0, 0)
      }
      if (start_midway == TRUE) {
        for (x in 1:M) {
          stt_all[1, 2:5] <- stt_list[[x]][max(which(stt_list[[x]][, "Time"] >=
                                                       totaltime)), 2:4]
        }
      }
      for (i in 2:nrow(stt_all)) {
        the_age <- stt_all[i, "Time"]
        store_richness_time_slice = matrix(nrow = M, ncol = 3)
        colnames(store_richness_time_slice) <- c("I", "A", "C")
        for (x in 1:(M / num_guild)) {
          testit::assert(x >= 1)
          testit::assert(x <= length(stt_list))
          testit::assert(is.matrix(stt_list[[x]]))
          testit::assert("Time" %in% colnames(stt_list[[x]]))
          testit::assert(!all(is.na(stt_list[[x]][, "Time"])))
          testit::assert(!all(is.infinite(stt_list[[x]][, "Time"])))
          testit::assert(!is.na(the_age))
          row_index <- max(which(stt_list[[x]][, "Time"] >= the_age))
          testit::assert(!is.na(row_index))
          testit::assert(row_index >= 1)
          testit::assert(row_index <= nrow(stt_list[[x]]))
          store_richness_time_slice[x, ] <- stt_list[[x]][max(which(stt_list[[x]][, "Time"] >= the_age)), 2:4]
        }
        count_time_slice <- store_richness_time_slice[, 1] +
          store_richness_time_slice[, 2] +
          store_richness_time_slice[,3]
        present_time_slice <- rep(0, M)
        present_time_slice[which(count_time_slice > 0)] <- 1
        store_richness_time_slice <- cbind(store_richness_time_slice,
                                           present_time_slice)
        stt_all[i, c(2:5)] <- apply(store_richness_time_slice, 2, sum)
      }

      if (number_type2_cols > 0) {
        stop("Two type simulation cannot be run with divdepmodel = 'GW'")
      } else {
        island_list[[1]] <- list(island_age = totaltime,
                                 not_present = number_not_present,
                                 stt_all = stt_all)
      }
      if (number_present > 0) {
        for (i in 1:number_present) {
          island_list[[1 + i]] <- full_list[[present[i]]]
          island_list[[1 + i]]$stt_table <- NULL
        }
      }
      if (number_present == 0) {
        island_list <- list()
        island_list[[1]] <- list(island_age = totaltime,
                                 not_present = M,
                                 stt_all = stt_all)
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
  }
  return(several_islands)
}
