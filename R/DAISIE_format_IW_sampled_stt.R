#' Formats island-wide simulation output into standard
#' DAISIE list output
#'
#' @inheritParams default_params_doc
#'
#' @return List with IW DAISIE simulation output
#' @keywords internal
DAISIE_format_IW_sampled_stt <- function(island_replicates,
                                         total_time,
                                         M,
                                         sample_freq,
                                         verbose) {

  several_islands <- list()
  for (rep in 1:length(island_replicates)) {
    the_island <- island_replicates[[rep]]
    stt_all <- matrix(ncol = 4, nrow = sample_freq + 1)
    colnames(stt_all) <- c("Time", "nI", "nA", "nC")
    stt_all[, "Time"] <- rev(seq(from = 0,
                                 to = total_time,
                                 length.out = sample_freq + 1))
    immig_spec <- the_island$stt_table[1, 2]
    ana_spec <- the_island$stt_table[1, 3]
    stt_all[1, 2:4] <- c(immig_spec,
                         ana_spec,
                         0)
    the_stt <- the_island$stt_table
    for (i in 2:nrow(stt_all)) {
      the_age <- stt_all[i, "Time"]
      stt_all[i, 2:4] <- the_stt[max(which(the_stt[, "Time"] >= the_age)), 2:4]
    }
    island_list <- list()
    if (sum(the_stt[nrow(the_stt), 2:4]) == 0) {
      island_list[[1]] <- list(
        island_age = total_time,
        not_present = M,
        stt_all = stt_all
      )
    } else {
      taxon_list_size <- length(the_island$taxon_list)
      island_list[[1]] <- list(
        island_age = total_time,
        not_present = M - taxon_list_size,
        stt_all = stt_all
      )
      if (taxon_list_size != 0) {
        for (y in seq_len(taxon_list_size)) {
          island_list[[y + 1]] <- the_island$taxon_list[[y]]
        }
      }
    }

    island_list <- add_brt_table(island_list)
    several_islands[[rep]] <- island_list

    if (verbose == TRUE) {
      message(
        "Island being formatted: ", rep, "/", length(island_replicates)
      )
    }
  }
  return(several_islands)
}
