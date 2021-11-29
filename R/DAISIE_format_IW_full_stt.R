#' Formats clade-specific simulation output into standard
#' DAISIE list output with complete STT table
#'
#' @inheritParams default_params_doc
#'
#' @return List with IW DAISIE simulation output
DAISIE_format_IW_full_stt <- function(island_replicates,
                                      total_time,
                                      M,
                                      verbose) {
  several_islands <- list()

  for (rep in seq_along(island_replicates)) {

    stt_all <- island_replicates[[rep]]$stt_table
    island_list <- list()

    if (sum(stt_all[nrow(stt_all), 2:4]) == 0) {
      island_list[[1]] <- list(
        island_age = total_time,
        not_present = M,
        stt_all = stt_all
      )
    } else {
      taxon_list_size <- length(island_replicates[[rep]]$taxon_list)
      island_list[[1]] <- list(
        island_age = total_time,
        not_present = M - taxon_list_size,
        stt_all = stt_all
      )
      if (taxon_list_size != 0) {
        for (y in seq_len(taxon_list_size)) {
          island_list[[y + 1]] <- island_replicates[[rep]]$taxon_list[[y]]
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
