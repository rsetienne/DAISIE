
TRAISIE_sim_CS <- function(total_time,
                           mainland,
                           trait_pars,
                           replicates,
                           sample_freq = 100,
                           cond = 0,
                           verbose = TRUE,
                           files_to_write = 0,
                           num_observed_states,
                           num_hidden_states) {
  island_replicates <- list()

  for (rep in 1:replicates) {
    island_replicates[[rep]] <- list()
    full_list <- list()
    if (cond == 0) {
      number_present <- -1
    } else {
      number_present <- 0
    }

    while (number_present < cond) {

      # counts per mainland group (e.g., c(M1, M2, M3, ...))
      counts <- unlist(mainland)
      G      <- length(counts)
      cum    <- cumsum(counts)
      n_tot  <- sum(counts)

      full_list <- vector("list", n_tot)

      for (m_spec in seq_len(n_tot)) {

        # which group does this mainland species belong to?
        g <- which(m_spec <= cum)[1]

        # one-hot root state
        root <- rep(0L, num_observed_states)

        n <- num_observed_states * num_hidden_states
        if (g %in% 1:n) {
          gg <- ceiling(g / num_hidden_states)
        }

        root[gg] <- 1L
        mainland_root <- rep(0L, n)
        mainland_root[g] <- 1L

        # run model with group-specific mainland vector
        full_list[[m_spec]] <- TRAISIE_sim_core(
          time = total_time,
          mainland = as.list(mainland_root),          # e.g., list(1,0,0,...)
          trait_pars = trait_pars,
          num_observed_states = num_observed_states,
          num_hidden_states = num_hidden_states
        )

        if (!is.null(full_list[[m_spec]])) {
          full_list[[m_spec]]$root_state <- root
        }
      }
      stac_vec <- unlist(full_list)[which(names(unlist(full_list)) == "stac")]
      present <- which(stac_vec != 0)
      number_present <- length(present)
    }
    island_replicates[[rep]] <- full_list
    if (verbose == TRUE) {
      message("Island replicate ", rep)
    }
  }

  if (files_to_write == 0) {
    island_replicates <- DAISIE_format_CS(island_replicates =
                                            island_replicates,
                                          time = total_time,
                                          M = mainland[[1]],
                                          sample_freq = sample_freq,
                                          verbose = verbose)
  }

  if (files_to_write > 0) {
    rm(island_replicates)
    for (filenum in 1:files_to_write) {
      chunks <- ceiling(seq_along(1:replicates) / files_to_write)
      start <- min(which(chunks == filenum))
      end <- max(which(chunks == filenum))
      load(paste0("DAISIE_sims", start, "-", end, ".Rdata"))
      island_replicates <- DAISIE_format_CS(island_replicates =
                                              island_replicates,
                                            time = total_time,
                                            M = mainland[[1]],
                                            sample_freq = sample_freq,
                                            verbose = verbose)
      save(start, end, island_replicates,
           file = paste0("DAISIE_sims_formatted", start, "-", end, ".Rdata"))
    }
  }

  return(island_replicates)
}
