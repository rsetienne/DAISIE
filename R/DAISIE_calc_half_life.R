#' Calculates the exact half-life from a DAISIE simulation that has reached
#' growth or decayed to half way between the initial conditions and the 
#' equilibrium
#'
#' @param island_replicates output from DAISIE_sim_core
#' @param pars the parameters used for the simulation
#'
#' @return a list of half-lives
#' @export
DAISIE_calc_half_life <- function(island_replicates, pars, divdepmodel) {
  if (divdepmodel == "IW") {
    N0 <- list()
    for (i in 1:length(island_replicates)) {
      N0[[i]] <- sum(island_replicates$stt_tables[[i]][1, 2:4])
    }
    spec_half <- list()
    for (i in 1:length(island_replicates)) {
      spec_half[[i]] <- N0[[i]] - ((N0[[i]] - pars[3]) / 2)
      spec_half[[i]] <- round(spec_half[[i]], digits = 0)
    }
    total_stt <- list()
    for (i in 1:length(island_replicates)) {
      total_stt[[i]] <- apply(island_replicates[[i]]$stt_table[, 2:4], 1, sum)
    }
    for (i in 1:length(island_replicates)) {
      if (all(total_stt[[i]] == 0)) {
        stop("The island has zero species through time and therefore has no half-
           life.")
      }
    }
    row_spec_half <- list()
    for (i in 1:length(island_replicates)) {
      if (length(which(total_stt[[i]] == spec_half[[i]])) == 0) {
        print(paste("Half-life is yet to be reached use for replicate", i, sep = ""))
        row_spec_half[[i]] <- NA
      } else {
      row_spec_half[[i]] <- min(which(total_stt[[i]] == spec_half[[i]]))
      }
    }
    half_life <- list()
    for (i in 1:length(island_replicates)) {
      if (is.na(row_spec_half[[i]])) {
        half_life[[i]] <- NA
      } else {
        half_life[[i]] <- as.numeric(island_replicates[[i]]$stt_table[1, 1]) -
          (island_replicates[[i]]$stt_table[[row_spec_half[[i]], 1]])  
      }
    }
  }
  if (divdepmodel == "CS") {
    island_replicates <- DAISIE_format_CS(island_replicates = island_replicates, time = time, M = M, sample_freq = 100000, island_type = island_type)
    N0 <- list()
    for (i in 1:length(island_replicates)) {
      N0[[i]] <- sum(island_replicates[[i]][[1]]$stt_all[1, 2:4])
    }
    spec_half <- list()
    for (i in 1:length(island_replicates)) {
      spec_half[[i]] <- N0[[i]] - ((N0[[i]] - pars[3]) / 2)
      spec_half[[i]] <- round(spec_half[[i]], digits = 0)
    }
    total_stt <- list()
    for (i in 1:length(island_replicates)) {
      total_stt[[i]] <- apply(island_replicates[[i]][[1]]$stt_all[, 2:4], 1, sum)
    }
    for (i in 1:length(island_replicates)) {
      if (all(total_stt[[i]] == 0)) {
        stop("The island has zero species through time and therefore has no half-
           life.")
      }
    }
    row_spec_half <- list()
    for (i in 1:length(island_replicates)) {
      if (length(which(total_stt[[i]] == spec_half[[i]])) == 0) {
        print(paste("Half-life is yet to be reached use for replicate", i, sep = ""))
        row_spec_half[[i]] <- NA
      }else {
        row_spec_half[[i]] <- min(which(total_stt[[i]] == spec_half[[i]]))
      }
    }
    half_life <- list()
    for (i in 1:length(island_replicates)) {
      if (is.na(row_spec_half[[i]])) {
        half_life[[i]] <- NA
      } else {
        half_life[[i]] <- as.numeric(island_replicates[[i]][[1]]$stt_all[1, 1]) -
          (island_replicates[[i]][[1]]$stt_all[[row_spec_half[[i]], 1]])
    }
    }
  }
  return(half_life)
}