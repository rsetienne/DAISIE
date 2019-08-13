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

#' Calculates the half-life from a DAISIE simulation using the exponential 
#' model of Diamond (1972)
#' 
#' @param island_replicates output from DAISIE_sim_core
#' @param time time for simulation
#' @param pars the parameters used for the simulation
#'
#' @return a list of half-lives
#' @export
DAISIE_calc_half_life_exp_model <- function(island_replicates, time, pars) {
  if (divdepmodel == "IW") {
    N0 <- list()
    for (i in 1:length(island_replicates)) {
      N0[[i]] <- sum(island_replicates$stt_tables[[i]][1, 2:4])
    }
    last_row <- list()
    species_at_present <- list()
    half_life <- list()
    K <- pars[3]
    for (i in 1:length(island_replicates)) {
      last_row[[i]] <- nrow(island_replicates[[i]]$stt_table)
      species_at_present[[i]] <- sum(island_replicates[[i]]$stt_table[last_row[[i]], 2:4])
      half_life[[i]] <- time / (-log((species_at_present[[i]] - K) / 
                                       (N0[[i]] - K)))
    }
  }
  if (divdepmodel == "CS") {
    stt_list <- list()
    for (rep in 1:length(island_replicates)) {
      stt_list[[rep]] <- list()
      for (spec in 1:length(island_replicates[[1]])) {
        stt_list[[rep]][[spec]] <- island_replicates[[rep]][[spec]]$stt_table
      }
    }
    stt_all <- lapply(stt_list, function(x){
      stt_all <- lapply(x, as.data.frame)
      stt_all <- do.call(rbind, stt_all)
      stt_all <- split(stt_all, stt_all$Time)
      stt_all <- lapply(stt_all, function(z){
        (apply(z[, 2:4], 2, sum))
      })
    })
    stt_all <- lapply(stt_all, function(z) {
      apply(matrix(ncol = 4, nrow = length(z), c(names(z), do.call(rbind, z))),
            2, as.numeric)
    })
    stt_all <- lapply(stt_all, function(z) {
      apply(z, 2, rev)
    })
    for (i in 1:length(stt_all)) {
      colnames(stt_all[[i]]) <- c("Time", "nI", "nA", "nC")
    }
    N0 <- list()
    for (i in 1:length(island_replicates)) {
      N0[[i]] <- sum(island_replicates$stt_tables[[i]][1, 2:4])
    }
    last_row <- list()
    species_at_present <- list()
    half_life <- list()
    K <- pars[3]
    for (i in 1:length(island_replicates)) {
      last_row[[i]] <- nrow(island_replicates[[i]]$stt_table)
      species_at_present[[i]] <- sum(island_replicates[[i]]$stt_table[last_row[[i]], 2:4])
      half_life[[i]] <- time / (-log((species_at_present[[i]] - K) /
                                       (N0[[i]] - K)))
    }
  }
  return(half_life)
}

#' Estimates half-life based on a DAISIE simulation
#'
#' @param simulation DAISIE_sim output
#' @param pars Numeric vector with simulation parameters
#'
#' @author Joshua Lambert
#'
#' @return List of numerics with half-lives
#' @export
#'
DAISIE_estimate_half_life <- function(simulation, pars) {
  N0 <- list()
  for (i in 1:length(simulation)) {
    N0[[i]] <- sum(simulation[[i]][[1]]$stt_tables[1, 2:4])
  }
  spec_half <- list()
  for (i in 1:length(simulation)) {
    spec_half[[i]] <- N0[[i]] - ((N0[[i]] - pars[3]) / 2)
    spec_half[[i]] <- round(spec_half[[i]], digits = 0)
  }
  total_stt <- list()
  for (i in 1:length(simulation)) {
    total_stt[[i]] <- apply(simulation[[i]][[1]]$stt_all[, 2:4], 1, sum)
  }
  time_seq <- simulation[[1]][[1]]$stt_all[, 1]
  splines <- list()
  for (i in 1:length(simulation)) {
    splines[[i]] <- smooth.spline(x = rev(time_seq),
                                  y = total_stt[[i]],
                                  df = length(time_seq))
  }
  half_life_estimate <- list()
  for (i in 1:length(simulation)) {
    half_life_estimate[[i]] <- predict(splines[[i]], spec_half[[i]])
    half_life_estimate[[i]] <- half_life_estimate[[i]]$y
  }
  return(half_life_estimate)
}