#' Formats clade-specific simulation output into standard
#' DAISIE list output
#'
#' @inheritParams default_params_doc
#'
#' @return List with CS DAISIE simulation output
DAISIE_format_CS_full_stt <- function(island_replicates,
                                      time,
                                      M,
                                      verbose = TRUE,
                                      trait_pars = NULL
) {
  totaltime <- time
  several_islands <- list()
  for (rep in seq_along(island_replicates)) {
    full_list <- island_replicates[[rep]]
    stac_vec <- unlist(full_list)[which(names(unlist(full_list)) == "stac")]
    number_not_present <- length(which(stac_vec == 0))
    present <- which(stac_vec != 0)
    number_present <- length(present)
    type_vec <- unlist(full_list)[which(names(unlist(full_list)) == "type1or2")]
    prop_type2_pool <- length(which(type_vec == 2)) / M
    number_type2_cols <- length(which(match(which(stac_vec != 0),
                                            which(type_vec == 2)) > 0))
    number_type1_cols <- number_present - number_type2_cols
    island_list <- list()
    for (i in 1:(number_present + 1)) {
      island_list[[i]] <- list()
    }
    ### all species
    stt_list <- list()
    if(is.null(trait_pars)){
      for (i in 1:M) {
        stt_list[[i]] <- full_list[[i]]$stt_table
      } 
    }else{
      for (i in 1:(M + trait_pars$M2)) {
        stt_list[[i]] <- full_list[[i]]$stt_table
      } 
    }
    
    #### Keep full STT ####
    stt_all <- create_full_CS_stt(
      stt_list = stt_list,
      stac_vec = stac_vec,
      totaltime = totaltime,
      trait_pars = trait_pars
    )
    ####  two trait states####
    if(!is.null(trait_pars)){
      immig_spec <- c()
      ana_spec <- c()
      immig_spec2 <- c()
      ana_spec2 <- c()
      for (i in 1:(M + trait_pars$M2)) {
        immig_spec[i] <- sum(full_list[[i]]$stt_table[1, 2])
        ana_spec[i] <- sum(full_list[[i]]$stt_table[1, 3])
        immig_spec2[i] <- sum(full_list[[i]]$stt_table[1, 5])
        ana_spec2[i] <- sum(full_list[[i]]$stt_table[1, 6])
      }
      immig_spec <- sum(immig_spec)
      ana_spec <- sum(ana_spec)
      immig_spec2 <- sum(immig_spec2)
      ana_spec2 <- sum(ana_spec2)
      init_present <- immig_spec + ana_spec + immig_spec2 + ana_spec2
      stt_all[1, 2:8] <- c(immig_spec, ana_spec, 0, immig_spec2, ana_spec2, 0, init_present)
    }else{
    #### Oceanic vs nonoceanic ####
      
      immig_spec <- c()
      ana_spec <- c()
      for (i in 1:M) {
        immig_spec[i] <- sum(full_list[[i]]$stt_table[1, 2])
        ana_spec[i] <- sum(full_list[[i]]$stt_table[1, 3])
      }
      immig_spec <- sum(immig_spec)
      ana_spec <- sum(ana_spec)
      init_present <- immig_spec + ana_spec
      stt_all[1, 2:5] <- c(immig_spec, ana_spec, 0, init_present)
    }
    
    
    #### 2 type ####
    if (number_type2_cols > 0) {
      # Type 1
      stt_list_type1 <- list()
      for (i in 1:max(which(type_vec == 1))) {
        stt_list_type1[[i]] <- full_list[[i]]$stt_table
      }
      stt_type1 <- create_full_CS_stt(
        stt_list = stt_list_type1,
        stac_vec = stac_vec,
        totaltime = totaltime,
        trait_pars = trait_pars
      )


      ######################################################### list type2
      type2len <- length(which(type_vec == 2))
      stt_list_type2 <- list()
      for (i in 1:type2len) {
        stt_list_type2[[i]] <- full_list[[which(type_vec == 2)[i]]]$stt_table
      }

      stt_type2 <- create_full_CS_stt(
        stt_list = stt_list_type2,
        stac_vec = stac_vec,
        totaltime = totaltime,
        trait_pars = trait_pars
      )

      island_list[[1]] <- list(island_age = totaltime,
                               not_present_type1 = DDD::roundn(
                                 M * (1 - prop_type2_pool)) -
                                 (number_type1_cols),
                               not_present_type2 = DDD::roundn(
                                 M * prop_type2_pool) - number_type2_cols,
                               stt_all = stt_all,
                               stt_type1 = stt_type1,
                               stt_type2 = stt_type2)
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
      island_list[[2]] <- list(
        branching_times = totaltime,
        stac = 0,
        missing_species = 0
      )
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

#' Unsampled CS full STT
#'
#' @param stt_list List of full stt tables as
#' returned by DAISIE_sim_core functions
#' @param totaltime Numeric double with total time of simulation.
#' @param stac_vec Vector with status of species on island.
#' @param trait_pars A named list containing diversification rates considering 
#' two trait states created by \code{\link{create_trait_pars}}:
#' \itemize{
#'   \item{[1]:A numeric with the per capita transition rate with state1}
#'   \item{[2]:A numeric with the per capita immigration rate with state2}
#'   \item{[3]:A numeric with the per capita extinction rate with state2}
#'   \item{[4]:A numeric with the per capita anagenesis rate with state2}
#'   \item{[5]:A numeric with the per capita cladogenesis rate with state2}
#'   \item{[6]:A numeric with the per capita transition rate with state2} 
#'   \item{[7]:A numeric with the number of species with trait state 2 on mainland} 
#' }
#'
#' @return 1 complete, unsampled STT table from all clades in an island of a
#' CS model as generated by DAISIE_sim_core functions.
#' @author Pedro Neves, Joshua Lambert, Shu Xie, Giovanni Laudanno
create_full_CS_stt <- function(stt_list, stac_vec, totaltime, trait_pars = NULL) {
  
  if(!is.null(trait_pars)){
        return(
          create_full_CS_stt_trait(
            stt_list = stt_list,
            stac_vec = stac_vec,
            totaltime = totaltime,
            trait_pars = trait_pars
          )
        )
      }
  # Return empty island, if empty
  present <- which(stac_vec != 0)
  
  # Checks if stt has only 2 rows and is empty at present (nothing happened)
  second_line_stts <- lapply(stt_list, "[", 2,)
  zeros_second_line <- sapply(second_line_stts, sum) == 0
  
  
  filled_stt_lists <- stt_list[!zeros_second_line]
  
  # Calculate 'present' and append to filled_stt_list
  # no_time_stts <- lapply(filled_stt_lists, "[", , 2:4)
  num_indep_colonists <- list()
  for (i in seq_along(filled_stt_lists)) {
    num_indep_colonists[[i]] <- filled_stt_lists[[i]][, 2] +
      filled_stt_lists[[i]][, 3] +
      filled_stt_lists[[i]][, 4]
    
    num_indep_colonists[[i]][which(num_indep_colonists[[i]] > 0)] <- 1
    filled_stt_lists[[i]] <- cbind(
      filled_stt_lists[[i]],
      present = num_indep_colonists[[i]]
    )
  }
  
  
  # If no colonization ever happened, just return 0 values
  if (length(filled_stt_lists) == 0) {
    times <- c(totaltime, 0)
    nI <- c(0, 0)
    nA <- c(0, 0)
    nC <- c(0, 0)
    diff_present <- c(0, 0)
  } else {
    
    deltas_matrix <- lapply(filled_stt_lists, FUN = diff)
    for (i in seq_along(filled_stt_lists)) {
      if (any(filled_stt_lists[[i]][1, ] !=
              c("Time" = totaltime, "nI" = 0, "nA" = 0, "nC" = 0, "present" = 0))) {
        deltas_matrix[[i]] <- rbind(
          filled_stt_lists[[i]][1, ],
          deltas_matrix[[i]]
        )
      } else {
        deltas_matrix[[i]] <- rbind(
          c("Time" = totaltime, "nI" = 0, "nA" = 0, "nC" = 0, "present" = 0),
          deltas_matrix[[i]]
        )
      }
    }
    
    times_list <- lapply(filled_stt_lists, "[", , 1) # nolint
    times <- unlist(times_list)
    
    nI_list <- lapply(deltas_matrix, "[", , 2) # nolint
    nA_list <- lapply(deltas_matrix, "[", , 3) # nolint
    nC_list <- lapply(deltas_matrix, "[", , 4) # nolint
    present_list <- lapply(deltas_matrix, "[", , 5) # nolint
    
    nI <- unlist(nI_list)
    nA <- unlist(nA_list)
    nC <- unlist(nC_list)
    diff_present <- unlist(present_list)
  }
  
  full_stt <- data.frame(
    times = times,
    nI = nI,
    nA = nA,
    nC = nC,
    present = diff_present
  )
  ordered_diffs <- full_stt[order(full_stt$times, decreasing = TRUE), ]
  
  complete_stt_table <- mapply(ordered_diffs[2:5], FUN = cumsum)
  complete_stt_table <- cbind(ordered_diffs$times, complete_stt_table)
  colnames(complete_stt_table) <- c("Time", "nI", "nA", "nC", "present")
  
  while (complete_stt_table[1, 1] == complete_stt_table[2, 1]) {
    complete_stt_table <- complete_stt_table[-1, ]
  }
  
  stt <- complete_stt_table
  # Remove final duplicate lines, if any
  while (
    all(stt[nrow(stt) - 1, ] == stt[nrow(stt), ])
  ) {
    stt <- stt[1:(nrow(stt) - 1), ]
  }
  return(stt)
}

create_full_CS_stt_trait <- function(stt_list, stac_vec, totaltime, trait_pars) {
  # Return empty island, if empty
  present <- which(stac_vec != 0)
  
  # Checks if stt has only 2 rows and is empty at present (nothing happened)
  second_line_stts <- lapply(stt_list, "[", 2,)
  zeros_second_line <- sapply(second_line_stts, sum) == 0
  filled_stt_lists <- stt_list[!zeros_second_line]
  
  # Calculate 'present' and append to filled_stt_list
  # no_time_stts <- lapply(filled_stt_lists, "[", , 2:4)
  num_indep_colonists <- list()
  for (i in seq_along(filled_stt_lists)) {
    num_indep_colonists[[i]] <- filled_stt_lists[[i]][, 2] +
      filled_stt_lists[[i]][, 3] +
      filled_stt_lists[[i]][, 4] +
      filled_stt_lists[[i]][, 5] + 
      filled_stt_lists[[i]][, 6] +
      filled_stt_lists[[i]][, 7]
    
    num_indep_colonists[[i]][which(num_indep_colonists[[i]] > 0)] <- 1
    filled_stt_lists[[i]] <- cbind(
      filled_stt_lists[[i]],
      present = num_indep_colonists[[i]]
    )
  }
  # If no colonization ever happened, just return 0 values
  if (length(filled_stt_lists) == 0) {
    times <- c(totaltime, 0)
    nI <- c(0, 0)
    nA <- c(0, 0)
    nC <- c(0, 0)
    nI2 <- c(0, 0)
    nA2 <- c(0, 0)
    nC2 <- c(0, 0)
    diff_present <- c(0, 0)
  } else {
    
    deltas_matrix <- lapply(filled_stt_lists, FUN = diff)
    for (i in seq_along(filled_stt_lists)) {
      if (any(filled_stt_lists[[i]][1, ] !=
              c("Time" = totaltime, "nI" = 0, "nA" = 0, "nC" = 0, "nI2" = 0, "nA2" = 0, "nC2" = 0, "present" = 0))) {
        deltas_matrix[[i]] <- rbind(
          filled_stt_lists[[i]][1, ],
          deltas_matrix[[i]]
        )
      } else {
        deltas_matrix[[i]] <- rbind(
          c("Time" = totaltime, "nI" = 0, "nA" = 0, "nC" = 0, "nI2" = 0, "nA2" = 0, "nC2" = 0, "present" = 0),
          deltas_matrix[[i]]
        )
      }
    }
    
    times_list <- lapply(filled_stt_lists, "[", , 1) # nolint
    all_times <- unlist(times_list)
    times <- all_times
    
    nI_list <- lapply(deltas_matrix, "[", , 2) # nolint
    nA_list <- lapply(deltas_matrix, "[", , 3) # nolint
    nC_list <- lapply(deltas_matrix, "[", , 4) # nolint
    nI2_list <- lapply(deltas_matrix, "[", , 5) # nolint
    nA2_list <- lapply(deltas_matrix, "[", , 6) # nolint
    nC2_list <- lapply(deltas_matrix, "[", , 7) # nolint
    present_list <- lapply(deltas_matrix, "[", , 8) #nolint
    
    nI <- unlist(nI_list)
    nA <- unlist(nA_list)
    nC <- unlist(nC_list)
    nI2 <- unlist(nI2_list)
    nA2 <- unlist(nA2_list)
    nC2 <- unlist(nC2_list)
    diff_present <- unlist(present_list)
  }
  
  full_stt <- data.frame(
    times = times,
    nI = nI,
    nA = nA,
    nC = nC,
    nI2 = nI2,
    nA2 = nA2,
    nC2 = nC2,
    present = diff_present
  )
  ordered_diffs <- full_stt[order(full_stt$times, decreasing = TRUE), ]
  
  complete_stt_table <- mapply(ordered_diffs[2:8], FUN = cumsum)
  complete_stt_table <- cbind(ordered_diffs$times, complete_stt_table)
  colnames(complete_stt_table) <- c("Time", "nI", "nA", "nC", "nI2", "nA2", "nC2", "present")
  
  while (complete_stt_table[1, 1] == complete_stt_table[2, 1]) {
    complete_stt_table <- complete_stt_table[-1, ]
  }
  
  stt <- complete_stt_table
  # Remove final duplicate lines, if any
  while (
    all(stt[nrow(stt) - 1, ] == stt[nrow(stt), ])
  ) {
    stt <- stt[1:(nrow(stt) - 1), ]
  }
  stt
}

