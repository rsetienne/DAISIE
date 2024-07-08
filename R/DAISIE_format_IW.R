#' Formats island-wide simulation output into
#' standard DAISIE list output
#'
#' @inheritParams default_params_doc
#'
#' @return List with CS DAISIE simulation output
#'
#' @noRd
DAISIE_format_IW <- function(island_replicates,
                             time,
                             M,
                             sample_freq = 25,
                             verbose = TRUE,
                             trait_pars = NULL) {

  if(!is.null(trait_pars)){
    return(
      DAISIE_format_IW_trait(
        island_replicates = island_replicates,
        time = time,
        M = M,
        sample_freq = sample_freq,
        verbose = verbose,
        trait_pars = trait_pars
      )
    )
  }

  total_time <- time
  testit::assert(
    !is.na(sample_freq) && !is.null(sample_freq) && sample_freq >= 1
  )

  if (is.infinite(sample_freq)) {
    several_islands <- DAISIE_format_IW_full_stt(
      island_replicates,
      total_time = total_time,
      M = M,
      verbose = verbose
    )
  } else {
    several_islands <- DAISIE_format_IW_sampled_stt(
      island_replicates = island_replicates,
      total_time = total_time,
      M = M,
      sample_freq = sample_freq,
      verbose = verbose
    )
  }

  return(several_islands)

}

DAISIE_format_IW_trait <- function(island_replicates,
                                   time,
                                   M,
                                   sample_freq,
                                   verbose = FALSE,
                                   trait_pars = NULL)
{

  total_time <- time
  several_islands = list()
  for(rep in 1:length(island_replicates))
  {
    the_island = island_replicates[[rep]]
    M1 <- M
    M2 <- trait_pars$M2
    Mtotal <- M1 + M2
    stt_all = matrix(ncol = 7,nrow = sample_freq + 1)
    colnames(stt_all) = c("Time","nI","nA","nC","nI2","nA2","nC2")

    stt_all[,"Time"] = rev(seq(from = 0,to = total_time, length.out = sample_freq + 1))
    stt_all[1, 2:7] = c(0, 0, 0, 0, 0, 0)


    the_stt = the_island$stt_table

    for(i in 2:nrow(stt_all))
    {
      the_age = stt_all[i,"Time"]
      stt_all[i,2:7] = the_stt[max(which(the_stt[,"Time"] >= the_age)),2:7]
    }
    island_list = list()

    if(sum(the_stt[nrow(the_stt),2:7]) == 0)
    {

      island_list[[1]] = list(
        island_age = total_time,
        not_present = Mtotal,
        stt_all = stt_all
      )

    } else {

      island_list[[1]] = list(island_age = total_time,
                              not_present = Mtotal - length(the_island$taxon_list),
                              stt_all = stt_all)

      for(y in 1:length(the_island$taxon_list))
      {
        island_list[[y + 1]] = the_island$taxon_list[[y]]
      }
    }

    island_list <- add_brt_table(island_list)

    several_islands[[rep]] = island_list
    if (verbose == TRUE) {
      message(
        "Island being formatted: ", rep, "/", length(island_replicates)
      )
    }
  }
  return(several_islands)
}

add_brt_table_old <- function(island, full_table = FALSE) {
  island_age <- island[[1]]$island_age
  island_top <- island[[1]]
  if (length(island) == 1) {
    brts_table <- matrix(ncol = 5, nrow = 1)
    brts_table[1, ] <-  c(island_age, 0, 0, NA, NA)
    island[[1]]$brts_table <- brts_table
  } else {
    island_top <- island[[1]]
    island[[1]] <- NULL
    btimes <- list()
    for (i in 1:length(island)) {
      btimes[[i]] <- island[[i]]$branching_times[-1]
    }
    island <- island[rev(order(sapply(btimes, "[", 1)))]
    il <- unlist(island)
    stac1s <- which(il[which(names(il) == "stac")] == "1")
    stac5s <- which(il[which(names(il) == "stac")] == "5")
    stac1_5s <- sort(c(stac1s, stac5s))
    if (length(stac1_5s) != 0) {
      if (length(stac1_5s) == length(island)) {
        brts_table <- matrix(ncol = 5, nrow = 1)
        brts_table[1, ] <- c(island_age, 0, 0, NA, NA)
        island_no_stac1or5 <- NULL
      } else {
        island_no_stac1or5 <- island[-stac1_5s]
      }
    }
    if (length(stac1_5s) == 0) {
      island_no_stac1or5 <- island
    }
    if (length(island_no_stac1or5) != 0) {
      btimes <- list()
      for (i in 1:length(island_no_stac1or5)) {
        btimes[[i]] <- island_no_stac1or5[[i]]$branching_times[-1]
      }
      brts <- rev(sort(unlist(btimes)))
      brts_IWK <- NULL
      pos1 <- 0
      j <- 1
      for (i in 1:length(btimes)) {
        the_brts <- btimes[[i]]
        the_stac <- island_no_stac1or5[[i]]$stac
        pos2 <- pos1 + length(the_brts)
        ff <- matrix(ncol = 5, nrow = pos2 - pos1)
        ff[1:(pos2 - pos1), 1] <- the_brts
        ff[1:(pos2 - pos1), 2] <- i
        ff[1:(pos2 - pos1), 3] <- seq(1, length(the_brts))
        ff[1:(pos2 - pos1), 4] <- (the_stac == 2) +
          (the_stac == 3) + (the_stac == 4) * 0
        ff[1:(pos2 - pos1), 5] <- NA
        brts_IWK <- rbind(brts_IWK,ff)
        pos1 <- pos2
        j <- j + 1
        if( !is.null(island[[i]]$all_colonisations) & full_table == TRUE) {
          for (k in 1:length(island[[i]]$all_colonisations)) {
            the_brts <- island[[i]]$all_colonisations[[k]]$event_times[-1]
            pos2 <- pos1 + length(the_brts)
            ff <- matrix(ncol = 5, nrow = pos2 - pos1)
            ff[1:(pos2 - pos1), 1] <- the_brts
            ff[1:(pos2 - pos1), 2] <- j
            ff[1:(pos2 - pos1), 3] <- seq(1, length(the_brts))
            ff[1:(pos2 - pos1), 4] <- NA
            ff[1:(pos2 - pos1), 5] <- j - 1
            brts_IWK <- rbind(brts_IWK,ff)
            pos1 <- pos2
            j <- j + 1
          }
        }
      }
      brts_table <- brts_IWK[rev(order(brts_IWK[, 1])), ]
      brts_table <- rbind(c(island_age, 0, 0, NA, NA), brts_table)
    }
    island_top$brts_table <- brts_table
    if (length(stac1_5s) != 0) {
      for (i in 1:length(stac1_5s)) {
        island[[length(island) + 1]] <- island[[stac1_5s[i]]]
        island[[stac1_5s[i]]] <- NULL
        stac1_5s <- stac1_5s - 1
      }
    }
    island <- append(list(island_top), island)
  }
  colnames(island[[1]]$brts_table) <- c("brt", "clade", "event", "endemic", "col")
  return(island)
}

add_brt_table <- function(island, full_table = TRUE) {
  island_age <- island[[1]]$island_age
  island_top <- island[[1]]
  if (length(island) == 1) {
    brts_table <- matrix(ncol = 5, nrow = 1)
    brts_table[1, ] <-  c(island_age, 0, 0, NA, NA)
    island[[1]]$brts_table <- brts_table
  } else {
    island_top <- island[[1]]
    island[[1]] <- NULL
    btimes <- list()
    for (i in 1:length(island)) {
      btimes[[i]] <- island[[i]]$branching_times[-1]
    }
    island <- island[rev(order(sapply(btimes, "[", 1)))]
    il <- unlist(island)
    stac1s <- which(il[which(names(il) == "stac")] == "1")
    stac5s <- which(il[which(names(il) == "stac")] == "5")
    stac1_5s <- sort(c(stac1s, stac5s))
    if (length(stac1_5s) != 0) {
      if (length(stac1_5s) == length(island)) {
        brts_table <- matrix(ncol = 5, nrow = 1)
        brts_table[1, ] <- c(island_age, 0, 0, NA, NA)
        island_no_stac1or5 <- NULL
      } else {
        island_no_stac1or5 <- island[-stac1_5s]
      }
    }
    if (length(stac1_5s) == 0) {
      island_no_stac1or5 <- island
    }
    if (length(island_no_stac1or5) != 0) {
      btimes <- list()
      for (i in 1:length(island_no_stac1or5)) {
        btimes[[i]] <- island_no_stac1or5[[i]]$branching_times[-1]
      }
      brts <- rev(sort(unlist(btimes)))
      brts_IWK <- NULL
      pos1 <- 0
      for (i in 1:length(btimes)) {
        the_stac <- island_no_stac1or5[[i]]$stac
        if(!is.null(island[[i]]$all_colonisations) & full_table == TRUE) {
          #if(length(island[[i]]$all_colonisations) > 0) print(i)
          for (k in 1:length(island[[i]]$all_colonisations)) {
            the_brts <- island[[i]]$all_colonisations[[k]]$event_times[-1]
            pos2 <- pos1 + length(the_brts)
            ff <- matrix(ncol = 5, nrow = pos2 - pos1)
            ff[1:(pos2 - pos1), 1] <- the_brts
            ff[1:(pos2 - pos1), 2] <- i
            ff[1:(pos2 - pos1), 3] <- seq(1, length(the_brts))
            ff[1:(pos2 - pos1), 4] <- (the_stac == 2) + (the_stac == 3)
            ff[1:(pos2 - pos1), 5] <- k
            brts_IWK <- rbind(brts_IWK,ff)
            pos1 <- pos2
          }
        } else {
          the_brts <- btimes[[i]]
          pos2 <- pos1 + length(the_brts)
          ff <- matrix(ncol = 5, nrow = pos2 - pos1)
          ff[1:(pos2 - pos1), 1] <- the_brts
          ff[1:(pos2 - pos1), 2] <- i
          ff[1:(pos2 - pos1), 3] <- seq(1, length(the_brts))
          ff[1:(pos2 - pos1), 4] <- (the_stac == 2) + (the_stac == 3)
          ff[1:(pos2 - pos1), 5] <- 1
          brts_IWK <- rbind(brts_IWK,ff)
          pos1 <- pos2
        }
      }
      brts_table <- brts_IWK[rev(order(brts_IWK[, 1])), ]
      brts_table <- rbind(c(island_age, 0, 0, NA, NA), brts_table)
    }
    island_top$brts_table <- brts_table
    if (length(stac1_5s) != 0) {
      for (i in 1:length(stac1_5s)) {
        island[[length(island) + 1]] <- island[[stac1_5s[i]]]
        island[[stac1_5s[i]]] <- NULL
        stac1_5s <- stac1_5s - 1
      }
    }
    island <- append(list(island_top), island)
  }
  colnames(island[[1]]$brts_table) <- c("brt", "clade", "event", "endemic", "col")
  return(island)
}
