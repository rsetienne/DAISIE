#' Updates state of island given sampled event
#'
#' Makes the event happen by updating island species matrix and species IDs.
#' What event happens is determined by the sampling in the algorithm.
#'
#' @inheritParams default_params_doc
#'
#' @return The updated state of the system, which is a list with the
#' \code{island_spec} matrix, an integer \code{maxspecID} with the most recent
#' ID of species and the \code{stt_table}, a matrix with the current species
#' through time table.
#'
#' @seealso \link{DAISIE_sim_core_time_dependent}
DAISIE_sim_update_state_time_dependent <- function(timeval,
                                                   totaltime,
                                                   possible_event,
                                                   maxspecID,
                                                   mainland_spec,
                                                   island_spec,
                                                   stt_table,
                                                   rates,
                                                   max_rates) {
  if (is.null(max_rates)) {
    max_rates <- list(
      ext_max_rate = rates$ext_rate,
      immig_max_rate = rates$immig_rate,
      ana_max_rate = rates$ana_rate,
      clado_max_rate = rates$clado_rate
    )
  }

  event <- "false"

  ##########################################
  #IMMIGRATION
  if (possible_event == 1) {
    if (rates$immig_rate / max_rates$immig_max_rate >= stats::runif(1)) {
      event <- "real"
      colonist <- DDD::sample2(mainland_spec, 1)
      if (length(island_spec[, 1]) != 0) {
        isitthere <- which(island_spec[, 1] == colonist)
      } else {
        isitthere <- c()
      }
      if (length(isitthere) == 0) {
        island_spec <- rbind(island_spec, c(colonist, colonist, timeval, "I", NA, NA, NA))
      }
      if (length(isitthere) != 0) {
        island_spec[isitthere, ] <- c(colonist, colonist, timeval, "I", NA, NA, NA)
      }
    }
  }
  ##########################################
  #EXTINCTION
  if (possible_event == 2) {
    if (rates$ext_rate / max_rates$ext_max_rate >= stats::runif(1)) {
      event <- "real"
      extinct <- DDD::sample2(1:length(island_spec[, 1]), 1)
      #this chooses the row of species data to remove
      typeofspecies <- island_spec[extinct, 4]
      if (typeofspecies == "I") {
        island_spec <- island_spec[-extinct, ]
      }
      #remove immigrant
      if (typeofspecies == "A") {
        island_spec <- island_spec[-extinct, ]
      }
      #remove anagenetic
      if (typeofspecies == "C") {
        #remove cladogenetic
        #first find species with same ancestor AND arrival totaltime
        sisters <- intersect(which(island_spec[, 2] == island_spec[extinct, 2]), which(island_spec[, 3] == island_spec[extinct, 3]))
        survivors <- sisters[which(sisters != extinct)]
        if (length(sisters) == 2) {
          #survivors status becomes anagenetic
          island_spec[survivors, 4] <- "A"
          island_spec[survivors, c(5, 6)] <- c(NA, NA)
          island_spec[survivors, 7] <- "Clado_extinct"
          island_spec <- island_spec[-extinct, ]
        }

        if (length(sisters) >= 3) {
          numberofsplits <- nchar(island_spec[extinct, 5])
          mostrecentspl <- substring(island_spec[extinct, 5], numberofsplits)

          if (mostrecentspl == "B") {
            sistermostrecentspl <- "A"
          }
          if (mostrecentspl == "A") {
            sistermostrecentspl <- "B"
          }
          motiftofind <- paste(substring(island_spec[extinct, 5], 1, numberofsplits - 1), sistermostrecentspl, sep = "")
          possiblesister <- survivors[which(substring(island_spec[survivors, 5], 1, numberofsplits) == motiftofind)]
          #different rules depending on whether a B or A is removed. B going extinct is simpler because it only
          #carries a record of the most recent speciation
          if (mostrecentspl == "A") {
            #change the splitting date of the sister species so that it inherits the early splitting that used to belong to A.
            # Bug fix here thanks to Nadiah Krisensen: max -> min
            tochange <- possiblesister[which(island_spec[possiblesister, 6] == min(as.numeric(island_spec[possiblesister, 6])))]
            island_spec[tochange, 6] <- island_spec[extinct, 6]
          }
          #remove the offending A/B from these species
          island_spec[possiblesister, 5] <- paste(substring(island_spec[possiblesister, 5], 1, numberofsplits - 1),
                                                  substring(island_spec[possiblesister, 5], numberofsplits + 1,
                                                            nchar(island_spec[possiblesister, 5])), sep = "")
          island_spec <- island_spec[-extinct, ]
        }
      }
      island_spec <- rbind(island_spec)
    }
  }
  ##########################################
  #ANAGENESIS
  if (possible_event == 3) {
    event <- "real"
    immi_specs <- which(island_spec[, 4] == "I")
    #we only allow immigrants to undergo anagenesis
    if(length(immi_specs) == 1) {
      anagenesis <- immi_specs
    }
    if (length(immi_specs) > 1) {
      anagenesis <- DDD::sample2(immi_specs, 1)
    }
    maxspecID <- maxspecID + 1
    island_spec[anagenesis, 4] <- "A"
    island_spec[anagenesis, 1] <- maxspecID
    island_spec[anagenesis, 7] <- "Immig_parent"
  }
  ##########################################
  #CLADOGENESIS - this splits species into two new species - both of which receive
  if (possible_event == 4) {
    if (rates$clado_rate / max_rates$clado_max_rate >= stats::runif(1)) {
      event <- "real"
      tosplit <- DDD::sample2(1:length(island_spec[, 1]), 1)
      #if the species that speciates is cladogenetic
      if (island_spec[tosplit, 4] == "C") {
        #for daughter A
        island_spec[tosplit, 4] <- "C"
        island_spec[tosplit, 1] <- maxspecID + 1
        oldstatus <- island_spec[tosplit, 5]
        island_spec[tosplit, 5] <- paste(oldstatus, "A", sep = "")
        #island_spec[tosplit,6] = timeval
        island_spec[tosplit, 7] <- NA
        #for daughter B
        island_spec <- rbind(island_spec, c(maxspecID + 2, island_spec[tosplit, 2], island_spec[tosplit, 3],
                                            "C", paste(oldstatus, "B", sep = ""), timeval, NA))
        maxspecID <- maxspecID + 2
      } else {
        #if the species that speciates is not cladogenetic
        #for daughter A
        island_spec[tosplit, 4] <- "C"
        island_spec[tosplit, 1] <- maxspecID + 1
        island_spec[tosplit, 5] <- "A"
        island_spec[tosplit, 6] <- island_spec[tosplit, 3]
        island_spec[tosplit, 7] <- NA
        #for daughter B
        island_spec <- rbind(island_spec, c(maxspecID + 2, island_spec[tosplit, 2], island_spec[tosplit, 3], "C", "B", timeval, NA))
        maxspecID <- maxspecID + 2
      }
    }
  }
  if (event == "real") {
    stt_table <- rbind(stt_table,
                       c(totaltime - timeval,
                         length(which(island_spec[, 4] == "I")),
                         length(which(island_spec[, 4] == "A")),
                         length(which(island_spec[, 4] == "C"))))
  }

  updated_state <- list(island_spec = island_spec,
                        maxspecID = maxspecID,
                        stt_table = stt_table)
  return(updated_state)
}
