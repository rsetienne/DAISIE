#' Updates state of island given sampled event for a constant rate case.
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
#' @keywords internal
#'
#' @seealso \link{DAISIE_sim_core_cr},
#' \link{DAISIE_sim_update_state_cr}
DAISIE_sim_update_state_cr <- function(timeval,
                                       total_time,
                                       possible_event,
                                       maxspecID,
                                       mainland_spec,
                                       island_spec,
                                       stt_table) {


  ##########################################
  #IMMIGRATION
  if (possible_event == 1) {
    colonist <- DDD::sample2(mainland_spec, 1)
    if (length(island_spec[, 1]) != 0) {
      isitthere <- which(island_spec[, 1] == colonist)
    } else {
      isitthere <- c()
    }
    testit::assert(length(isitthere) <= 1)
    if (length(isitthere) == 0) {
      island_spec <- rbind(island_spec, c(colonist, colonist, timeval, "I", NA, NA, NA))
    }
    if (length(isitthere) != 0) {
      island_spec[isitthere, ] <- c(colonist, colonist, timeval, "I", NA, NA, NA)
    }
  }
  ##########################################
  #EXTINCTION
  if (possible_event == 2) {
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
        #first find species with same ancestor AND arrival total_time
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
            # Bug fix here thanks to Nadiah Kristensen: max -> min
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
  ##########################################
  #ANAGENESIS
  if (possible_event == 3) {
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
  stt_table <- rbind(stt_table,
                     c(total_time - timeval,
                       length(which(island_spec[, 4] == "I")),
                       length(which(island_spec[, 4] == "A")),
                       length(which(island_spec[, 4] == "C"))))

  updated_state <- list(island_spec = island_spec,
                        maxspecID = maxspecID,
                        stt_table = stt_table)
  return(updated_state)
}
