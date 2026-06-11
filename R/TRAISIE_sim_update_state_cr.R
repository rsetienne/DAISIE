#' @keywords internal
TRAISIE_sim_update_state_cr <- function(timeval,
                                        total_time,
                                        possible_event,
                                        maxspecID,
                                        island_spec,
                                        stt_table,
                                        trait_pars,
                                        num_observed_states,
                                        num_hidden_states,
                                        mainland) {

  n <- num_observed_states * num_hidden_states
  p <- trait_pars$p

  ##########################################
  #IMMIGRATION
  for (i in 0:(n - 1)) {

    # Calculate the corresponding trait state
    trait_state <- i + 1  # Maps event to trait state dynamically

    # IMMIGRATION (4*i + 1)
    if (possible_event == (4 * i + 1)) {
      mainland_total <- sum(unlist(mainland))

      # Determine which trait block the species belongs to
      block <- ceiling(trait_state / num_hidden_states)

      # Compute start and end indices for that block
      start_idx <- if (block == 1) 1 else sum(unlist(mainland[1:(block - 1)])) + 1
      end_idx   <- sum(unlist(mainland[1:block]))

      # Sample a colonist from the block
      colonist <- DDD::sample2(start_idx:end_idx, 1)




      if (length(island_spec[, 1]) != 0) {
        isitthere = which(island_spec[, 1] == colonist)
      } else {
        isitthere = c()
      }

      # Check if species is already present
      testit::assert(length(isitthere) <= 1)
      if (length(isitthere) == 0) {
        island_spec = rbind(island_spec, c(colonist, colonist, timeval,
                                           "I", NA, NA, NA, trait_state, NA))
      }

      if (length(isitthere) != 0) {
        island_spec[isitthere, ] = c(colonist, colonist, timeval,
                                     "I", NA, NA, NA, trait_state, NA)
      }
    }

    # EXTINCTION (4*i + 2)
    if (possible_event == (4 * i + 2)) {
      island_spec_state = which(island_spec[, 8] == as.character(trait_state))


      extinct = DDD::sample2(island_spec_state, 1)

      typeofspecies = island_spec[extinct, 4]

      if (typeofspecies == "I") {
        island_spec = island_spec[-extinct, ]
      }

      if (typeofspecies == "A") {
        island_spec = island_spec[-extinct, ]
      }

      if (typeofspecies == "C") {
        sisters = intersect(which(island_spec[, 2] == island_spec[extinct, 2]),
                            which(island_spec[, 3] == island_spec[extinct, 3]))
        survivors = sisters[which(sisters != extinct)]

        if (length(sisters) == 2) {
          island_spec[survivors, 4] = "A"
          island_spec[survivors, c(5, 6)] = c(NA, NA)
          island_spec[survivors, 7] = "Clado_extinct"
          island_spec = island_spec[-extinct, ]
        }

        if (length(sisters) >= 3) {
          numberofsplits = nchar(island_spec[extinct, 5])
          mostrecentspl = substring(island_spec[extinct, 5], numberofsplits)

          sistermostrecentspl = ifelse(mostrecentspl == "B", "A", "B")
          motiftofind = paste(substring(island_spec[extinct, 5], 1, numberofsplits - 1), sistermostrecentspl, sep = "")
          possiblesister = survivors[which(substring(island_spec[survivors, 5], 1, numberofsplits) == motiftofind)]

          if (mostrecentspl == "A") {
            tochange = possiblesister[which(island_spec[possiblesister, 6] == min(as.numeric(island_spec[possiblesister, 6])))]
            island_spec[tochange, 6] = island_spec[extinct, 6]
          }

          island_spec[possiblesister, 5] = paste(substring(island_spec[possiblesister, 5], 1, numberofsplits - 1),
                                                 substring(island_spec[possiblesister, 5], numberofsplits + 1,
                                                           nchar(island_spec[possiblesister, 5])), sep = "")
          island_spec = island_spec[-extinct, ]
        }
      }
      island_spec = rbind(island_spec)
    }

    # ANAGENESIS (4*i + 3)
    if (possible_event == (4 * i + 3)) {
      immi_specs = intersect(which(island_spec[, 4] == "I"), which(island_spec[, 8] == as.character(trait_state)))

      if (length(immi_specs) == 0) next
      if (length(immi_specs) == 1) {
        anagenesis = immi_specs
      }

      if (length(immi_specs) > 1) {
        anagenesis = DDD::sample2(immi_specs, 1)
      }

      # Step 1: get the value of column 2 for the anagenesis row
      col2_value <- island_spec[anagenesis, 2]

      # Step 2: select all rows that have the same mailand ancestor
      rows_same_col2 <- island_spec[island_spec[, 2] == col2_value, , drop = FALSE]

      # Step 3: check if the immigrant species sampled (anagenesis) is unique of if it is a recolonist

      if (length(unique(rows_same_col2[, 3])) == 1) {

        maxspecID = maxspecID + 1
        island_spec[anagenesis, 4] = "A"
        island_spec[anagenesis, 1] = maxspecID
        island_spec[anagenesis, 7] = "Immig_parent"
        if (!is.null(trait_pars)) {
          island_spec[anagenesis, 8] = as.character(trait_state)
        }
      } else {
        maxspecID = maxspecID + 1
        island_spec[anagenesis, 4] = "A"
        island_spec[anagenesis, 1] = maxspecID
        island_spec[anagenesis, 7] = "Immig_parent"
        if (!is.null(trait_pars)) {
          island_spec[anagenesis, 8] = as.character(trait_state)
        }

        ### check why this
        island_spec = rbind(island_spec, island_spec[anagenesis, ])
        island_spec = island_spec[-anagenesis, ]

      }


    }

    # CLADOGENESIS (4*i + 4)
    if (possible_event == (4 * i + 4)) {

      island_spec_state = which(island_spec[, 8] == as.character(trait_state))

      tosplit = DDD::sample2(island_spec_state, 1)


      if (island_spec[tosplit, 4] == "C") {
        island_spec[tosplit, 4] = "C"
        island_spec[tosplit, 1] = maxspecID + 1
        oldstatus = island_spec[tosplit, 5]
        island_spec[tosplit, 5] = paste(oldstatus, "A", sep = "")
        island_spec[tosplit, 7] = NA
        island_spec[tosplit, 8] = as.character(trait_state)
        oldsplit = island_spec[tosplit, 9]
        split_time <- total_time - as.numeric(timeval)
        island_spec[tosplit, 9] = paste(as.character(oldsplit), as.character(split_time), sep = " ")

        island_spec = rbind(island_spec, c(maxspecID + 2, island_spec[tosplit, 2], island_spec[tosplit, 3],
                                           "C", paste(oldstatus, "B", sep = ""), timeval, NA, trait_state, paste(as.character(oldsplit), as.character(split_time), sep = " ")))
        maxspecID = maxspecID + 2
      } else {
        island_spec[tosplit, 4] = "C"
        island_spec[tosplit, 1] = maxspecID + 1
        island_spec[tosplit, 5] = "A"
        island_spec[tosplit, 6] = island_spec[tosplit, 3]
        island_spec[tosplit, 7] = NA
        island_spec[tosplit, 8] = as.character(trait_state)
        split_time <- total_time - as.numeric(timeval)
        oldsplit = island_spec[tosplit, 9]
        island_spec[tosplit, 9] = paste(as.character(oldsplit), as.character(split_time), sep = " ")

        island_spec = rbind(island_spec, c(maxspecID + 2, island_spec[tosplit, 2], island_spec[tosplit, 3],
                                           "C", "B", timeval, NA, trait_state, paste(as.character(oldsplit), as.character(split_time), sep = " ")))
        maxspecID = maxspecID + 2
      }




    }

  }
  ## trait change

  #######transition rate

  # assume events 1:4n are the birth/death/speciation ones…

  # TRAIT CHANGE
  for (i in 0:(n - 1)) {
    for (j in 1:n) if (j != i + 1) {
      # now each (i,j) pair maps to 4*n + (i * n + j)
      event_idx <- 4 * n + i * n + j

      if (possible_event == event_idx) {
        # pick one species in trait state (i+1)
        island_spec_state1 <- which(island_spec[, 8] == as.character(i + 1))
        if (length(island_spec_state1) > 0) {
          totrans <- DDD::sample2(island_spec_state1, 1)
          # optional: convert immigrant → endemic if p==1, etc.
          if (p == 1 && island_spec[totrans, 4] == "I") {
            island_spec[totrans, 4] <- "A"
          }
          # finally update the trait

          island_spec[totrans, 8] <- as.character(j)
        }
      }
    }

  }


  if (total_time >= timeval) {
    # Initialize a vector to hold counts for each state


    # Loop through each trait state (1 to n)
    for (state in 1:n) {
      # Count the number of species in each category (I, A, C) for each state


      nI <- length(intersect(which(island_spec[, 4] == "I"), which(island_spec[, 8] == as.character(state))))  # Immigrant species in state
      nA <- length(intersect(which(island_spec[, 4] == "A"), which(island_spec[, 8] == as.character(state))))  # Ancestor species in state
      nC <- length(intersect(which(island_spec[, 4] == "C"), which(island_spec[, 8] == as.character(state))))  # Colonist species in state



      # Append the counts for this state to the counts vector

    }

    # Update the stt_table with the new counts for each state
    stt_table <- rbind(stt_table,
                       c(total_time - timeval,
                         length(which(island_spec[, 4] == "I")),
                         length(which(island_spec[, 4] == "A")),
                         length(which(island_spec[, 4] == "C"))))
  }

  # Add a final row of zeros

  updated_state <- list(island_spec = island_spec,
                        maxspecID = maxspecID,
                        stt_table = stt_table)
  return(updated_state)
}
