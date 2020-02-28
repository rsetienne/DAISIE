#' Counts the number of species
#'
#' @param datalistelement stub
#'
#' @return A numeric
#' @export
countspecies <- function(datalistelement) {
  N <- length(datalistelement$branching_times) -
    1 + datalistelement$missing_species
  N
}

#' Counts the number of type 1 species
#'
#' @param datalistelement stub
#'
#' @return  something
#' @export
counttype1 <- function(datalistelement) {
  N1 <- 0
  if (length(datalistelement$type1or2) > 0) {
    N1 <- (datalistelement$type1or2 == 1)
    N1
  }
}

#' Title
#'
#' @param datalistelement stub
#'
#' @return  something
#' @export
countspeciestype1 <- function(datalistelement) {
  N1 <- 0
  if (length(datalistelement$type1or2) > 0) {
    if (datalistelement$type1or2 == 1) {
      N1 <- length(datalistelement$branching_times) -
        1 + datalistelement$missing_species
    }
  }
}

#' Title
#'
#' @param datalistelement stub
#'
#' @return  something
#' @export
countimmi <- function(datalistelement) {
  datalistelement$stac != 2
}

#' Title
#'
#' @param datalistelement stub
#' @param stac stub
#'
#' @return  something
#' @export
countstac <- function(datalistelement, stac) {
  return(datalistelement$stac == stac)
}

fconstr13 <- function(x, pars1, x_E, age) {
  lac <- pars1[1]
  laa <- pars1[5]
  ga <- pars1[4]
  A <- x - lac
  C <- ga + laa + 2 * lac
  ff <- (1 + A / C * (1 - exp(-C * age))) * exp(-A * age) - (1 - x_E)
  return(ff)
}

fconstr15 <- function(x, pars1, x_E, x_I, age) {
  lac <- pars1[1]
  laa <- pars1[5]
  A <- x - lac
  B_c <- -1 / age * log(1 - x_I)
  ga <- B_c - x - laa - lac
  C <- ga + laa + 2 * lac
  ff <- (1 + A / C * (1 - exp(-C * age))) * exp(-A * age) - (1 - x_E)
  return(ff)
}

calcMN <- function(datalist, pars1) {
  N <- sum(unlist(lapply(datalist, countspecies)))
  if (is.null(datalist[[1]]$not_present)) {
    M <- datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2 +
      length(datalist) - 1
    if (!is.na(pars1[6])) {
      if (is.na(pars1[11])) {
        M <- datalist[[1]]$not_present_type1 +
          sum(unlist(lapply(datalist, counttype1)))
      } else {
        M <- M - max(0, DDD::roundn(pars1[11] * M))
      }
      N <- sum(unlist(lapply(datalist, countspeciestype1)))
    }
  } else {
    M <- datalist[[1]]$not_present + length(datalist) - 1
  }
  return(c(M, N))
}

DAISIE_eq <- function(datalist, pars1, pars2) {
  eqmodel <- pars2[5]
  ddep <- pars2[2]
  MN <- calcMN(datalist, pars1)
  M <- MN[1]
  N <- MN[2]
  I <- sum(unlist(lapply(datalist, countimmi)))
  rNM <- N / M
  rIM <- I / (M - I)
  rIN <- I / (N - I)
  clado <- pars1[1] * ((1 - N / pars1[3]) ^ (ddep == 1 || ddep == 11)) *
    (exp(-N / pars1[3])) ^ (ddep == 2 || ddep == 21)
  ana <- pars1[5]
  # Equilibrium based on deterministic model in terms of N
  if (eqmodel == 1) {
    immi <- pars1[4] * ((1 - N / pars1[3]) ^ (ddep == 11)) *
      (exp(-N / pars1[3])) ^ (ddep == 21)
    ext <- clado + immi * (1 / rNM - 1)
    pars1[2] <- ext
  }
  # Equilibrium model based on deterministic model in terms of E and I
  if (eqmodel == 2) { # Only eq for N
    ext <- pars1[2]
    immitot <- 1 / (1 / rNM * 1 / (ext - clado) - 1 / (ana + clado + ext))
    immi <- immitot / ((1 - N / pars1[3]) ^ (ddep == 11) *
                         (exp(-N / pars1[3])) ^ (ddep == 21))
    pars1[4] <- immi
  }
  if (eqmodel == 3) { # Only eq for E
    immi <- pars1[4] * ((1 - N / pars1[3]) ^ (ddep == 11)) *
      (exp(-N / pars1[3])) ^ (ddep == 21)
    ext <- clado + (ana + 2 * clado) * rIN
    pars1[2] <- ext
  }
  if (eqmodel == 4) { # Only eq for I
    ext <- pars1[2]
    immitot <- (ext + ana + clado) * rIM
    immi <- immitot / ((1 - N / pars1[3]) ^ (ddep == 11) *
                         (exp(-N / pars1[3])) ^ (ddep == 21))
    pars1[4] <- immi
  }
  if (eqmodel == 5) { # Eq for E and I
    ext <- clado + (ana + 2 * clado) * rIN
    immitot <- (ext + ana + clado) * rIM
    immi <- immitot / ((1 - N / pars1[3]) ^ (ddep == 11) *
                         (exp(-N / pars1[3])) ^ (ddep == 21))
    pars1[2] <- ext
    pars1[4] <- immi
  }
  # Within x_E of equilibrium for E - diversity-dependence not implemented
  if (eqmodel == 13) {
    x_E <- pars2[10]
    x_I <- pars2[11]
    age <- datalist[[1]]$island_age
    pars1[2] <- stats::uniroot(f = fconstr13,
                               interval = c(pars1[1] + 1E-6, pars1[1] + 10),
                               pars1 = pars1,
                               x_E = x_E,
                               age = age)$root
    ga_c <- -1 / age * log(1 - x_I) - pars1[1] - pars1[2] - pars1[5]
    if (pars1[4] < ga_c) {
      cat("The non-endemics do not satisfy the equilibrium criterion for
                these parameters.\n")
    }
  }
  #Within x_E and x_I of equilibrium for both E and
  #I - diversity-dependence not implemented
  if (eqmodel == 15) {
    x_E <- pars2[10]
    x_I <- pars2[11]
    age <- datalist[[1]]$island_age
    pars1[2] <- stats::uniroot(f = fconstr15, interval = c(pars1[1] + 1E-6, pars1[1] + 10), pars1 = pars1, x_E = x_E, x_I = x_I, age = age)$root
    pars1[4] <- -1 / age * log(1 - x_I) - pars1[1] - pars1[2] - pars1[5]
  }
  return(pars1)
}

antidiagSums <- function(mat) {
  dime <- dim(mat)
  out <- rep(0, sum(dime) - 1)
  nr <- nrow(mat)
  nc <- ncol(mat)
  for (i in 1:(nr + nc - 1)) {
    rownums <- min(i, nr):max(1, i - nc + 1)
    colnums <- max(1, i - nr + 1):min(i, nc)
    for (j in seq_along(rownums)) {
      out[i] <- out[i] + mat[rownums[j], colnums[j]]
    }
  }
  return(out)
}

#' Translate user-friendly ontogeny codes to numerics
#'
#' @inherit DAISIE_sim_time_dependent
#'
#' @return Numeric, 0 for null-ontogeny, 1 for beta function
#' @export
#' @examples translated_ontogeny <- translate_island_ontogeny("const")
translate_island_ontogeny <- function(island_ontogeny) {

  if (island_ontogeny == "const" || island_ontogeny == 0) {
    island_ontogeny <- 0
  }

  if (island_ontogeny == "beta" || island_ontogeny == 1) {
    island_ontogeny <- 1
  }
  return(island_ontogeny)
}

#' Translate user-friendly sea-level codes to numerics
#'
#' @inherit DAISIE_sim_time_dependent
#'
#' @return Numeric, 0 for null-sea-level, 1 for sine function
#' @export
#' @examples tanslated_sea_level <- translate_sea_level("const")
translate_sea_level <- function(sea_level) {

  if (sea_level == "const" || sea_level == 0) {
    sea_level <- 0
  }

  if (sea_level == "sine" || sea_level == 1) {
    sea_level <- 1
  }
  return(sea_level)
}


#' Determine if list has only numerical values.
#'
#'
#' @param x Object to determine
#'
#' @return Boolean indicating if object is list with only numerical values
#' @note do not forget: NAs are removed from a list!
#' @examples
#'   testit::assert(
#'     DAISIE:::is_numeric_list(
#'       x = list(char = "character", numerical = 1)
#'     ) == FALSE
#'   )
#'
#'   testit::assert(
#'     DAISIE:::is_numeric_list(
#'       x = list(numerical_1 = 1, numerical_2 = 2)
#'     ) == TRUE
#'   )
is_numeric_list <- function(x) {
  is.list(x) && is.numeric(unlist(x))
}

#' Calculates the species on the island initially when \code{nonoceanic_pars[1]
#' != 0}
#'
#' @param prob_samp probability of a species being sampled
#' from the mainland pool
#' @param prob_nonend probability of a species sampled being non-endemic
#' @param mainland_n number of species in the mainland pool
#'
#' @return A list of non-endemic species, endemic species and the new
#' mainland species pool
#' @export
#'
#' @examples DAISIE_nonoceanic_spec(
#' prob_samp = 0.1,
#' prob_nonend = 0.9,
#' mainland_n = 1000)
DAISIE_nonoceanic_spec <- function(prob_samp, prob_nonend, mainland_n) {
  testit::assert(prob_samp <= 1)
  testit::assert(prob_samp >= 0)
  testit::assert(prob_nonend <= 1)
  testit::assert(prob_nonend  >= 0)
  testit::assert(length(mainland_n) > 0)
  if (prob_samp != 0) {
    prob_not_samp <- 1 - prob_samp
    prob_nonend <- prob_samp * prob_nonend
    prob_end <- 1 - (prob_not_samp + prob_nonend)
    num_native_spec <- sample(1:3, length(1:mainland_n),
                              replace = TRUE,
                              prob = c(prob_not_samp, prob_nonend, prob_end))
    init_nonend_spec_vec <- sample(1:mainland_n,
                                   length(which(num_native_spec == 2)),
                                   replace = FALSE)
    new_source_pool <- setdiff(1:mainland_n, init_nonend_spec_vec)
    init_end_spec_vec <- sample(new_source_pool,
                                length(which(num_native_spec == 3)),
                                replace = FALSE)
    mainland_spec <- setdiff(1:mainland_n, init_end_spec_vec)
    testit::assert(sum(length(which(num_native_spec == 1)),
                       length(which(num_native_spec == 2)),
                       length(which(num_native_spec == 3)))
                   == sum(mainland_n))
    init_nonend_spec <- length(init_nonend_spec_vec)
    init_end_spec <- length(init_end_spec_vec)
    if (length(mainland_spec) == 0) {
      mainland_spec <- 0
    }
  } else {
    init_nonend_spec <- 0
    init_end_spec <- 0
    init_nonend_spec_vec <- integer(0)
    init_end_spec_vec <- integer(0)
    
    if(mainland_n != 0){
      mainland_spec <- seq(1, mainland_n, 1)
    }else{
      mainland_spec = c()
    }
    # mainland_spec <- seq(1, mainland_n, 1)
  }
  return(list(init_nonend_spec = init_nonend_spec,
              init_end_spec = init_end_spec,
              init_nonend_spec_vec = init_nonend_spec_vec,
              init_end_spec_vec = init_end_spec_vec,
              mainland_spec = mainland_spec))
}

#' Update internal Gillespie bookeeping objects
#'
#' @param stt_table A species=through-time table.
#' @param totaltime Simulated amount of time.
#' @param timeval Current time of simulation.
#' @param mainland_spec A vector with the numeric IDs of the mainland species
#' (i.e. potential colonizers).
#' @param island_spec A matrix with the species on the island (state of the
#' system at each time point).
#'
#' @return A named list with the updated input arguments at time of simulation.
#'
#' @noRd
DAISIE_spec_tables <- function(stt_table,
                               totaltime,
                               timeval,
                               nonoceanic_sample,
                               island_spec) {
  init_nonend_spec <- nonoceanic_sample$init_nonend_spec
  init_end_spec <- nonoceanic_sample$init_end_spec
  mainland_spec <- nonoceanic_sample$mainland_spec
  stt_table[1, ] <- c(totaltime,
                      init_nonend_spec,
                      init_end_spec,
                      0)
  if (init_nonend_spec != 0) {
    for (i in seq_along(1:init_nonend_spec)) {
      island_spec <- rbind(island_spec,
                           c(nonoceanic_sample$init_nonend_spec_vec[i],
                             nonoceanic_sample$init_nonend_spec_vec[i],
                             timeval,
                             "I",
                             NA,
                             NA,
                             NA))
    }
  }
  if (init_end_spec != 0) {
    for (j in seq_along(1:init_end_spec)) {
      island_spec <- rbind(island_spec,
                           c(nonoceanic_sample$init_end_spec_vec[j],
                             nonoceanic_sample$init_end_spec_vec[j],
                             timeval,
                             "A",
                             NA,
                             NA,
                             NA))
    }
  }
  return(list(stt_table = stt_table,
              init_nonend_spec = init_nonend_spec,
              init_end_spec = init_end_spec,
              mainland_spec = mainland_spec,
              island_spec = island_spec))
}

#' Create an empty phylogeny
#' @param age Time of phylogeny
#' @author Giovanni Laudanno
#' @export
create_singleton_phylo <- function(age) {
  tr <- list(edge = matrix(c(2, 1), 1, 2), tip.label = "t1", Nnode = 1L)
  class(tr) <- "phylo"
  tr$edge.length <- age
  tr$tip.label <- "stem"
  tr
}

#' #' Unsampled CS full STT
#' #'
#' #' @param stt_list List of full stt tables as
#' #' returned by DAISIE_sim_core functions
#' #' @param totaltime Numeric double with total time of simulation.
#' #' @param stac_vec Vector with status of species on island.
#' #' @param trait_pars A named list containing diversification rates considering 
#' #' two trait states created by \code{\link{create_trait_pars}}:
#' #' \itemize{
#' #'   \item{[1]:A numeric with the per capita transition rate with state1}
#' #'   \item{[2]:A numeric with the per capita immigration rate with state2}
#' #'   \item{[3]:A numeric with the per capita extinction rate with state2}
#' #'   \item{[4]:A numeric with the per capita anagenesis rate with state2}
#' #'   \item{[5]:A numeric with the per capita cladogenesis rate with state2}
#' #'   \item{[6]:A numeric with the per capita transition rate with state2} 
#' #'   \item{[7]:A numeric with the number of species with trait state 2 on mainland} 
#' #' }
#' #'
#' #' @return 1 complete, unsampled STT table from all clades in an island of a
#' #' CS model as generated by DAISIE_sim_core functions.
#' #' @author Pedro Neves, Joshua Lambert, Shu Xie, Giovanni Laudanno
#' create_full_CS_stt <- function(stt_list, stac_vec, totaltime, trait_pars = NULL) {
#'   
#'   if(!is.null(trait_pars)){
#'     return(
#'       create_full_CS_stt_trait(
#'         stt_list = stt_list,
#'         stac_vec = stac_vec,
#'         totaltime = totaltime,
#'         trait_pars = trait_pars
#'       )
#'     )
#'   }
#'   # Return empty island, if empty
#'   present <- which(stac_vec != 0)
#' 
#'   # Checks if stt has only 2 rows and is empty at present (nothing happened)
#'   second_line_stts <- lapply(stt_list, "[", 2,)
#'   zeros_second_line <- sapply(second_line_stts, sum) == 0
#' 
#' 
#'   filled_stt_lists <- stt_list[!zeros_second_line]
#'   # If no colonization ever happened, just return 0 values
#'   if (length(filled_stt_lists) == 0) {
#'     times <- c(totaltime, 0)
#'     nI <- c(0, 0)
#'     nA <- c(0, 0)
#'     nC <- c(0, 0)
#'     diff_present <- c(0, 0)
#'   } else {
#' 
#'     deltas_matrix <- lapply(filled_stt_lists, FUN = diff)
#'     for (i in seq_along(filled_stt_lists)) {
#'       if (any(filled_stt_lists[[i]][1, ] !=
#'               c("Time" = totaltime, "nI" = 0, "nA" = 0, "nC" = 0))) {
#'         deltas_matrix[[i]] <- rbind(
#'           filled_stt_lists[[i]][1, ],
#'           deltas_matrix[[i]]
#'         )
#'       } else {
#'         deltas_matrix[[i]] <- rbind(
#'           c("Time" = totaltime, "nI" = 0, "nA" = 0, "nC" = 0),
#'           deltas_matrix[[i]]
#'         )
#'       }
#'     }
#' 
#'     times_list <- lapply(filled_stt_lists, "[", , 1) # nolint
#'     all_times <- unlist(times_list)
#'     times <- all_times
#' 
#'     nI_list <- lapply(deltas_matrix, "[", , 2) # nolint
#'     nA_list <- lapply(deltas_matrix, "[", , 3) # nolint
#'     nC_list <- lapply(deltas_matrix, "[", , 4) # nolint
#' 
#'     nI <- unlist(nI_list)
#'     nA <- unlist(nA_list)
#'     nC <- unlist(nC_list)
#'     diff_present <- nI + nA + nC
#'   }
#' 
#'   full_stt <- data.frame(
#'     times = times,
#'     nI = nI,
#'     nA = nA,
#'     nC = nC,
#'     present = diff_present
#'   )
#'   ordered_diffs <- full_stt[order(full_stt$times, decreasing = TRUE), ]
#' 
#'   complete_stt_table <- mapply(ordered_diffs[2:5], FUN = cumsum)
#'   complete_stt_table <- cbind(ordered_diffs$times, complete_stt_table)
#'   colnames(complete_stt_table) <- c("Time", "nI", "nA", "nC", "present")
#' 
#'   while (complete_stt_table[1, 1] == complete_stt_table[2, 1]) {
#'     complete_stt_table <- complete_stt_table[-1, ]
#'   }
#' 
#'   stt <- complete_stt_table
#'   # Remove final duplicate lines, if any
#'   while (
#'     all(stt[nrow(stt) - 1, ] == stt[nrow(stt), ])
#'   ) {
#'     stt <- stt[1:(nrow(stt) - 1), ]
#'   }
#'   stt
#' }

# create_full_CS_stt_trait <- function(stt_list, stac_vec, totaltime, trait_pars) {
#   # Return empty island, if empty
#   present <- which(stac_vec != 0)
#   
#   # Checks if stt has only 2 rows and is empty at present (nothing happened)
#   second_line_stts <- lapply(stt_list, "[", 2,)
#   zeros_second_line <- sapply(second_line_stts, sum) == 0
#   
#   
#   filled_stt_lists <- stt_list[!zeros_second_line]
#   # If no colonization ever happened, just return 0 values
#   if (length(filled_stt_lists) == 0) {
#     times <- c(totaltime, 0)
#     nI <- c(0, 0)
#     nA <- c(0, 0)
#     nC <- c(0, 0)
#     nI2 <- c(0, 0)
#     nA2 <- c(0, 0)
#     nC2 <- c(0, 0)
#     diff_present <- c(0, 0)
#   } else {
#     
#     deltas_matrix <- lapply(filled_stt_lists, FUN = diff)
#     for (i in seq_along(filled_stt_lists)) {
#       if (any(filled_stt_lists[[i]][1, ] !=
#               c("Time" = totaltime, "nI" = 0, "nA" = 0, "nC" = 0, "nI2" = 0, "nA2" = 0, "nC2" = 0))) {
#         deltas_matrix[[i]] <- rbind(
#           filled_stt_lists[[i]][1, ],
#           deltas_matrix[[i]]
#         )
#       } else {
#         deltas_matrix[[i]] <- rbind(
#           c("Time" = totaltime, "nI" = 0, "nA" = 0, "nC" = 0, "nI2" = 0, "nA2" = 0, "nC2" = 0),
#           deltas_matrix[[i]]
#         )
#       }
#     }
#     
#     times_list <- lapply(filled_stt_lists, "[", , 1) # nolint
#     all_times <- unlist(times_list)
#     times <- all_times
#     
#     nI_list <- lapply(deltas_matrix, "[", , 2) # nolint
#     nA_list <- lapply(deltas_matrix, "[", , 3) # nolint
#     nC_list <- lapply(deltas_matrix, "[", , 4) # nolint
#     nI2_list <- lapply(deltas_matrix, "[", , 5) # nolint
#     nA2_list <- lapply(deltas_matrix, "[", , 6) # nolint
#     nC2_list <- lapply(deltas_matrix, "[", , 7) # nolint
#     
#     nI <- unlist(nI_list)
#     nA <- unlist(nA_list)
#     nC <- unlist(nC_list)
#     nI2 <- unlist(nI2_list)
#     nA2 <- unlist(nA2_list)
#     nC2 <- unlist(nC2_list)
#     diff_present <- nI + nA + nC + nI2 + nA2 + nC2
#   }
#   
#   full_stt <- data.frame(
#     times = times,
#     nI = nI,
#     nA = nA,
#     nC = nC,
#     nI2 = nI2,
#     nA2 = nA2,
#     nC2 = nC2,
#     present = diff_present
#   )
#   ordered_diffs <- full_stt[order(full_stt$times, decreasing = TRUE), ]
#   
#   complete_stt_table <- mapply(ordered_diffs[2:8], FUN = cumsum)
#   complete_stt_table <- cbind(ordered_diffs$times, complete_stt_table)
#   colnames(complete_stt_table) <- c("Time", "nI", "nA", "nC", "nI2", "nA2", "nC2", "present")
#   
#   while (complete_stt_table[1, 1] == complete_stt_table[2, 1]) {
#     complete_stt_table <- complete_stt_table[-1, ]
#   }
#   
#   stt <- complete_stt_table
#   # Remove final duplicate lines, if any
#   while (
#     all(stt[nrow(stt) - 1, ] == stt[nrow(stt), ])
#   ) {
#     stt <- stt[1:(nrow(stt) - 1), ]
#   }
#   stt
# }
