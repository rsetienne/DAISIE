#' Counts the number of species
#'
#' @param datalistelement stub
#'
#' @return A numeric
#' @export
countspecies <- function(datalistelement) {
  N <- length(datalistelement$branching_times) -
    1 + datalistelement$missing_species
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

#' Checks whether an input is odd
#'
#' @param x Object to determine
#'
#' @return A boolean indicating if object is odd
#' @examples DAISIE:::is_odd(5)
is_odd <- function(x) {
  if (!is.numeric(x) || length(x) > 1) {
    stop("'x' must be a single numeric")
  }
  if (!x %% 1 == 0) {
    stop("'x' must be an integer")
  }
  res <- x %% 2
  if (res != 0) {
    out <- TRUE
  } else {
    out <- FALSE
  }
  return(out)
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

quantiles <- function(probdist, probs) {
  result <- NULL
  cdf <- cumsum(probdist[2, ])
  for (i in seq_along(probs)) {
    n <- max(which(cdf <= probs[i]))
    x <- probdist[1, n]
    if (cdf[n] == probs[i]) {
      result[i] <- x
    } else
      if (n < length(cdf)) {
        result[i] <- ((x + 1) * (probs[i] - cdf[n]) + x * (cdf[n + 1] - probs[i])) / (cdf[n + 1] - cdf[n])
      } else {
        result[i] <- x
      }
  }
  names(result) <- probs
  return(result)
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
#' @inherit DAISIE_sim
#'
#' @return Numeric, 0 for null-ontogeny, 1 for beta function
#' @export
#' @examples translate_island_ontogeny("const")
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
#' @inherit DAISIE_sim
#'
#' @return Numeric, 0 for null-sea-level, 1 for sine function
#' @export
#' @examples translate_sea_level("const")
translate_sea_level <- function(sea_level) {

  if (sea_level == "const" || sea_level == 0) {
    sea_level <- 0
  }

  if (sea_level == "sine" || sea_level == 1) {
    sea_level <- 1
  }
  return(sea_level)
}

order_pars1 <- function(pars1) {
  np <- names(pars1)
  correct_order <- c("max_area", "proportional_peak_t",
                     "peak_sharpness", "total_island_age",
                     "lac", "mu_min", "mu_max", "K0", "gam", "laa")
  if (!is.null(np)) {
    pars1ff <- pars1
    pars1ff[1] <- pars1[which(names(pars1) == "max_area")]
    pars1ff[2] <- pars1[which(names(pars1) == "proportional_peak_t")]
    pars1ff[3] <- pars1[which(names(pars1) == "peak_sharpness")]
    pars1ff[4] <- pars1[which(names(pars1) == "total_island_age")]
    pars1ff[5] <- pars1[which(names(pars1) == "lac")]
    pars1ff[6] <- pars1[which(names(pars1) == "mu_min")]
    pars1ff[7] <- pars1[which(names(pars1) == "mu_max")]
    pars1ff[8] <- pars1[which(names(pars1) == "K0")]
    pars1ff[9] <- pars1[which(names(pars1) == "gam")]
    pars1ff[10] <- pars1[which(names(pars1) == "laa")]
    pars1 <- pars1ff
    names(pars1) <- correct_order
  }
  return(pars1)
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

#' Calculates the species on the island initially when \code{island_type =
#' 'noncoeanic'}
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
  return(list(init_nonend_spec_vec, init_end_spec_vec, mainland_spec))
}

DAISIE_nonoceanic_stt_table <- function(stt_table,
                                        totaltime,
                                        timeval,
                                        init_nonend_spec_vec,
                                        init_end_spec_vec,
                                        mainland_spec,
                                        island_spec) {
  stt_table[1, ] <- c(totaltime,
                      length(init_nonend_spec_vec),
                      length(init_end_spec_vec),
                      0)
  if (length(init_nonend_spec_vec) == 0) {
    init_nonend_spec <- 0
  } else {
    init_nonend_spec <- length(init_nonend_spec_vec)
  }
  if (length(init_end_spec_vec) == 0) {
    init_end_spec <- 0
  } else {
    init_end_spec <- length(init_end_spec_vec)
  }
  if (length(mainland_spec) == 0) {
    mainland_spec <- 0
  }
  if (length(init_nonend_spec_vec) == 1 &&
      init_nonend_spec_vec != 0 || length(init_nonend_spec_vec) > 1) {
    for (i in seq_along(init_nonend_spec_vec)) {
      island_spec <- rbind(island_spec,
                           c(init_nonend_spec_vec[i],
                             init_nonend_spec_vec[i],
                             timeval,
                             "I",
                             NA,
                             NA,
                             NA))
    }
  }
  if (length(init_end_spec_vec) == 1 &&
      init_end_spec_vec != 0 || length(init_end_spec_vec) > 1) {
    for (j in seq_along(init_end_spec_vec)) {
      island_spec <- rbind(island_spec,
                           c(init_end_spec_vec[j],
                             init_end_spec_vec[j],
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
#' Create a full-blown DAISIE parameter structure
#' @param time something
#' @param M something
#' @param pars something
#' @param replicates something
#' @export
create_daisie_pars <- function(time, M, pars, replicates) {
  # testit::assert(time > 0)
  if (length(M) > 1) {
    stop("'M' must be one non-zero and positive value")
  }
  if (length(time) > 1) {
    stop("'time' must be one non-zero and positive value")
  }
  if (length(pars) < 5) {
    stop("'pars' must have a length of at least 5")
  }
  if (time <= 0) {
    stop("'time' must be non-zero and positive")
  }
  if (M <= 0) {
    stop("'M' must be non-zero and positive")
  }
  if (replicates <= 0) {
    stop("'replicates' must be non-zero and positive")
  }
  if (pars[1] < 0 || pars[2] < 0 || pars[3] < 0 || pars[4] < 0 || pars[5] < 0) {
    stop("'pars' must be non-zero and positive")
  }
  list(time = time,
       M = M,
       pars = pars,
       replicates = replicates
  )
}
#' Create a sunction to test full-blown DAISIE parameter structure
#' @export
create_test_daisie_pars <- function() {
  create_daisie_pars(time = 3,
                     M = 1,
                     pars = c(2.5, 2.6, Inf, 0.01, 1.0),
                     replicates = 1)

}

#' Determines which rate set to use in the shift-rates simulation
#'
#' @param timeval numeric current time during a simulation.
#' @param totaltime simulated amount of time
#' @param shift_times a vector of numerics with times of shifts,
#' times are time back from the present.
#'
#' @return a list of Booleans and a numeric
#'
#' @examples periods <- DAISIE:::land_bridge_periods(timeval = 0.5,
#'                               totaltime = 10,
#'                               shift_times = c(3, 6))
land_bridge_periods <- function(timeval,
                                totaltime,
                                shift_times) {
  testit::assert(is.numeric(timeval))
  testit::assert(is.numeric(totaltime))
  testit::assert(is.numeric(shift_times) || is.null(shift_times))

  if (is.null(shift_times)) {
    return(list(
      present = FALSE,
      shift_time = "no_shift"
    ))
  }

  testit::assert(totaltime >= max(shift_times))
  shift_times <- totaltime - shift_times
  shift_times <- sort(shift_times)
  list_length <- length(shift_times) %/% 2 + length(shift_times) %% 2
  testit::assert(is.numeric(list_length) && length(list_length) > 0)
  if (length(shift_times) == 1) {
    land_bridge_periods <- list(c(shift_times, totaltime))
    island_periods <- list(c(0, shift_times))
  } else if (length(shift_times) == 2) {
    land_bridge_periods <- list(c(shift_times))
    island_periods <- list(c(0, shift_times[1]), c(shift_times[2], totaltime))
  } else if (is_odd(length(shift_times))) {
    land_bridge_periods <- unname(split(
      shift_times,
      as.numeric(gl(length(shift_times), 2, length(shift_times)))
    ))
    land_bridge_periods[[length(land_bridge_periods)]][2] <- totaltime
    island_periods <- c(0, shift_times)
    island_periods <- unname(split(
      island_periods,
      as.numeric(gl(length(island_periods), 2, length(island_periods)))
    ))
  } else {
    land_bridge_periods <- unname(split(
      shift_times,
      as.numeric(gl(length(shift_times), 2, length(shift_times)))
    ))
    island_periods <- c(0, shift_times)
    island_periods <- unname(split(
      island_periods,
      as.numeric(gl(length(island_periods), 2, length(island_periods)))
    ))
    island_periods[[length(island_periods)]][2] <- totaltime
  }
  testit::assert(is.list(land_bridge_periods))
  testit::assert(is.list(island_periods))
  land_bridge_eval <- c()
  island_eval <- c()
  for (i in 1:length(land_bridge_periods)) {
    if (timeval >= land_bridge_periods[[i]][1] &
        timeval < land_bridge_periods[[i]][2]) {
      land_bridge_eval[i] <- TRUE
    } else {
      land_bridge_eval[i] <- FALSE
    }
  }
  for (i in 1:length(island_periods)) {
    if (timeval >= island_periods[[i]][1] &
        timeval < island_periods[[i]][2]) {
      island_eval[i] <- TRUE
    } else {
      island_eval[i] <- FALSE
    }
  }
  if (any(land_bridge_eval) == TRUE) {
    return(list(
      present = TRUE,
      shift_time = land_bridge_periods[which(land_bridge_eval == TRUE)][[1]][1]
    ))
  } else {
    return(list(
      present = FALSE,
      shift_time = island_periods[which(island_eval == TRUE)][[1]][1]
    ))
  }
}

#' Unsampled CS full STT
#'
#' @param stt_list List of full stt tables as
#' returned by \code{\link{DAISIE_sim_core}}
#' @param totaltime Numeric double with total time of simulation.
#' @param stac_vec Vector with status of species on island.
#'
#' @return 1 Complete, unsampled STT table from all clades in an island of a
#' CS model as generated by \code{\link{DAISIE_sim_core}}.
#' @author Pedro Neves, Joshua Lambert, Shu Xie, Giovanni Laudanno

create_full_CS_stt <- function(stt_list, stac_vec, totaltime) {
  # Return empty island, if empty
  present <- which(stac_vec != 0)
  number_present <- length(present)
  if (number_present == 0) {
    stt <- matrix(c(totaltime, rep(0, 9)), ncol = 5)
    colnames(stt) <- c("Time", "nI", "nA", "nC", "present")
  } else {

    small_stts <- lapply(stt_list, nrow) == 2
    second_line_stts <- lapply(stt_list, "[", 2,)
    zeros_second_line <- sapply(second_line_stts, sum) == 0

    comparison <- zeros_second_line == small_stts
    testit::assert(all(comparison))

    filled_stt_lists <- stt_list[!zeros_second_line]


    deltas_matrix <- lapply(filled_stt_lists, FUN = diff)
    for (i in seq_along(filled_stt_lists)) {
      if (any(filled_stt_lists[[i]][1, ] !=
              c("Time" = totaltime, "nI" = 0, "nA" = 0, "nC" = 0))) {
        deltas_matrix[[i]] <- rbind(
          filled_stt_lists[[i]][1, ],
          deltas_matrix[[i]]
        )
      } else {
        deltas_matrix[[i]] <- rbind(
          c("Time" = totaltime, "nI" = 0, "nA" = 0, "nC" = 0),
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

    nI <- unlist(nI_list)
    nA <- unlist(nA_list)
    nC <- unlist(nC_list)
    diff_present <- nI + nA + nC

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
  }
  stt
}

#' Sets standard metaparameters to defaults when NULL
#'
#' @inheritParams DAISIE_sim_core
#'
#' @return List with standard metaparameters
set_default_pars <- function(island_ontogeny,
                             sea_level,
                             hyper_pars,
                             dist_pars,
                             ext_pars,
                             totaltime,
                             pars) {
  if (island_ontogeny == 0 && sea_level == 0) {
    area_pars <- create_area_pars(
      max_area = 1,
      proportional_peak_t = 0,
      peak_sharpness = 0,
      total_island_age = totaltime,
      sea_level_amplitude = 0,
      sea_level_frequency = 0
    )
  }
  if (is.null(hyper_pars)) {
    hyper_pars <- create_hyper_pars(d_0 = 0, x = 0, alpha = 0, beta = 0)
  }
  if (is.null(dist_pars)) {
    dist_pars <- create_dist_pars(D = 1)
  }
  if (is.null(ext_pars)) {
    ext_pars <- pars[2]
  }
  out <- list(
    area_pars = area_pars,
    hyper_pars = hyper_pars,
    dist_pars = dist_pars,
    ext_pars = ext_pars
  )
  return(out)
}
