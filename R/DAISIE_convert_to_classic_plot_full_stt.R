#' Prepare input for DAISIE_stt
#'
#' @param simulation_outputs A list with matrices? of simulation produced by
#' DAISIE_sim.
#' @seealso \code{\link{DAISIE_plot_stt}}, \code{\link{DAISIE_plot_sims}}
#' @examples
#' utils::data("islands_1type_1000reps", package = "DAISIE")
#' simulation_outuputs <- DAISIE::DAISIE_convert_to_classic_plot(
#' islands_1type_1000reps
#' )
#'
#'
#' @return a list with wrangled data to be used for plotting STT plots with
#' DAISIE_plot_stt
#' @export
DAISIE_convert_to_classic_plot_full_stt <- function(simulation_outputs) {
  if (!DAISIE::is_simulation_outputs(simulation_outputs)) {
    stop(
      "'simulation_outputs' should be a set of simulation outputs. \n",
      "Actual value: ", simulation_outputs
    )
  }
  replicates <- length(simulation_outputs)
  ### STT ALL species
  extract_stts <- lapply(X = simulation_outputs, FUN = mylist mylist[[1]]$stt_all))
  s_freq <- lapply(X = extract_stts, FUN = nrow())

  complete_arr <- array(dim = c(s_freq, 6, replicates))
  for (x in 1:replicates) {
    sum_endemics <- simulation_outputs[[x]][[1]]$stt_all[, "nA"] +
      simulation_outputs[[x]][[1]]$stt_all[, "nC"]
    total <- simulation_outputs[[x]][[1]]$stt_all[, "nA"] +
      simulation_outputs[[x]][[1]]$stt_all[, "nC"] +
      simulation_outputs[[x]][[1]]$stt_all[, "nI"]
    complete_arr[, , x] <- cbind(simulation_outputs[[x]][[1]]$stt_all[, c("Time", "nI", "nA", "nC")],
                                 sum_endemics,
                                 total)
  }
  stt_average_all <- apply(complete_arr, c(1, 2), stats::median)
  testit::assert(stt_average_all ==
                   DAISIE::DAISIE_extract_stt_median(simulation_outputs))
  stt_q0.025_all <- apply(complete_arr, c(1, 2), stats::quantile, 0.025)
  stt_q0.25_all <- apply(complete_arr, c(1, 2), stats::quantile, 0.25)
  stt_q0.75_all <- apply(complete_arr, c(1, 2), stats::quantile, 0.75)
  stt_q0.975_all <- apply(complete_arr, c(1, 2), stats::quantile, 0.975)
  colnames(stt_average_all) <- c("Time", "nI", "nA", "nC", "Endemic", "Total")
  colnames(stt_q0.025_all) <- c("Time", "nI", "nA", "nC", "Endemic", "Total")
  colnames(stt_q0.25_all) <- c("Time", "nI", "nA", "nC", "Endemic", "Total")
  colnames(stt_q0.75_all) <- c("Time", "nI", "nA", "nC", "Endemic", "Total")
  colnames(stt_q0.975_all) <- c("Time", "nI", "nA", "nC", "Endemic", "Total")
  all_species <- list(
    stt_average = stt_average_all,
    stt_q0.025 = stt_q0.025_all,
    stt_q0.25 = stt_q0.25_all,
    stt_q0.75 = stt_q0.75_all,
    stt_q0.975 = stt_q0.975_all
  )
  if (is.null(simulation_outputs[[1]][[1]]$stt_type1) == FALSE) {
    ### STT TYPE1
    s_freq <- length(simulation_outputs[[1]][[1]]$stt_type1[, 1])
    complete_arr <- array(dim = c(s_freq, 7, replicates))
    for (x in 1:replicates) {
      sum_endemics <- simulation_outputs[[x]][[1]]$stt_type1[, "nA"] +
        simulation_outputs[[x]][[1]]$stt_type1[, "nC"]
      total <- simulation_outputs[[x]][[1]]$stt_type1[, "nA"] +
        simulation_outputs[[x]][[1]]$stt_type1[, "nC"] +
        simulation_outputs[[x]][[1]]$stt_type1[, "nI"]
      complete_arr[, , x] <- cbind(simulation_outputs[[x]][[1]]$stt_type1,
                                   sum_endemics,
                                   total)
    }
    stt_average_type1 <- apply(complete_arr, c(1, 2), stats::median)
    stt_q0.025_type1 <- apply(complete_arr, c(1, 2), stats::quantile, 0.025)
    stt_q0.25_type1 <- apply(complete_arr, c(1, 2), stats::quantile, 0.25)
    stt_q0.75_type1 <- apply(complete_arr, c(1, 2), stats::quantile, 0.75)
    stt_q0.975_type1 <- apply(complete_arr, c(1, 2), stats::quantile, 0.975)
    colnames(stt_average_type1) <- c(
      "Time",
      "nI",
      "nA",
      "nC",
      "present",
      "Endemic",
      "Total"
    )
    colnames(stt_q0.025_type1) <- c(
      "Time",
      "nI",
      "nA",
      "nC",
      "present",
      "Endemic",
      "Total"
    )
    colnames(stt_q0.25_type1) <- c(
      "Time",
      "nI",
      "nA",
      "nC",
      "present",
      "Endemic",
      "Total"
    )
    colnames(stt_q0.75_type1) <- c(
      "Time",
      "nI",
      "nA",
      "nC",
      "present",
      "Endemic",
      "Total"
    )
    colnames(stt_q0.975_type1) <- c(
      "Time",
      "nI",
      "nA",
      "nC",
      "present",
      "Endemic",
      "Total"
    )
    type1_species <- list(
      stt_average = stt_average_type1,
      stt_q0.025 = stt_q0.025_type1,
      stt_q0.25 = stt_q0.25_type1,
      stt_q0.75 = stt_q0.75_type1,
      stt_q0.975 = stt_q0.975_type1
    )
    ### STT TYPE2
    s_freq <- length(simulation_outputs[[1]][[1]]$stt_type2[, 1])
    complete_arr <- array(dim = c(s_freq, 7, replicates))
    for (x in 1:replicates) {
      sum_endemics <- simulation_outputs[[x]][[1]]$stt_type2[, "nA"] +
        simulation_outputs[[x]][[1]]$stt_type2[, "nC"]
      total <- simulation_outputs[[x]][[1]]$stt_type2[, "nA"] +
        simulation_outputs[[x]][[1]]$stt_type2[, "nC"] +
        simulation_outputs[[x]][[1]]$stt_type2[, "nI"]
      complete_arr[, , x] <- cbind(
        simulation_outputs[[x]][[1]]$stt_type2,
        sum_endemics,
        total
      )
    }
    stt_average_type2 <- apply(complete_arr, c(1, 2), stats::median)
    stt_q0.025_type2 <- apply(complete_arr, c(1, 2), stats::quantile, 0.025)
    stt_q0.25_type2 <- apply(complete_arr, c(1, 2), stats::quantile, 0.25)
    stt_q0.75_type2 <- apply(complete_arr, c(1, 2), stats::quantile, 0.75)
    stt_q0.975_type2 <- apply(complete_arr, c(1, 2), stats::quantile, 0.975)
    colnames(stt_average_type2) <- c(
      "Time",
      "nI",
      "nA",
      "nC",
      "present",
      "Endemic",
      "Total"
    )
    colnames(stt_q0.025_type2) <- c(
      "Time",
      "nI",
      "nA",
      "nC",
      "present",
      "Endemic",
      "Total"
    )
    colnames(stt_q0.25_type2) <- c(
      "Time",
      "nI",
      "nA",
      "nC",
      "present",
      "Endemic",
      "Total"
    )
    colnames(stt_q0.75_type2) <- c(
      "Time",
      "nI",
      "nA",
      "nC",
      "present",
      "Endemic",
      "Total"
    )
    colnames(stt_q0.975_type2) <- c(
      "Time",
      "nI",
      "nA",
      "nC",
      "present",
      "Endemic",
      "Total"
    )
    type2_species <- list(
      stt_average = stt_average_type2,
      stt_q0.025 = stt_q0.025_type2,
      stt_q0.25 = stt_q0.25_type2,
      stt_q0.75 = stt_q0.75_type2,
      stt_q0.975 = stt_q0.975_type2
    )
    return(list(
      all_species = all_species,
      type1_species = type1_species,
      type2_species = type2_species
    )
    )
  } else {
    return(list(
      all_species = all_species,
      type1_species = NULL,
      type2_species = NULL)
    )
  }
}
