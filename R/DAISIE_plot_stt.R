#' Prepare input for DAISIE_stt
#'
#' @inheritParams default_params_doc
#'
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
DAISIE_convert_to_classic_plot <- function(simulation_outputs, trait_pars = NULL) {
  if (!DAISIE::is_simulation_outputs(simulation_outputs)) {
    stop(
      "'simulation_outputs' should be a set of simulation outputs. \n",
      "Actual value: ", simulation_outputs
    )
  }
  replicates <- length(simulation_outputs)
  ### STT ALL species
  s_freq <- length(simulation_outputs[[1]][[1]]$stt_all[, 1])
  complete_arr <- array(dim = c(s_freq, 6, replicates))
  for (x in 1:replicates) {
    if(is.null(trait_pars)){
      sum_endemics <- simulation_outputs[[x]][[1]]$stt_all[, "nA"] +
        simulation_outputs[[x]][[1]]$stt_all[, "nC"]
      total <- simulation_outputs[[x]][[1]]$stt_all[, "nA"] +
        simulation_outputs[[x]][[1]]$stt_all[, "nC"] +
        simulation_outputs[[x]][[1]]$stt_all[, "nI"]
      complete_arr[, , x] <- cbind(simulation_outputs[[x]][[1]]$stt_all[, c("Time", "nI", "nA", "nC")],
                                   sum_endemics,
                                   total)
    }else{
      sum_endemics <- simulation_outputs[[x]][[1]]$stt_all[, "nA"] +
        simulation_outputs[[x]][[1]]$stt_all[, "nC"] +
        simulation_outputs[[x]][[1]]$stt_all[, "nA2"] +
        simulation_outputs[[x]][[1]]$stt_all[, "nC2"]
      total <- simulation_outputs[[x]][[1]]$stt_all[, "nA"] +
        simulation_outputs[[x]][[1]]$stt_all[, "nC"] +
        simulation_outputs[[x]][[1]]$stt_all[, "nI"] +
        simulation_outputs[[x]][[1]]$stt_all[, "nA2"] +
        simulation_outputs[[x]][[1]]$stt_all[, "nC2"] +
        simulation_outputs[[x]][[1]]$stt_all[, "nI2"]
      nI <- simulation_outputs[[x]][[1]]$stt_all[, "nI"] +
        simulation_outputs[[x]][[1]]$stt_all[, "nI2"]
      nA <- simulation_outputs[[x]][[1]]$stt_all[, "nA"] +
        simulation_outputs[[x]][[1]]$stt_all[, "nA2"]
      nC <- simulation_outputs[[x]][[1]]$stt_all[, "nC"] +
        simulation_outputs[[x]][[1]]$stt_all[, "nC2"]
      complete_arr[,,x]<-cbind(simulation_outputs[[x]][[1]]$stt_all[, 'Time'],
                               nI,
                               nA,
                               nC,
                               sum_endemics,
                               total)
    }
    
  }
  stt_average_all <- apply(complete_arr, c(1, 2), stats::median)
  testit::assert(stt_average_all ==
                   DAISIE::DAISIE_extract_stt_median(island_replicates = simulation_outputs,trait_pars = trait_pars))
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

#' Create the Species-Through-Time plot. This is used to visualize
#' the output of DAISIE_sim functions
#'
#' @inheritParams default_params_doc
#'
#' @export
DAISIE_plot_stt <- function(
  plot_plus_one = TRUE,
  time,
  plot_lists = plot_lists,
  type = type
) {
  # Plot the y axis iff plus one
  y_axis_type <- "n"
  y_axis_label <- "No of species"
  if (plot_plus_one == TRUE) {
    y_axis_type <- "s"
    y_axis_label <- "No of species + 1"
  }
  stt <- plot_lists[[type]]
  if (is.null(stt)) {
    return()
  }
  suppressWarnings(
    graphics::plot(
      NULL, NULL, xlim = rev(c(0, time)), ylim = c(1, max(stt$stt_q0.975)),
      ylab = y_axis_label,
      bty = "l", xaxs = "i", xlab = "Time before present",
      main = "Species-through-time - All species",
      log = "y", cex.lab = 1.2, cex.main = 1.2, cex.axis = 1.2,
      yaxt = y_axis_type
    )
  )
  graphics::polygon(c(stt$stt_average[, "Time"], rev(stt$stt_average[, "Time"])), c(stt$stt_q0.025[, "Total"] +
                                                                              1, rev(stt$stt_q0.975[, "Total"] + 1)), col = "light grey", border = NA)
  graphics::polygon(c(stt$stt_average[, "Time"], rev(stt$stt_average[, "Time"])), c(stt$stt_q0.25[, "Total"] +
                                                                              1, rev(stt$stt_q0.75[, "Total"] + 1)), col = "dark grey", border = NA)
  graphics::lines(stt$stt_average[, "Time"], stt$stt_average[, "Total"] + 1, lwd = 2)
  graphics::lines(stt$stt_average[, "Time"], stt$stt_average[, "nI"] + 1, lwd = 2, col = "cyan3")
  graphics::lines(stt$stt_average[, "Time"], stt$stt_average[, "Endemic"] + 1, lwd = 2, col = "dodgerblue1")
  legend_names <- c("Total", "Non-endemic", "Endemic")
  legend_colors <- c("black", "cyan3", "dodgerblue1")
  graphics::legend(
    time, max(stt$stt_q0.975), legend_names, lty = 1, lwd = 2,
    col = legend_colors, cex = 1.2, border = NA, bty = "n"
  )
  if (plot_plus_one == FALSE) {
    y_axis_values <- c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
    graphics::axis(2, at = y_axis_values, labels = y_axis_values - 1)
  }
}
