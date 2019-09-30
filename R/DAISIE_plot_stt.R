#' Prepare input for DAISIE_stt
#'
#' @inheritParams DAISIE_plot_sims
#' @param simulation_outputs A list with matrices? of simulation produced by
#' DAISIE_sim.
#' @param Tpars A named list containing diversification rates considering two trait states:
#' \itemize{
#'   \item{[1]:A numeric with the per capita transition rate with state1}
#'   \item{[2]:A numeric with the per capita immigration rate with state2}
#'   \item{[3]:A numeric with the per capita extinction rate with state2}
#'   \item{[4]:A numeric with the per capita anagenesis rate with state2}
#'   \item{[5]:A numeric with the per capita cladogenesis rate with state2}
#'   \item{[6]:A numeric with the per capita transition rate with state2}
#'   \item{[7]:A numeric with the number of species with trait state 2 on mainland}
#' }
#' @seealso \code{\link{DAISIE_plot_stt}}, \code{\link{DAISIE_plot_sims}}
#' @examples
#' utils::data("islands_1type_1000reps", package = "DAISIE")
#' simulation_outuputs <- DAISIE:::DAISIE_convert_to_classic_plot(
#' islands_1type_1000reps
#' )
#'
#'
#' @return a list with wrangled data to be used for plotting STT plots with
#' DAISIE_plot_stt
#' @export
DAISIE_convert_to_classic_plot <- function(simulation_outputs,
                                           Tpars = NULL) {
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
    testit::assert("nA" %in% colnames(simulation_outputs[[x]][[1]]$stt_all))
    testit::assert("nI" %in% colnames(simulation_outputs[[x]][[1]]$stt_all))
    testit::assert("nC" %in% colnames(simulation_outputs[[x]][[1]]$stt_all))
    testit::assert(x >= 1)
    testit::assert(x <= length(simulation_outputs))
    testit::assert(length(simulation_outputs[[x]]) >= 1)
    testit::assert("stt_all" %in% names(simulation_outputs[[x]][[1]]))
    testit::assert(is.matrix(simulation_outputs[[x]][[1]]$stt_all))


    if(is.null(Tpars)){
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
  # testit::assert(stt_average_all == DAISIE::DAISIE_extract_stt_median(simulation_outputs))
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
        "Total")
      colnames(stt_q0.025_type1) <- c(
        "Time",
        "nI",
        "nA",
        "nC",
        "present",
        "Endemic",
        "Total")
      colnames(stt_q0.25_type1) <- c(
        "Time",
        "nI",
        "nA",
        "nC",
        "present",
        "Endemic",
        "Total")
      colnames(stt_q0.75_type1) <- c(
        "Time",
        "nI",
        "nA",
        "nC",
        "present",
        "Endemic",
        "Total")
      colnames(stt_q0.975_type1) <- c(
        "Time",
        "nI",
        "nA",
        "nC",
        "present",
        "Endemic",
        "Total")
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
        type2_species = type2_species)
      )
    } else {
      return(list(
        all_species = all_species,
        type1_species = NULL,
        type2_species = NULL)
      )
    }
  # }else{
  #   stop("Don't consider two species types with two trait states")
  # }
  # }else{
  #
  #   if (is.null(simulation_outputs[[1]][[1]]$stt_type1) == FALSE){
  #     stop("Don't consider two species types with two trait states")
  #   }
  #   ### STT state1
  #   s_freq <- length(simulation_outputs[[1]][[1]]$stt_all[, 1])
  #   complete_arr <- array(dim = c(s_freq, 7, replicates))
  #
  #   for (x in 1:replicates) {
  #     sum_endemics <- simulation_outputs[[x]][[1]]$stt_all[, "nA"] +
  #       simulation_outputs[[x]][[1]]$stt_all[, "nC"]
  #     total <- simulation_outputs[[x]][[1]]$stt_all[, "nA"] +
  #       simulation_outputs[[x]][[1]]$stt_all[, "nC"] +
  #       simulation_outputs[[x]][[1]]$stt_all[, "nI"]
  #     nI <- simulation_outputs[[x]][[1]]$stt_all[,"nI"]
  #     nA <- simulation_outputs[[x]][[1]]$stt_all[,"nA"]
  #     nC <- simulation_outputs[[x]][[1]]$stt_all[,"nC"]
  #     complete_arr[, , x] <- cbind(simulation_outputs[[x]][[1]]$stt_all[,'Time'],nI,nA,nC,sum_endemics,total)
  #   }
  #
  #
  #   stt_average_state <- apply(complete_arr, c(1, 2), stats::median)
  #   stt_q0.025_state <- apply(complete_arr, c(1, 2), stats::quantile, 0.025)
  #   stt_q0.25_state <- apply(complete_arr, c(1, 2), stats::quantile, 0.25)
  #   stt_q0.75_state <- apply(complete_arr, c(1, 2), stats::quantile, 0.75)
  #   stt_q0.975_state <- apply(complete_arr, c(1, 2), stats::quantile, 0.975)
  #
  #   colnames(stt_average_state) <- c(
  #     "Time",
  #     "nI",
  #     "nA",
  #     "nC",
  #     "present",
  #     "Endemic",
  #     "Total"
  #   )
  #   colnames(stt_q0.025_state) <- c(
  #     "Time",
  #     "nI",
  #     "nA",
  #     "nC",
  #     "present",
  #     "Endemic",
  #     "Total"
  #   )
  #   colnames(stt_q0.25_state) <- c(
  #     "Time",
  #     "nI",
  #     "nA",
  #     "nC",
  #     "present",
  #     "Endemic",
  #     "Total"
  #   )
  #   colnames(stt_q0.75_state) <- c(
  #     "Time",
  #     "nI",
  #     "nA",
  #     "nC",
  #     "present",
  #     "Endemic",
  #     "Total"
  #   )
  #   colnames(stt_q0.975_state) <- c(
  #     "Time",
  #     "nI",
  #     "nA",
  #     "nC",
  #     "present",
  #     "Endemic",
  #     "Total"
  #   )
  #
  #   state0_species <- list(
  #     stt_average = stt_average_state,
  #     stt_q0.025 = stt_q0.025_state,
  #     stt_q0.25 = stt_q0.25_state,
  #     stt_q0.75 = stt_q0.75_state,
  #     stt_q0.975 = stt_q0.975_state
  #   )
  #
  #   ### STT state2
  #   s_freq <- length(simulation_outputs[[1]][[1]]$stt_all[, 1])
  #   complete_arr <- array(dim = c(s_freq, 7, replicates))
  #
  #   for (x in 1:replicates) {
  #     sum_endemics <- simulation_outputs[[x]][[1]]$stt_all[, "nA2"] +
  #       simulation_outputs[[x]][[1]]$stt_all[, "nC2"]
  #     total <- simulation_outputs[[x]][[1]]$stt_all[, "nA2"] +
  #       simulation_outputs[[x]][[1]]$stt_all[, "nC2"] +
  #       simulation_outputs[[x]][[1]]$stt_all[, "nI2"]
  #     nI2 <- simulation_outputs[[x]][[1]]$stt_all[,"nI2"]
  #     nA2 <- simulation_outputs[[x]][[1]]$stt_all[,"nA2"]
  #     nC2 <- simulation_outputs[[x]][[1]]$stt_all[,"nC2"]
  #     complete_arr[, , x] <- cbind(simulation_outputs[[x]][[1]]$stt_all[,'Time'],nI2,nA2,nC2,sum_endemics,total)
  #
  #   }
  #
  #   stt_average_state1 <- apply(complete_arr, c(1, 2), stats::median)
  #   stt_q0.025_state1 <- apply(complete_arr, c(1, 2), stats::quantile, 0.025)
  #   stt_q0.25_state1 <- apply(complete_arr, c(1, 2), stats::quantile, 0.25)
  #   stt_q0.75_state1 <- apply(complete_arr, c(1, 2), stats::quantile, 0.75)
  #   stt_q0.975_state1 <- apply(complete_arr, c(1, 2), stats::quantile, 0.975)
  #
  #   colnames(stt_average_state1) <- c(
  #     "Time",
  #     "nI1",
  #     "nA1",
  #     "nC1",
  #     "present",
  #     "Endemic",
  #     "Total"
  #   )
  #   colnames(stt_q0.025_state1) <- c(
  #     "Time",
  #     "nI1",
  #     "nA1",
  #     "nC1",
  #     "present",
  #     "Endemic",
  #     "Total"
  #   )
  #   colnames(stt_q0.25_state1) <- c(
  #     "Time",
  #     "nI1",
  #     "nA1",
  #     "nC1",
  #     "present",
  #     "Endemic",
  #     "Total"
  #   )
  #   colnames(stt_q0.75_state1) <- c(
  #     "Time",
  #     "nI1",
  #     "nA1",
  #     "nC1",
  #     "present",
  #     "Endemic",
  #     "Total"
  #   )
  #   colnames(stt_q0.975_state1) <- c(
  #     "Time",
  #     "nI1",
  #     "nA1",
  #     "nC1",
  #     "present",
  #     "Endemic",
  #     "Total"
  #   )
  #
  #   state1_species <- list(
  #     stt_average = stt_average_state1,
  #     stt_q0.025 = stt_q0.025_state1,
  #     stt_q0.25 = stt_q0.25_state1,
  #     stt_q0.75 = stt_q0.75_state1,
  #     stt_q0.975 = stt_q0.975_state1
  #   )
  #   return(list(
  #     all_species = all_species,
  #     state0_species = state0_species,
  #     state1_species = state1_species)
  #   )
  # }
}



#' Create the Species-Through-Time plot. This is used to visualize
#' the output of \code{\link{DAISIE_sim}}
#'
#' @param plot_plus_one Boolean to indicate to plot all values plus one.
#'   Set to \code{TRUE} for default behavior.
#'   Set to \code{FALSE} to plot all values without adding one.
#'   Only works when there is one type of species
#' @param plot_lists List of lists containing average and quantile species
#'   through time.
#' @param type String to indicate if stt of all species or all possible stt
#'   should be plotted. Default is \code{"all_species"}.
#' @param time the time span simulated
#' @param Tpars A named list containing diversification rates considering two trait states:
#' \itemize{
#'   \item{[1]:A numeric with the per capita transition rate with state1}
#'   \item{[2]:A numeric with the per capita immigration rate with state2}
#'   \item{[3]:A numeric with the per capita extinction rate with state2}
#'   \item{[4]:A numeric with the per capita anagenesis rate with state2}
#'   \item{[5]:A numeric with the per capita cladogenesis rate with state2}
#'   \item{[6]:A numeric with the per capita transition rate with state2}
#'   \item{[7]:A numeric with the number of species with trait state 2 on mainland}
#' }
#' @export
DAISIE_plot_stt <- function(
  plot_plus_one = TRUE,
  time,
  plot_lists = plot_lists,
  type = type,
  Tpars = NULL
) {
  # Plot the y axis iff plus one
  y_axis_type <- 'n'
  y_axis_label <- "No of species"
  if (plot_plus_one == TRUE) {
    y_axis_type <- 's'
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

