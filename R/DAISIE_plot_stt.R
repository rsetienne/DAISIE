#' Prepare input for DAISIE_stt
#'
#' @inheritParams DAISIE_plot_sims 
#' @param simulation_outputs A list with matrices? of simulation produced by
#' DAISIE_sim.  
#' @seealso \code{\link{DAISIE_plot_stt}}, \code{\link{DAISIE_plot_sims}}
#' @examples 
#' utils::data("islands_1type_1000reps", package = "DAISIE")
#' DAISIE:::DAISIE_convert_to_classic_plot(islands_1type_1000reps)
#'
#'
#' @return a list with wrangled data to be used for plotting STT plots with
#' DAISIE_plot_stt
#' 
DAISIE_convert_to_classic_plot <- function(simulation_outputs) {
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
    sum_endemics<-simulation_outputs[[x]][[1]]$stt_all[,"nA0"]+simulation_outputs[[x]][[1]]$stt_all[,"nC0"]+simulation_outputs[[x]][[1]]$stt_all[,"nA1"]+simulation_outputs[[x]][[1]]$stt_all[,"nC1"]
    total<-simulation_outputs[[x]][[1]]$stt_all[,"nA0"]+simulation_outputs[[x]][[1]]$stt_all[,"nC0"]+simulation_outputs[[x]][[1]]$stt_all[,"nI0"]+simulation_outputs[[x]][[1]]$stt_all[,"nA1"]+simulation_outputs[[x]][[1]]$stt_all[,"nC1"]+simulation_outputs[[x]][[1]]$stt_all[,"nI1"]
    nItotal<-simulation_outputs[[x]][[1]]$stt_all[,"nI0"]+simulation_outputs[[x]][[1]]$stt_all[,"nI1"]
    nAtotal<-simulation_outputs[[x]][[1]]$stt_all[,"nA0"]+simulation_outputs[[x]][[1]]$stt_all[,"nA1"]
    nCtotal<-simulation_outputs[[x]][[1]]$stt_all[,"nC0"]+simulation_outputs[[x]][[1]]$stt_all[,"nC1"]
    complete_arr[,,x]<-cbind(simulation_outputs[[x]][[1]]$stt_all[,'Time'],nItotal,nAtotal,nCtotal,sum_endemics,total)
  }
  
  stt_average_all <- apply(complete_arr, c(1, 2), stats::median)
  testit::assert(stt_average_all == DAISIE::DAISIE_extract_stt_median(simulation_outputs))
  stt_q0.025_all <- apply(complete_arr, c(1, 2), stats::quantile, 0.025)
  stt_q0.25_all <- apply(complete_arr, c(1, 2), stats::quantile, 0.25)
  stt_q0.75_all <- apply(complete_arr, c(1, 2), stats::quantile, 0.75)
  stt_q0.975_all <- apply(complete_arr, c(1, 2), stats::quantile, 0.975)
  
  colnames(stt_average_all)<-c("Time","nItotal","nAtotal","nCtotal","Endemic","Total")
  colnames(stt_q0.025_all)<-c("Time","nItotal","nAtotal","nCtotal","Endemic","Total")
  colnames(stt_q0.25_all)<-c("Time","nItotal","nAtotal","nCtotal","Endemic","Total")
  colnames(stt_q0.75_all)<-c("Time","nItotal","nAtotal","nCtotal","Endemic","Total")
  colnames(stt_q0.975_all)<-c("Time","nItotal","nAtotal","nCtotal","Endemic","Total")
  
  all_species <- list(
    stt_average = stt_average_all,
    stt_q0.025 = stt_q0.025_all,
    stt_q0.25 = stt_q0.25_all,
    stt_q0.75 = stt_q0.75_all,
    stt_q0.975 = stt_q0.975_all
  )
  
 # if (is.null(simulation_outputs[[1]][[1]]$stt_type1) == FALSE) {
    
    ### STT state0
    s_freq <- length(simulation_outputs[[1]][[1]]$stt_all[, 1])
    complete_arr <- array(dim = c(s_freq, 7, replicates))
    
    for (x in 1:replicates) {
      sum_endemics <- simulation_outputs[[x]][[1]]$stt_all[, "nA0"] + simulation_outputs[[x]][[1]]$stt_all[, 
                                                                                                            "nC0"]
      total <- simulation_outputs[[x]][[1]]$stt_all[, "nA0"] + simulation_outputs[[x]][[1]]$stt_all[, 
                                                                                                     "nC0"] + simulation_outputs[[x]][[1]]$stt_all[, "nI0"]
      nI0 <- simulation_outputs[[x]][[1]]$stt_all[,"nI0"]
      nA0 <- simulation_outputs[[x]][[1]]$stt_all[,"nA0"]
      nC0 <- simulation_outputs[[x]][[1]]$stt_all[,"nC0"]
      complete_arr[, , x] <- cbind(simulation_outputs[[x]][[1]]$stt_all[,'Time'],nI0,nA0,nC0,sum_endemics,total)
    }
    
    
    stt_average_state0 <- apply(complete_arr, c(1, 2), stats::median)
    stt_q0.025_state0 <- apply(complete_arr, c(1, 2), stats::quantile, 0.025)
    stt_q0.25_state0 <- apply(complete_arr, c(1, 2), stats::quantile, 0.25)
    stt_q0.75_state0 <- apply(complete_arr, c(1, 2), stats::quantile, 0.75)
    stt_q0.975_state0 <- apply(complete_arr, c(1, 2), stats::quantile, 0.975)
    
    colnames(stt_average_state0) <- c(
      "Time", 
      "nI0", 
      "nA0", 
      "nC0", 
      "present", 
      "Endemic", 
      "Total"
    )
    colnames(stt_q0.025_state0) <- c(
      "Time", 
      "nI0", 
      "nA0", 
      "nC0", 
      "present", 
      "Endemic", 
      "Total"
    )
    colnames(stt_q0.25_state0) <- c(
      "Time", 
      "nI0", 
      "nA0", 
      "nC0", 
      "present", 
      "Endemic", 
      "Total"
    )
    colnames(stt_q0.75_state0) <- c(
      "Time",
      "nI0",
      "nA0",
      "nC0",
      "present",
      "Endemic", 
      "Total"
    )
    colnames(stt_q0.975_state0) <- c(
      "Time",
      "nI0",
      "nA0",
      "nC0",
      "present",
      "Endemic",
      "Total"
    )
    
    state0_species <- list(
      stt_average = stt_average_state0,
      stt_q0.025 = stt_q0.025_state0,
      stt_q0.25 = stt_q0.25_state0,
      stt_q0.75 = stt_q0.75_state0,
      stt_q0.975 = stt_q0.975_state0
    )
    
    ### STT state1
    s_freq <- length(simulation_outputs[[1]][[1]]$stt_all[, 1])
    complete_arr <- array(dim = c(s_freq, 7, replicates))
    
    for (x in 1:replicates) {
      sum_endemics <- simulation_outputs[[x]][[1]]$stt_all[, "nA1"] +
        simulation_outputs[[x]][[1]]$stt_all[,"nC"]
      total <- simulation_outputs[[x]][[1]]$stt_all[, "nA"] +
        simulation_outputs[[x]][[1]]$stt_all[, 
                                              "nC"] +
        simulation_outputs[[x]][[1]]$stt_all[, "nI"]
      complete_arr[, , x] <- cbind(
        simulation_outputs[[x]][[1]]$stt_all,
        sum_endemics,
        total
      )
    }
    
    stt_average_state1 <- apply(complete_arr, c(1, 2), stats::median)
    stt_q0.025_state1 <- apply(complete_arr, c(1, 2), stats::quantile, 0.025)
    stt_q0.25_state1 <- apply(complete_arr, c(1, 2), stats::quantile, 0.25)
    stt_q0.75_state1 <- apply(complete_arr, c(1, 2), stats::quantile, 0.75)
    stt_q0.975_state1 <- apply(complete_arr, c(1, 2), stats::quantile, 0.975)
    
    colnames(stt_average_state1) <- c(
      "Time",
      "nI1",
      "nA1",
      "nC1",
      "present",
      "Endemic",
      "Total"
    )
    colnames(stt_q0.025_state1) <- c(
      "Time", 
      "nI1", 
      "nA1", 
      "nC1", 
      "present", 
      "Endemic", 
      "Total"
    )
    colnames(stt_q0.25_state1) <- c(
      "Time",
      "nI1",
      "nA1",
      "nC1",
      "present",
      "Endemic",
      "Total"
    )
    colnames(stt_q0.75_state1) <- c(
      "Time",
      "nI1",
      "nA1",
      "nC1",
      "present",
      "Endemic",
      "Total"
    )
    colnames(stt_q0.975_state1) <- c(
      "Time",
      "nI1",
      "nA1",
      "nC1",
      "present",
      "Endemic",
      "Total"
    )
    
    state1_species <- list(
      stt_average = stt_average_state1,
      stt_q0.025 = stt_q0.025_state1,
      stt_q0.25 = stt_q0.25_state1,
      stt_q0.75 = stt_q0.75_state1,
      stt_q0.975 = stt_q0.975_state1
    )
    
    return(list(
      all_species = all_species,
      state0_species = state0_species,
      state1_species = state1_species
      )
    )
} 
  # } else {
  #   return(list(
  #     all_species = all_species,
  #     state0_species = NULL,
  #     state1_species = NULL
  #   )
  #   )
  # }

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
DAISIE_plot_stt <- function(
  plot_plus_one = TRUE,
  time,
  plot_lists = plot_lists,
  type = type
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
  graphics::lines(stt$stt_average[, "Time"], stt$stt_average[, "nItotal"] + 1, lwd = 2, col = "cyan3")
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

