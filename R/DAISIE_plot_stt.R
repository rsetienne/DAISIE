#' Prepare input for DAISIE_stt
#'
#' @inheritParams DAISIE_plot_sims 
#'
#' @return
#' 
DAISIE_prepare_data_plotting <- function(island_replicates) {
  replicates <- length(island_replicates)
  
  
  ### STT ALL species
  s_freq <- length(island_replicates[[1]][[1]]$stt_all[, 1])
  complete_arr <- array(dim = c(s_freq, 6, replicates))
  
  for (x in 1:replicates) {
    sum_endemics <- island_replicates[[x]][[1]]$stt_all[, "nA"] + island_replicates[[x]][[1]]$stt_all[, 
                                                                                                      "nC"]
    total <- island_replicates[[x]][[1]]$stt_all[, "nA"] + island_replicates[[x]][[1]]$stt_all[, 
                                                                                               "nC"] + island_replicates[[x]][[1]]$stt_all[, "nI"]
    complete_arr[, , x] <- cbind(island_replicates[[x]][[1]]$stt_all[, c("Time", "nI", "nA", "nC")], 
                                 sum_endemics, total)
  }
  
  stt_average_all <- apply(complete_arr, c(1, 2), stats::median)
  testit::assert(stt_average_all == DAISIE::DAISIE_extract_stt_median(island_replicates))
  stt_q0.025_all <- apply(complete_arr, c(1, 2), stats::quantile, 0.025)
  stt_q0.25_all <- apply(complete_arr, c(1, 2), stats::quantile, 0.25)
  stt_q0.75_all <- apply(complete_arr, c(1, 2), stats::quantile, 0.75)
  stt_q0.975_all <- apply(complete_arr, c(1, 2), stats::quantile, 0.975)
  
  colnames(stt_average_all) <- c("Time", "nI", "nA", "nC", "Endemic", "Total")
  colnames(stt_q0.025_all) <- c("Time", "nI", "nA", "nC", "Endemic", "Total")
  colnames(stt_q0.25_all) <- c("Time", "nI", "nA", "nC", "Endemic", "Total")
  colnames(stt_q0.75_all) <- c("Time", "nI", "nA", "nC", "Endemic", "Total")
  colnames(stt_q0.975_all) <- c("Time", "nI", "nA", "nC", "Endemic", "Total")
  
  all_species_out <- list(
    stt_average_all = stt_average_all,
    stt_q0.025_all = stt_q0.025_all,
    stt_q0.25_all = stt_q0.25_all,
    stt_q0.75_all = stt_q0.75_all,
    stt_q0.975_all = stt_q0.975_all
  )
  
  if (is.null(island_replicates[[1]][[1]]$stt_type1) == FALSE) {
    
    ### STT TYPE1
    s_freq <- length(island_replicates[[1]][[1]]$stt_type1[, 1])
    complete_arr <- array(dim = c(s_freq, 7, replicates))
    
    for (x in 1:replicates) {
      sum_endemics <- island_replicates[[x]][[1]]$stt_type1[, "nA"] + island_replicates[[x]][[1]]$stt_type1[, 
                                                                                                            "nC"]
      total <- island_replicates[[x]][[1]]$stt_type1[, "nA"] + island_replicates[[x]][[1]]$stt_type1[, 
                                                                                                     "nC"] + island_replicates[[x]][[1]]$stt_type1[, "nI"]
      complete_arr[, , x] <- cbind(island_replicates[[x]][[1]]$stt_type1, sum_endemics, total)
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
    
    type1_species_out <- list(
      stt_average_type1 = stt_average_type1,
      stt_q0.025_type1 = stt_q0.025_type1,
      stt_q0.25_type1 = stt_q0.25_type1,
      stt_q0.75_type1 = stt_q0.75_type1,
      stt_q0.975_type1 = stt_q0.975_type1
    )
    
    ### STT TYPE2
    s_freq <- length(island_replicates[[1]][[1]]$stt_type2[, 1])
    complete_arr <- array(dim = c(s_freq, 7, replicates))
    
    for (x in 1:replicates) {
      sum_endemics <- island_replicates[[x]][[1]]$stt_type2[, "nA"] +
        island_replicates[[x]][[1]]$stt_type2[,"nC"]
      total <- island_replicates[[x]][[1]]$stt_type2[, "nA"] +
        island_replicates[[x]][[1]]$stt_type2[, 
                                              "nC"] +
        island_replicates[[x]][[1]]$stt_type2[, "nI"]
      complete_arr[, , x] <- cbind(island_replicates[[x]][[1]]$stt_type2, sum_endemics, total)
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
    
    type2_species_out <- list(
      stt_average_type2 = stt_average_type2,
      stt_q0.025_type2 = stt_q0.025_type2,
      stt_q0.25_type2 = stt_q0.25_type2,
      stt_q0.75_type2 = stt_q0.75_type2,
      stt_q0.975_type2 = stt_q0.975_type2
    )
    
    return(list(
      all_species_out = all_species_out,
      type1_species_out = type1_species_out,
      type2_species_out = type2_species_out)
    )
    
  } else {
    return(list(
      all_species_out = all_species_out,
      type1_species_out = NULL,
      type2_species_out = NULL)
    )
  }
}

#' Create the Species-Through-Time plot. This is used to visualize
#' the output of \code{\link{DAISIE_sim}}
#' @param plot_plus_one Boolean to indicate to plot all values plus one.
#'   Set to \code{TRUE} for default behavior.
#'   Set to \code{FALSE} to plot all values without adding one.
#'   Only works when there is one type of species
#' @param time the time span simulated
#' @param stt_q0.025 STT 2.5\% quantile. 
#'   Must be a data frame with columns name \code{Time}, \code{nI}, 
#'   \code{Endemic} and \code{Total}
#' @param stt_q0.25 STT 25\% quantile
#'   Must be a data frame with columns name \code{Time}, \code{nI}, 
#'   \code{Endemic} and \code{Total}
#' @param stt_average STT average
#'   Must be a data frame with columns name \code{Time}, \code{nI}, 
#'   \code{Endemic} and \code{Total}
#' @param stt_q0.75 STT 75\% quantile
#'   Must be a data frame with columns name \code{Time}, \code{nI}, 
#'   \code{Endemic} and \code{Total}
#' @param stt_q0.975 STT 97.5\% quantile
#'   Must be a data frame with columns name \code{Time}, \code{nI}, 
#'   \code{Endemic} and \code{Total}
DAISIE_plot_stt <- function(
  plot_plus_one = FALSE,
  time,
  stt_q0.025,
  stt_q0.25,
  stt_average,
  stt_q0.75,
  stt_q0.975
) {
  # Plot the y axis iff plus one
  y_axis_type <- 'n'
  y_axis_label <- "No of species" 
  if (plot_plus_one == TRUE) {
    y_axis_type <- 's'
    y_axis_label <- "No of species + 1" 
  }
  
  suppressWarnings(
    graphics::plot(NULL, NULL, xlim = rev(c(0, time)), ylim = c(1, max(stt_q0.975)), 
                   ylab = y_axis_label, 
                   bty = "l", xaxs = "i", xlab = "Time before present", 
                   main = "Species-through-time - All species", 
                   log = "y", cex.lab = 1.2, cex.main = 1.2, cex.axis = 1.2,
                   yaxt = y_axis_type
    )
  )
  graphics::polygon(c(stt_average[, "Time"], rev(stt_average[, "Time"])), c(stt_q0.025[, "Total"] + 
                                                                              1, rev(stt_q0.975[, "Total"] + 1)), col = "light grey", border = NA)
  graphics::polygon(c(stt_average[, "Time"], rev(stt_average[, "Time"])), c(stt_q0.25[, "Total"] + 
                                                                              1, rev(stt_q0.75[, "Total"] + 1)), col = "dark grey", border = NA)
  graphics::lines(stt_average[, "Time"], stt_average[, "Total"] + 1, lwd = 2)
  graphics::lines(stt_average[, "Time"], stt_average[, "nI"] + 1, lwd = 2, col = "cyan3")
  graphics::lines(stt_average[, "Time"], stt_average[, "Endemic"] + 1, lwd = 2, col = "dodgerblue1")
  
  legend_names <- c("Total", "Non-endemic", "Endemic")
  legend_colors <- c("black", "cyan3", "dodgerblue1")
  graphics::legend(
    time, max(stt_q0.975), legend_names, lty = 1, lwd = 2, 
    col = legend_colors, cex = 1.2, border = NA, bty = "n"
  )
  if (plot_plus_one == FALSE) {
    y_axis_values <- c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)  
    graphics::axis(2, at = y_axis_values, labels = y_axis_values - 1)
  }
}