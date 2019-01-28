#' Create the Species-Through-Time plot. This is used to visualize
#' the output of \code{\link{DAISIE_sim}}
#' @param plot_plus_one Boolean to indicate to plot all values plus one.
#'   Set to \code{TRUE} for default behavior.
#'   Set to \code{FALSE} to plot all values without adding one.
#'   Only works when there is one type of species
#' @param time the time span simulated
#' @param stt_q0.025_all STT 2.5\% quantile. 
#'   Must be a data frame with columns name \code{Time}, \code{nI}, 
#'   \code{Endemic} and \code{Total}
#' @param stt_q0.25_all STT 25\% quantile
#'   Must be a data frame with columns name \code{Time}, \code{nI}, 
#'   \code{Endemic} and \code{Total}
#' @param stt_average_all STT average
#'   Must be a data frame with columns name \code{Time}, \code{nI}, 
#'   \code{Endemic} and \code{Total}
#' @param stt_q0.75_all STT 75\% quantile
#'   Must be a data frame with columns name \code{Time}, \code{nI}, 
#'   \code{Endemic} and \code{Total}
#' @param stt_q0.975_all STT 97.5\% quantile
#'   Must be a data frame with columns name \code{Time}, \code{nI}, 
#'   \code{Endemic} and \code{Total}
DAISIE_plot_stt <- function(
  plot_plus_one = FALSE,
  time,
  stt_q0.025_all,
  stt_q0.25_all,
  stt_average_all,
  stt_q0.75_all,
  stt_q0.975_all
) {
  # Plot the y axis iff plus one
  y_axis_type <- 'n'
  y_axis_label <- "No of species" 
  if (plot_plus_one == TRUE) {
    y_axis_type <- 's'
    y_axis_label <- "No of species + 1" 
  }
  
  suppressWarnings(
    plot(NULL, NULL, xlim = rev(c(0, time)), ylim = c(1, max(stt_q0.975_all)), 
      ylab = y_axis_label, 
      bty = "l", xaxs = "i", xlab = "Time before present", 
      main = "Species-through-time - All species", 
      log = "y", cex.lab = 1.2, cex.main = 1.2, cex.axis = 1.2,
      yaxt = y_axis_type
    )
  )
  polygon(c(stt_average_all[, "Time"], rev(stt_average_all[, "Time"])), c(stt_q0.025_all[, "Total"] + 
      1, rev(stt_q0.975_all[, "Total"] + 1)), col = "light grey", border = NA)
  polygon(c(stt_average_all[, "Time"], rev(stt_average_all[, "Time"])), c(stt_q0.25_all[, "Total"] + 
      1, rev(stt_q0.75_all[, "Total"] + 1)), col = "dark grey", border = NA)
  lines(stt_average_all[, "Time"], stt_average_all[, "Total"] + 1, lwd = 2)
  lines(stt_average_all[, "Time"], stt_average_all[, "nI"] + 1, lwd = 2, col = "cyan3")
  lines(stt_average_all[, "Time"], stt_average_all[, "Endemic"] + 1, lwd = 2, col = "dodgerblue1")

  legend_names <- c("Total", "Non-endemic", "Endemic")
  legend_colors <- c("black", "cyan3", "dodgerblue1")
  legend(
    time, max(stt_q0.975_all), legend_names, lty = 1, lwd = 2, 
    col = legend_colors, cex = 1.2, border = NA, bty = "n"
  )
  if (plot_plus_one == FALSE) {
    y_axis_values <- c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)  
    axis(2, at = y_axis_values, labels = y_axis_values - 1)
  }
}