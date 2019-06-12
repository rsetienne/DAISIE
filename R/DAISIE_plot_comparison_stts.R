#' Plot STT and overlay additional STT curves. 
#'
#' @inheritParams DAISIE_plot_stt 
#' @param plot_lists_simulations List with simulation output after parsing by
#' \code{DAISIE_prepare_data_plotting}
#' @param plot_lists_simulations_MLE List with simulation output after parsing by
#' \code{DAISIE_prepare_data_plotting}, but obtained by simulating MLE output
#' @param type Character vector stating if \code{"all-species"}, \code{"type1"}
#' or \code{"type2"} should be plotted
#' 
#' @seealso \code{\link{DAISIE_plot_sims}}, \code{\link{DAISIE_plot_stt}}, 
#' \code{\link{DAISIE_prepare_data_plotting}}
#' 
#' @return Standard \code{\link{DAISIE_plot_stt}} with overlaid additional
#' STT curves for comparison.
#' @export 
#'
DAISIE_plot_comparison_stts <- function(
  time,
  plot_lists_simulations,
  plot_lists_simulations_MLE,
  type = type
) {
  
  
  
  
  y_axis_type <- 's'
  y_axis_label <- "No of species + 1"
  
  stt_simulations <- plot_lists_simulations[[type]]
  stt_simulations_MLE <- plot_lists_simulations_MLE # This must be a list with the 10 indep lines
  if (is.null(stt_simulations)) {
    return()
  }
  
  
  
  # Plot standard stt (start by opening empty canvas)
  suppressWarnings(
    graphics::plot(
      x = NULL,
      y = NULL,
      xlim = rev(c(0, time)),
      ylim = c(1, max(stt_simulations$stt_q0.975)),
      ylab = y_axis_label,
      bty = "l",
      xaxs = "i",
      xlab = "Time before present",
      main = "Species-through-time - All species",
      log = "y",
      cex.lab = 1.2,
      cex.main = 1.2,
      cex.axis = 1.2,
      yaxt = y_axis_type
    )
  )
  graphics::polygon(
    c(stt_simulations$stt_average[, "Time"],
      rev(stt_simulations$stt_average[, "Time"])),
    c(stt_simulations$stt_q0.025[, "Total"] + 1,
      rev(stt_simulations$stt_q0.975[, "Total"] + 1)),
    col = "light grey", border = NA
  )
  graphics::polygon(c(stt_simulations$stt_average[, "Time"], rev(stt_simulations$stt_average[, "Time"])), c(stt_simulations$stt_q0.25[, "Total"] + 
                                                                                                              1, rev(stt_simulations$stt_q0.75[, "Total"] + 1)), col = "dark grey", border = NA)
  graphics::lines(stt_simulations$stt_average[, "Time"], stt_simulations$stt_average[, "Total"] + 1, lwd = 2)
  graphics::lines(stt_simulations$stt_average[, "Time"], stt_simulations$stt_average[, "nI"] + 1, lwd = 2, col = "cyan3")
  graphics::lines(stt_simulations$stt_average[, "Time"], stt_simulations$stt_average[, "Endemic"] + 1, lwd = 2, col = "dodgerblue1")
  
  # Plot MLE obtained simulations
  for (replicate in seq_along(stt_simulations_MLE)) {
    graphics::lines(
      stt_simulations_MLE[[replicate]][, "Time"],
      stt_simulations_MLE[[replicate]][, "Total"] + 1,
      lwd = 1, col = "darkgreen"
    )
  }
  
  # Write legend
  legend_names <- c("Total", "Non-endemic", "Endemic")
  legend_colors <- c("black", "cyan3", "dodgerblue1")
  graphics::legend(
    time, max(stt_simulations$stt_q0.975), legend_names, lty = 1, lwd = 2, 
    col = legend_colors, cex = 1.2, border = NA, bty = "n"
  )
}
