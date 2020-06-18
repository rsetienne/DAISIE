#' Plot STT and overlay additional STT curves.
#'
#' @inheritParams default_params_doc
#'
#' @seealso \code{\link{DAISIE_plot_sims}}, \code{\link{DAISIE_plot_stt}},
#'   \code{\link{DAISIE_convert_to_classic_plot}}
#' @family plotting
#' @return Standard \code{\link{DAISIE_plot_stt}} with overlaid additional
#'   STT curves for comparison.
#' @export
#'
#' @author Pedro Neves
DAISIE_plot_comparison_stts <- function(
  time,
  plot_lists_simulations,
  plot_lists_simulations_MLE,
  type,
  kind_of_plot = "line"
) {
  valid_types <- c("all_species", "type1_species", "type2_species")
  if (all(type != valid_types)) {
    stop(
      "type should be 'all_species', 'type1_species' or 'type2_species'. \n",
      "Actual value: ", type
    )
  }

  valid_kinds <- c("line", "shade")
  if (all(kind_of_plot != valid_kinds)) {
    stop(
      "type should be 'line' or 'shade'. \n",
      "Actual value: ", kind_of_plot
    )
  }
  y_axis_type <- "s"
  y_axis_label <- "No of species + 1"

  stt_simulations <- plot_lists_simulations[[type]]
  stt_simulations_MLE <- plot_lists_simulations_MLE
  #This must be a list with the 10 indep lines
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

  if (kind_of_plot == "shade") {
    stt_simulations_MLE <- plot_lists_simulations_MLE[[type]]
    graphics::polygon(
      c(stt_simulations_MLE$stt_average[, "Time"],
        rev(stt_simulations_MLE$stt_average[, "Time"])),
      c(stt_simulations_MLE$stt_q0.025[, "Total"] + 1,
        rev(stt_simulations_MLE$stt_q0.975[, "Total"] + 1)),
      col = "light green", border = NA
    )
  }

  # graphics::polygon(
  #   c(stt_simulations$stt_average[, "Time"],
  #     rev(stt_simulations$stt_average[, "Time"])),
  #   c(stt_simulations$stt_q0.025[, "Total"] + 1,
  #     rev(stt_simulations$stt_q0.975[, "Total"] + 1)),
  #   col = "light grey", border = NA
  # )
  # graphics::polygon(
  #   c(
  #     stt_simulations$stt_average[, "Time"],
  #     rev(stt_simulations$stt_average[, "Time"])
  #   ), c(
  #     stt_simulations$stt_q0.25[, "Total"] +
  #       1, rev(stt_simulations$stt_q0.75[, "Total"] + 1)
  #   ), col = "dark grey", border = NA
  # )

  graphics::lines(
    stt_simulations$stt_average[, "Time"],
    stt_simulations$stt_average[, "Total"] + 1,
    lwd = 2
  )
  graphics::lines(
    stt_simulations$stt_average[, "Time"],
    stt_simulations$stt_average[, "nI"] + 1,
    lwd = 2,
    col = "cyan3"
  )
  graphics::lines(
    stt_simulations$stt_average[, "Time"],
    stt_simulations$stt_average[, "Endemic"] + 1,
    lwd = 2,
    col = "dodgerblue1"
  )

  # Plot MLE obtained simulations
  if (kind_of_plot == "line") {
    for (n_replicate in seq_along(stt_simulations_MLE)) {
      graphics::lines(
        stt_simulations_MLE[[n_replicate]][[1]]$stt_average[, "Time"],
        stt_simulations_MLE[[n_replicate]][[1]]$stt_average[, "Total"] + 1,
        lwd = 1, col = "darkgreen"
      )
    }
  }
  # } else if (kind_of_plot == "shade") {
  #   stt_simulations_MLE <- plot_lists_simulations_MLE[[type]]
  #   graphics::polygon(
  #     c(stt_simulations_MLE$stt_average[, "Time"],
  #       rev(stt_simulations_MLE$stt_average[, "Time"])),
  #     c(stt_simulations_MLE$stt_q0.025[, "Total"] + 1,
  #       rev(stt_simulations_MLE$stt_q0.975[, "Total"] + 1)),
  #     col = "light green", border = NA
  #   )
  #
  # Write legend
  legend_names <- c("Total", "Non-endemic", "Endemic")
  legend_colors <- c("black", "cyan3", "dodgerblue1")
  graphics::legend(
    time,
    max(stt_simulations$stt_q0.975),
    legend_names,
    lty = 1,
    lwd = 2,
    col = legend_colors,
    cex = 1.2,
    border = NA,
    bty = "n"
  )
}
