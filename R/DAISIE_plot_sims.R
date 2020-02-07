#' @title Plot island species-through-time (STT) plots
#' @description Produces STT plots. If only one type of species is present in
#' the simulated islands, STT is plotted for all species. If two types are
#' present, three plots are produced: STT for all, STT for type 1 and STT
#' for type 2.
#'
#' R plots with number of total, endemic and non-endemic STTs for different
#' types of species for the entire time span the islands were simulated.
#' 2.5-97.5th percentiles are plotted in light grey, 25-75th percentiles
#' plotted in dark grey.
#'
#' @param island_replicates Island replicates in DAISIE format (produced in
#'   \code{\link{DAISIE_sim_constant_rate}},
#'   \code{\link{DAISIE_sim_time_dependent}}, or
#'   \code{\link{DAISIE_sim_constant_rate_shift}} with \code{format = TRUE}
#'   option). Minimally, this must be a list, that has as much elements as
#'    replicates. Each element must be a list with the elements
#'    \code{island_age}, \code{not_present} and \code{stt_all}.
#'    \code{stt_all} must be a data frame with the column names \code{Time},
#'     \code{nI}, \code{nA}, \code{nC} and \code{present}.
#' @param plot_plus_one Boolean to indicate to plot all values plus one.
#'   Set to \code{TRUE} for default behavior.
#'   Set to \code{FALSE} to plot all values without adding one.
#' @param type String to indicate if stt of all species or all possible stt
#'   should be plotted. Default is \code{"all_species"}.
#' @param sample_freq Specifies the number of units times should be divided by
#' for plotting purposes. Larger values will lead to plots with higher
#' resolution, but will also run slower.
#' @return R plot.
#' @author Luis Valente
#' @seealso \code{\link{DAISIE_sim_constant_rate}},
#' \code{\link{DAISIE_sim_time_dependent}},
#' \code{\link{DAISIE_sim_constant_rate_shift}}, \code{\link{DAISIE_format_CS}}
#' @references Valente, L.M., A.B. Phillimore and R.S. Etienne (2015).
#' Equilibrium and non-equilibrium dynamics simultaneously operate in the
#' Galapagos islands. Ecology Letters 18: 844-852.
#' @keywords models
#' @export
#' @examples
#'
#'
#' ### Plot islands with single process (only one type of species)
#' utils::data(islands_1type_1000reps)
#' DAISIE_plot_sims(
#'   island_replicates = islands_1type_1000reps
#'   )
#'
#'
#' ### Plot island with type 1 and type 2
#' utils::data(islands_2types_1000reps)
#' DAISIE_plot_sims(
#'   island_replicates = islands_2types_1000reps
#'   )
#'
#'
#'
DAISIE_plot_sims <- function(
  island_replicates,
  plot_plus_one = TRUE,
  type = "all_species",
  sample_freq = 25
) {
  time <- max(island_replicates[[1]][[1]]$stt_all[, 1])
  if (sample_freq != Inf) {
    # Prepare dataset
    plot_lists <- DAISIE_convert_to_classic_plot(island_replicates)
  } else {
    stop("Plotting STT with sample_freq = Inf not yet available. \n")
  }
  if (type == "all") {
    types <- names(plot_lists)
  } else {
    types <- type
  }
  num_plots <- sum(!sapply(plot_lists[types], FUN = is.null))
  graphics::par(mfrow = c(1, num_plots))
  for (type_here in types) {
    DAISIE_plot_stt(
      plot_plus_one = plot_plus_one,
      time = time,
      plot_lists = plot_lists,
      type = type_here)
  }
}
