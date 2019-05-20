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
#'   \code{\link{DAISIE_sim}} with \code{format = TRUE} option). Minimally, this must be 
#'   a list, that has as much elements as replicates. Each element must be a
#'   list with the elements \code{island_age}, \code{not_present} 
#'   and \code{stt_all}. \code{stt_all} must be a data frame with
#'   the column names \code{Time}, \code{nI}, \code{nA}, \code{nC} 
#'   and \code{present}.
#' @param use_dev_new Boolean to indicate to use \code{dev.new} before plotting. 
#'   Set to \code{TRUE} for default behavior.
#'   Set to \code{FALSE} when plotting within a vignette
#' @param plot_plus_one Boolean to indicate to plot all values plus one.
#'   Set to \code{TRUE} for default behavior.
#'   Set to \code{FALSE} to plot all values without adding one.
#'   Only works when there is one type of species
#' @return R plot.
#' @author Luis Valente
#' @seealso \code{\link{DAISIE_sim}} \code{\link{DAISIE_format_CS}}
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
#'   island_replicates = islands_1type_1000reps,
#'   use_dev_new = FALSE
#'   )
#' 
#' 
#' ### Plot island with type 1 and type 2
#' utils::data(islands_2types_1000reps)
#' DAISIE_plot_sims(
#'   island_replicates = islands_2types_1000reps,
#'   use_dev_new = FALSE
#'   )
#' 
#' 
#' 
DAISIE_plot_sims <- function(
  island_replicates, 
  use_dev_new = TRUE,
  plot_plus_one = TRUE
) {
  time <- max(island_replicates[[1]][[1]]$stt_all[, 1])
  # Prepare dataset
  outs <- DAISIE_prepare_data_plotting(island_replicates)
  
  if (use_dev_new == TRUE) {
    grDevices::dev.new(width = 12, height = 4)
  }
  graphics::par(mfrow = c(1, 3))
  
  if (is.null(island_replicates[[1]][[1]]$stt_type1) == FALSE) {
    
    # All species
    DAISIE_plot_stt(
      plot_plus_one = plot_plus_one,
      time = time,
      stt_q0.025 = stt_q0.025_all,
      stt_q0.25 = stt_q0.25_all,
      stt_average = stt_average_all,
      stt_q0.75 = stt_q0.75_all,
      stt_q0.975 = stt_q0.975_all
    )
    
    # Type 1 species
    DAISIE_plot_stt(
      plot_plus_one = plot_plus_one,
      time = time,
      stt_q0.025 = stt_q0.025_type1,
      stt_q0.25 = stt_q0.25_type1,
      stt_average = stt_average_type1,
      stt_q0.75 = stt_q0.75_type1,
      stt_q0.975 = stt_q0.975_type1
    )
    
    # Type 2 species
    DAISIE_plot_stt(
      plot_plus_one = plot_plus_one,
      time = time,
      stt_q0.025 = stt_q0.025_type2,
      stt_q0.25 = stt_q0.25_type2,
      stt_average = stt_average_type2,
      stt_q0.75 = stt_q0.75_type2,
      stt_q0.975 = stt_q0.975_type2
    )
    
  } else {
    # Only plot all species
    if (use_dev_new == TRUE) {
      # Default behavior to open a new device, which hurts vignettes
      grDevices::dev.new(width = 6, height = 6)
    }
    
    graphics::par(mfrow = c(1, 1))
    
    DAISIE_plot_stt(
      plot_plus_one = plot_plus_one,
      time = time,
      stt_q0.025 = stt_q0.025_all,
      stt_q0.25 = stt_q0.25_all,
      stt_average = stt_average_all,
      stt_q0.75 = stt_q0.75_all,
      stt_q0.975 = stt_q0.975_all
    )
  }
}
