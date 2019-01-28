#' Plot island species-through-time (STT) plots
#' 
#' Produces STT plots. If only one type of species is present in the simulated
#' islands, STT is plotted for all species. If two types are present, three
#' plots are produced: STT for all, STT for type 1 and STT for type 2.
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
#' @examples
#'  
#' 
#' ### Plot islands with single process (only one type of species)
#' data(islands_1type_1000reps)
#' DAISIE_plot_sims(island_replicates = islands_1type_1000reps)
#' 
#' 
#' ### Plot island with type 1 and type 2
#' data(islands_2types_1000reps)
#' DAISIE_plot_sims(island_replicates = islands_2types_1000reps)
#' 
#' 
#' # Plot islands with single process
#' # Start counting from zero on the Y axis
#' DAISIE_plot_sims(
#' island_replicates = islands_1type_1000reps, 
#'   plot_plus_one = FALSE
#' )
DAISIE_plot_sims <- function(
  island_replicates, 
  use_dev_new = TRUE,
  plot_plus_one = TRUE
) {
  
  replicates <- length(island_replicates)
  time <- max(island_replicates[[1]][[1]]$stt_all[, 1])
  
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
  
  stt_average_all <- apply(complete_arr, c(1, 2), median)
  testit::assert(stt_average_all == DAISIE_extract_stt_median(island_replicates))
  stt_q0.025_all <- apply(complete_arr, c(1, 2), quantile, 0.025)
  stt_q0.25_all <- apply(complete_arr, c(1, 2), quantile, 0.25)
  stt_q0.75_all <- apply(complete_arr, c(1, 2), quantile, 0.75)
  stt_q0.975_all <- apply(complete_arr, c(1, 2), quantile, 0.975)

  colnames(stt_average_all) <- c("Time", "nI", "nA", "nC", "Endemic", "Total")
  colnames(stt_q0.025_all) <- c("Time", "nI", "nA", "nC", "Endemic", "Total")
  colnames(stt_q0.25_all) <- c("Time", "nI", "nA", "nC", "Endemic", "Total")
  colnames(stt_q0.75_all) <- c("Time", "nI", "nA", "nC", "Endemic", "Total")
  colnames(stt_q0.975_all) <- c("Time", "nI", "nA", "nC", "Endemic", "Total")
    
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
    
    
    stt_average_type1 <- apply(complete_arr, c(1, 2), median)
    stt_q0.025_type1 <- apply(complete_arr, c(1, 2), quantile, 0.025)
    stt_q0.25_type1 <- apply(complete_arr, c(1, 2), quantile, 0.25)
    stt_q0.75_type1 <- apply(complete_arr, c(1, 2), quantile, 0.75)
    stt_q0.975_type1 <- apply(complete_arr, c(1, 2), quantile, 0.975)
    
    colnames(stt_average_type1) <- c("Time", "nI", "nA", "nC", "present", "Endemic", "Total")
    colnames(stt_q0.025_type1) <- c("Time", "nI", "nA", "nC", "present", "Endemic", "Total")
    colnames(stt_q0.25_type1) <- c("Time", "nI", "nA", "nC", "present", "Endemic", "Total")
    colnames(stt_q0.75_type1) <- c("Time", "nI", "nA", "nC", "present", "Endemic", "Total")
    colnames(stt_q0.975_type1) <- c("Time", "nI", "nA", "nC", "present", "Endemic", "Total")
    
    ### STT TYPE2
    s_freq <- length(island_replicates[[1]][[1]]$stt_type2[, 1])
    complete_arr <- array(dim = c(s_freq, 7, replicates))
    
    for (x in 1:replicates) {
      sum_endemics <- island_replicates[[x]][[1]]$stt_type2[, "nA"] + island_replicates[[x]][[1]]$stt_type2[, 
                                                                                                            "nC"]
      total <- island_replicates[[x]][[1]]$stt_type2[, "nA"] + island_replicates[[x]][[1]]$stt_type2[, 
                                                                                                     "nC"] + island_replicates[[x]][[1]]$stt_type2[, "nI"]
      complete_arr[, , x] <- cbind(island_replicates[[x]][[1]]$stt_type2, sum_endemics, total)
    }
    
    stt_average_type2 <- apply(complete_arr, c(1, 2), median)
    stt_q0.025_type2 <- apply(complete_arr, c(1, 2), quantile, 0.025)
    stt_q0.25_type2 <- apply(complete_arr, c(1, 2), quantile, 0.25)
    stt_q0.75_type2 <- apply(complete_arr, c(1, 2), quantile, 0.75)
    stt_q0.975_type2 <- apply(complete_arr, c(1, 2), quantile, 0.975)
    
    colnames(stt_average_type2) <- c("Time", "nI", "nA", "nC", "present", "Endemic", "Total")
    colnames(stt_q0.025_type2) <- c("Time", "nI", "nA", "nC", "present", "Endemic", "Total")
    colnames(stt_q0.25_type2) <- c("Time", "nI", "nA", "nC", "present", "Endemic", "Total")
    colnames(stt_q0.75_type2) <- c("Time", "nI", "nA", "nC", "present", "Endemic", "Total")
    colnames(stt_q0.975_type2) <- c("Time", "nI", "nA", "nC", "present", "Endemic", "Total")
    
    if (use_dev_new == TRUE) {
      dev.new(width = 12, height = 4)
    }
    par(mfrow = c(1, 3))
        
    # Could use DAISIE_plot_stt here one day...
    suppressWarnings(plot(NULL, NULL, xlim = rev(c(0, time)), ylim = c(1, max(stt_q0.975_all)), ylab = "No of species + 1", 
        bty = "l", xaxs = "i", xlab = "Time before present", main = "Species-through-time - All species", 
        log = "y", cex.lab = 1.5, cex.main = 1.2, cex.axis = 1.2))
    polygon(c(stt_average_all[, "Time"], rev(stt_average_all[, "Time"])), c(stt_q0.025_all[, "Total"] + 
        1, rev(stt_q0.975_all[, "Total"] + 1)), col = "light grey", border = NA)
    polygon(c(stt_average_all[, "Time"], rev(stt_average_all[, "Time"])), c(stt_q0.25_all[, "Total"] + 
        1, rev(stt_q0.75_all[, "Total"] + 1)), col = "dark grey", border = NA)
    lines(stt_average_all[, "Time"], stt_average_all[, "Total"] + 1, lwd = 2)
    lines(stt_average_all[, "Time"], stt_average_all[, "nI"] + 1, lwd = 2, col = "cyan3")
    lines(stt_average_all[, "Time"], stt_average_all[, "Endemic"] + 1, lwd = 2, col = "dodgerblue1")

    legend(time, max(stt_q0.975_all), c("Total", "Non-endemic", "Endemic"), lty = 1, lwd = 2, col = c("black",
        "cyan3", "dodgerblue1"), cex = 1.2, border = NA, bty = "n")

    suppressWarnings(plot(NULL, NULL, xlim = rev(c(0, time)), ylim = c(1, max(stt_q0.975_type1)), 
        ylab = "No of species + 1", bty = "l", xaxs = "i", xlab = "Time before present", main = "STT Type 1 species", 
        log = "y", cex.lab = 1.5, cex.main = 1.2, cex.axis = 1.2))
    polygon(c(stt_average_type1[, "Time"], rev(stt_average_type1[, "Time"])), c(stt_q0.025_type1[, 
        "Total"] + 1, rev(stt_q0.975_type1[, "Total"] + 1)), col = "light grey", border = NA)
    polygon(c(stt_average_type1[, "Time"], rev(stt_average_type1[, "Time"])), c(stt_q0.25_type1[, 
        "Total"] + 1, rev(stt_q0.75_type1[, "Total"] + 1)), col = "dark grey", border = NA)
    lines(stt_average_type1[, "Time"], stt_average_type1[, "Total"] + 1, lwd = 2)
    lines(stt_average_type1[, "Time"], stt_average_type1[, "nI"] + 1, lwd = 2, col = "cyan3")
    lines(stt_average_type1[, "Time"], stt_average_type1[, "Endemic"] + 1, lwd = 2, col = "dodgerblue1")

    suppressWarnings(plot(NULL, NULL, xlim = rev(c(0, time)), ylim = c(1, max(stt_q0.975_type2)), 
        ylab = "No of species + 1", bty = "l", xaxs = "i", xlab = "Time before present", main = "STT Type 2 species", 
        log = "y", cex.lab = 1.5, cex.main = 1.2, cex.axis = 1.2))
    polygon(c(stt_average_type2[, "Time"], rev(stt_average_type2[, "Time"])), c(stt_q0.025_type2[, 
        "Total"] + 1, rev(stt_q0.975_type2[, "Total"] + 1)), col = "light grey", border = NA)
    polygon(c(stt_average_type2[, "Time"], rev(stt_average_type2[, "Time"])), c(stt_q0.25_type2[, 
        "Total"] + 1, rev(stt_q0.75_type2[, "Total"] + 1)), col = "dark grey", border = NA)
    lines(stt_average_type2[, "Time"], stt_average_type2[, "Total"] + 1, lwd = 2)
    lines(stt_average_type2[, "Time"], stt_average_type2[, "nI"] + 1, lwd = 2, col = "cyan3")
    lines(stt_average_type2[, "Time"], stt_average_type2[, "Endemic"] + 1, lwd = 2, col = "dodgerblue1")
        
    par(mfrow = c(1, 3))
    
    suppressWarnings(plot(NULL, NULL, xlim = rev(c(0, time)), ylim = c(1, max(stt_q0.975_all)), ylab = "No of species + 1", 
                          bty = "l", xaxs = "i", xlab = "Time before present", main = "Species-through-time - All species", 
                          log = "y", cex.lab = 1.5, cex.main = 1.2, cex.axis = 1.2))
    polygon(c(stt_average_all[, "Time"], rev(stt_average_all[, "Time"])), c(stt_q0.025_all[, "Total"] + 
                                                                              1, rev(stt_q0.975_all[, "Total"] + 1)), col = "light grey", border = NA)
    polygon(c(stt_average_all[, "Time"], rev(stt_average_all[, "Time"])), c(stt_q0.25_all[, "Total"] + 
                                                                              1, rev(stt_q0.75_all[, "Total"] + 1)), col = "dark grey", border = NA)
    lines(stt_average_all[, "Time"], stt_average_all[, "Total"] + 1, lwd = 2)
    lines(stt_average_all[, "Time"], stt_average_all[, "nI"] + 1, lwd = 2, col = "cyan3")
    lines(stt_average_all[, "Time"], stt_average_all[, "Endemic"] + 1, lwd = 2, col = "dodgerblue1")
    
    legend(time, max(stt_q0.975_all), c("Total", "Non-endemic", "Endemic"), lty = 1, lwd = 2, col = c("black", 
                                                                                                      "cyan3", "dodgerblue1"), cex = 1.2, border = NA, bty = "n")
    
    suppressWarnings(plot(NULL, NULL, xlim = rev(c(0, time)), ylim = c(1, max(stt_q0.975_type1)), 
                          ylab = "No of species + 1", bty = "l", xaxs = "i", xlab = "Time before present", main = "STT Type 1 species", 
                          log = "y", cex.lab = 1.5, cex.main = 1.2, cex.axis = 1.2))
    polygon(c(stt_average_type1[, "Time"], rev(stt_average_type1[, "Time"])), c(stt_q0.025_type1[, 
                                                                                                 "Total"] + 1, rev(stt_q0.975_type1[, "Total"] + 1)), col = "light grey", border = NA)
    polygon(c(stt_average_type1[, "Time"], rev(stt_average_type1[, "Time"])), c(stt_q0.25_type1[, 
                                                                                                "Total"] + 1, rev(stt_q0.75_type1[, "Total"] + 1)), col = "dark grey", border = NA)
    lines(stt_average_type1[, "Time"], stt_average_type1[, "Total"] + 1, lwd = 2)
    lines(stt_average_type1[, "Time"], stt_average_type1[, "nI"] + 1, lwd = 2, col = "cyan3")
    lines(stt_average_type1[, "Time"], stt_average_type1[, "Endemic"] + 1, lwd = 2, col = "dodgerblue1")
    
    suppressWarnings(plot(NULL, NULL, xlim = rev(c(0, time)), ylim = c(1, max(stt_q0.975_type2)), 
                          ylab = "No of species + 1", bty = "l", xaxs = "i", xlab = "Time before present", main = "STT Type 2 species", 
                          log = "y", cex.lab = 1.5, cex.main = 1.2, cex.axis = 1.2))
    polygon(c(stt_average_type2[, "Time"], rev(stt_average_type2[, "Time"])), c(stt_q0.025_type2[, 
                                                                                                 "Total"] + 1, rev(stt_q0.975_type2[, "Total"] + 1)), col = "light grey", border = NA)
    polygon(c(stt_average_type2[, "Time"], rev(stt_average_type2[, "Time"])), c(stt_q0.25_type2[, 
                                                                                                "Total"] + 1, rev(stt_q0.75_type2[, "Total"] + 1)), col = "dark grey", border = NA)
    lines(stt_average_type2[, "Time"], stt_average_type2[, "Total"] + 1, lwd = 2)
    lines(stt_average_type2[, "Time"], stt_average_type2[, "nI"] + 1, lwd = 2, col = "cyan3")
    lines(stt_average_type2[, "Time"], stt_average_type2[, "Endemic"] + 1, lwd = 2, col = "dodgerblue1")
    
  } else {
    # Default behavior to open a new device, which hurts vignettes
    if (use_dev_new == TRUE) {
      dev.new(width = 6, height = 6)
    }

    par(mfrow = c(1, 1))

    DAISIE_plot_stt(
      plot_plus_one = plot_plus_one,
      time = time,
      stt_q0.025_all = stt_q0.025_all,
      stt_q0.25_all = stt_q0.25_all,
      stt_average_all = stt_average_all,
      stt_q0.75_all = stt_q0.75_all,
      stt_q0.975_all = stt_q0.975_all
    )
   
  }
  
}