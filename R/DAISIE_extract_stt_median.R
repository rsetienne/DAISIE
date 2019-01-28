#' Extract the STT median from the output of \code{\link{DAISIE_sim}}
#' @param island_replicates the result of \code{\link{DAISIE_sim}}
#' @return a matrix (?) 
DAISIE_extract_stt_median <- function(
  island_replicates
) {
  replicates <- length(island_replicates)
  time <- max(island_replicates[[1]][[1]]$stt_all[, 1])
  
  ### STT ALL species
  s_freq <- length(island_replicates[[1]][[1]]$stt_all[, 1])
  complete_arr <- array(dim = c(s_freq, 6, replicates))
  
  for (x in 1:replicates) {
    sum_endemics <- island_replicates[[x]][[1]]$stt_all[, "nA"] + 
      island_replicates[[x]][[1]]$stt_all[, "nC"]
    total <- island_replicates[[x]][[1]]$stt_all[, "nA"] + 
             island_replicates[[x]][[1]]$stt_all[, "nC"] + 
             island_replicates[[x]][[1]]$stt_all[, "nI"]
    complete_arr[, , x] <- cbind(
      island_replicates[[x]][[1]]$stt_all[, c("Time", "nI", "nA", "nC")], 
      sum_endemics, 
      total
    )
  }
  stt_average_all <- apply(complete_arr, c(1, 2), median)
  colnames(stt_average_all) <- c("Time", "nI", "nA", "nC", "Endemic", "Total")
  stt_average_all
}