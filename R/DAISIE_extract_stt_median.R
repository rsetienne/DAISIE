#' Extract the STT median from the output of DAISIE_sim functions
#'
#' @inheritParams default_params_doc
#'
#' @return a matrix (?)
#' @export
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
      total)
  }
  stt_average_all <- apply(complete_arr, c(1, 2), stats::median)
  colnames(stt_average_all) <- c("Time", "nI", "nA", "nC", "Endemic", "Total")
  return(stt_average_all)
}
