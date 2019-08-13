#' Extract the STT median from the output of \code{\link{DAISIE_sim}}
#' @param island_replicates the result of \code{\link{DAISIE_sim}}
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
    sum_endemics <- island_replicates[[x]][[1]]$stt_all[, "nA0"] + 
      island_replicates[[x]][[1]]$stt_all[, "nC0"] + 
      island_replicates[[x]][[1]]$stt_all[,"nA1"] +
      island_replicates[[x]][[1]]$stt_all[,"nC1"]
    total <- island_replicates[[x]][[1]]$stt_all[, "nA0"] + 
             island_replicates[[x]][[1]]$stt_all[, "nC0"] + 
             island_replicates[[x]][[1]]$stt_all[, "nI0"] +
             island_replicates[[x]][[1]]$stt_all[, "nA1"] + 
             island_replicates[[x]][[1]]$stt_all[, "nC1"] + 
             island_replicates[[x]][[1]]$stt_all[, "nI1"]
    
    nItotal<-island_replicates[[x]][[1]]$stt_all[,"nI0"] +
             island_replicates[[x]][[1]]$stt_all[,"nI1"]
    nAtotal<-island_replicates[[x]][[1]]$stt_all[,"nA0"] +
             island_replicates[[x]][[1]]$stt_all[,"nA1"]
    nCtotal<-island_replicates[[x]][[1]]$stt_all[,"nC0"] +
             island_replicates[[x]][[1]]$stt_all[,"nC1"]
      
    complete_arr[, , x] <- cbind(
      island_replicates[[x]][[1]]$stt_all[, "Time"],
      nItotal,
      nAtotal,
      nCtotal,
      sum_endemics, 
      total
    )
  }
  stt_average_all <- apply(complete_arr, c(1, 2), stats::median)
  colnames(stt_average_all) <- c("Time", "nItotal", "nAtotal", "nCtotal", "Endemic", "Total")
  stt_average_all
}