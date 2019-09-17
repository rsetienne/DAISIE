#' Extract the STT median from the output of \code{\link{DAISIE_sim}}
#' @param island_replicates the result of \code{\link{DAISIE_sim}}
#' @param Tpars A named list containing diversification rates considering two trait states:
#' \itemize{
#'   \item{[1]:A numeric with the per capita transition rate with state1}
#'   \item{[2]:A numeric with the per capita immigration rate with state2}
#'   \item{[3]:A numeric with the per capita extinction rate with state2}
#'   \item{[4]:A numeric with the per capita anagenesis rate with state2}
#'   \item{[5]:A numeric with the per capita cladogenesis rate with state2}
#'   \item{[6]:A numeric with the per capita transition rate with state2} 
#'   \item{[7]:A numeric with the number of species with trait state 2 on mainland} 
#' }
#' @return a matrix (?)
#' @export
DAISIE_extract_stt_median <- function(
  island_replicates,
  Tpars = NULL
) {
  replicates <- length(island_replicates)
  time <- max(island_replicates[[1]][[1]]$stt_all[, 1])
  
  ### STT ALL species
  s_freq <- length(island_replicates[[1]][[1]]$stt_all[, 1])
  complete_arr <- array(dim = c(s_freq, 6, replicates))
  
  if(is.null(Tpars)){
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
  }else{
    for (x in 1:replicates) {
      sum_endemics <- island_replicates[[x]][[1]]$stt_all[, "nA"] + 
        island_replicates[[x]][[1]]$stt_all[, "nC"] + 
        island_replicates[[x]][[1]]$stt_all[,"nA2"] +
        island_replicates[[x]][[1]]$stt_all[,"nC2"]
      total <- island_replicates[[x]][[1]]$stt_all[, "nA"] + 
        island_replicates[[x]][[1]]$stt_all[, "nC"] + 
        island_replicates[[x]][[1]]$stt_all[, "nI"] +
        island_replicates[[x]][[1]]$stt_all[, "nA2"] + 
        island_replicates[[x]][[1]]$stt_all[, "nC2"] + 
        island_replicates[[x]][[1]]$stt_all[, "nI2"]
      
      nI<-island_replicates[[x]][[1]]$stt_all[,"nI"] +
        island_replicates[[x]][[1]]$stt_all[,"nI2"]
      nA<-island_replicates[[x]][[1]]$stt_all[,"nA"] +
        island_replicates[[x]][[1]]$stt_all[,"nA2"]
      nC<-island_replicates[[x]][[1]]$stt_all[,"nC"] +
        island_replicates[[x]][[1]]$stt_all[,"nC2"]
      
      complete_arr[, , x] <- cbind(
        island_replicates[[x]][[1]]$stt_all[, "Time"],
        nI,
        nA,
        nC,
        sum_endemics, 
        total
      )
    }
  }
  stt_average_all <- apply(complete_arr, c(1, 2), stats::median)
  colnames(stt_average_all) <- c("Time", "nI", "nA", "nC", "Endemic", "Total")
  return(stt_average_all)
}