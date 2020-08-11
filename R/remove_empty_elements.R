#' Remove superfluous lines in empty simulations
#'
#' @param sims List with output of DAISIE_sim
#'
#' @return Output of DAISIE_sim without empty elements
#' @keywords internal
remove_empty_elements <- function(sims) {
  for (i in seq_along(sims)) {
    if (sum(sims[[i]][[1]]$stt_all[nrow(sims[[i]][[1]]$stt_all), ]) == 0) {
      sims[[i]][[2]] <- NULL
    }
  }
  return(sims)
}
