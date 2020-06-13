#' Remove superfluous lines in empty simulations
#'
#' @param sims
#'
#' @return
#' @keywords Internal
#'
#' @examples
remove_empty_elements <- function(sims) {
  for (i in seq_along(sims)) {
    if (sum(sims[[i]][[1]]$stt_all[nrow(sims[[i]][[1]]$stt_all), ]) == 0) {
      sims[[i]][[2]] <- NULL
    }
  }
  return(sims)
}
