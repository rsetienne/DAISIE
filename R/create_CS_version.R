#' Creates the list object for CS_version argument in DAISIE_ML_CS
#'
#' @param model the CS model to run, options are \code{1} for single rate
#' DAISIE model, \code{2} for multi-rate DAISIE, or \code{0} for IW test
#' model
#' @param pick_parameter the parameter to relax (integrate over). Options are
#' \code{"cladogenesis"}, \code{"extinction"}, \code{"carrying_capacity"},
#' \code{"immigration"}, or \code{"anagenesis"}
#' @param distribution the distribution to weigh the likelihood, either
#' \code{"lognormal"} or \code{"gamma"}
#' @param sd standard deviation of the distribution
#' @param number_of_cores number of cores to use in the integration
#'
#' @return A list of five elements
#' @export
create_CS_version <- function(model = 1,
                              pick_parameter = NULL,
                              distribution = NULL,
                              sd = NULL,
                              number_of_cores = NULL) {
  CS_version <- list(model = model,
                     pick_parameter = pick_parameter,
                     distribution = distribution,
                     sd = sd,
                     number_of_cores = number_of_cores)
  return(CS_version)
}
