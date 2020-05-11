#' Creates the list object for CS_version argument in DAISIE_ML_CS
#'
#' @param model the CS model to run, options are \code{"single"},
#' \code{"multi"}, or \code{"test}
#' @param pick_parameter the parameter to relax (integrate over). Options are
#' \code{"cladogenesis"}, \code{"extinction"}, \code{"carrying_capacity"},
#' \code{"immigration"}, or \code{"anagenesis"}
#' @param distribution the distribution to weigh the likelihood, either
#' \code{"lognormal"} or \code{"gamma"}
#' @param sd standard deviation of the distribution
#' @param num_cores number of cores to use in the integration
#'
#' @return A list of five elements
#' @export
create_CS_version <- function(model = "single",
                              pick_parameter = NULL,
                              distribution = NULL,
                              sd = NULL,
                              num_cores = NULL) {
  CS_version <- list(model = model,
                     pick_parameter = pick_parameter,
                     distribution = distribution,
                     sd = sd,
                     num_cores = num_cores)
  return(CS_version)
}
