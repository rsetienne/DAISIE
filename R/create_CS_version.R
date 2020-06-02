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
#' @param multi_rate_optim_method A string describing the optimization method
#'  used to find the optimum of the integrand, either \code{"optimize"} or
#'  \code{"subplex"} or \code{"Nelder-Mead"}.
#' @return A list of five elements
#' \itemize{
#'   \item{model: the CS model to run, options are \code{1} for single rate
#'   DAISIE model, \code{2} for multi-rate DAISIE, or \code{0} for IW test
#'   model}
#'   \item{pick_parameter: the parameter to relax (integrate over). Options are
#' \code{"cladogenesis"}, \code{"extinction"}, \code{"carrying_capacity"},
#' \code{"immigration"}, or \code{"anagenesis"}}
#'   \item{distribution: the distribution to weigh the likelihood, either
#' \code{"lognormal"} or \code{"gamma"}}
#'   \item{sd: standard deviation of the distribution}
#'   \item{optimmethod: the method used to find the optimum of the integrand,
#'   either \code{"optimize"} or \code{"subplex"} or \code{"Nelder-Mead"}}
#' }
#' @export
create_CS_version <- function(model = 1,
                              pick_parameter = NULL,
                              distribution = NULL,
                              sd = NULL,
                              multi_rate_optim_method = NULL) {

  if (model != 1 && model != 2 && model != 3) {
    stop("model must be either 1, 2 or 3")
  }
  if (model == 2 && is.null(pick_parameter) ||
      model == 2 && is.null(distribution) ||
      model == 2 && is.null(sd) ||
      model == 2 && is.null(multi_rate_optim_method)) {
    stop("pick_parameter, distribution, sd and multi_rate_optim_method
         required for multi-rate model")
  }
  CS_version <- list(model = model,
                     pick_parameter = pick_parameter,
                     distribution = distribution,
                     sd = sd,
                     multi_rate_optim_method = multi_rate_optim_method)
  return(CS_version)
}
