#' @name DAISIE_DE_logpEC_general_cpp
#' @title Function to calculate the likelihood of observing an endemic lineage
#' with fixed colonization time. This is valid for infinite K according to the
#' DE equations.
#' @description This function calculates the log-likelihood of observing an
#' endemic lineage with fixed colonization time. This is valid for infinite K
#' according to the DE equations.
#' @inheritParams default_params_doc
#' @return the loglikelihood
#' @examples
#'
#' # Select a dataset from a DAISIE package
#'
#' data(Galapagos_datalist)
#' datalist <- Galapagos_datalist
#' brts <- datalist[[4]]$branching_times
#' missnumspec <- datalist[[4]]$missing_species
#'
#' # Define example parameters
#' pars1 <- c(0.2, 0.1, 0.05, 0.02, 0.03)
#'
#' # choose the method to solve the system of differential equations
#' log_likelihood <- DAISIE_DE_logpEC_general(brts = brts,
#'                                    missnumspec = missnumspec,
#'                                    pars1 = pars1,
#'                                    methode = "lsodes",
#'                                    reltolint = 1e-16,
#'                                    abstolint = 1e-16)
#' @noRd

DAISIE_DE_logpEC_general_cpp <- function(brts,
                                         missnumspec = 0,
                                         stac = 0,
                                         pars1,
                                         methode = "odeint::runge_kutta_cash_karp54",
                                         reltolint = 1e-15,
                                         abstolint = 1e-15) {

  if (!(stac %in% c(2, 3, 6))) {
    stop("stac must be 2, 3, or 6 for this function.")
  }

  lambda_c <- parameter[[1]]
  mu      <- parameter[[2]]
  gamma   <- parameter[[4]]
  lambda_a <- parameter[[5]]

  res <- .Call("DAISIE_DE_logpEC_general_rcpp",
               brts,
               missnumspec,
               lambda_c,
               lambda_a,
               mu,
               gamma,
               stac,
               methode,
               reltolint,
               abstolint)

  return(res)
}
