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
#' brts <- datalist[[5]]$branching_times
#' missnumspec <- datalist[[5]]$missing_species
#'
#' # Define example parameters
#' pars1 <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)
#'
#' # choose the method to solve the system of differential equations
#' log_likelihood <- DAISIE_DE_logpEC_general_cpp(brts = brts,
#'                                    stac = 2,
#'                                    missnumspec = missnumspec,
#'                                    pars1 = pars1,
#'                                    methode = "odeint::runge_kutta_cash_karp54",
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

  lambda_c <- pars1[[1]]
  mu      <- pars1[[2]]
  gamma   <- pars1[[4]]
  lambda_a <- pars1[[5]]

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
