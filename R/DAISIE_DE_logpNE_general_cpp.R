#' @name DAISIE_DE_logpNE_general_cpp
#' @title Function to calculate the likelihood of observing a non-endemic lineage
#' with fixed colonization time. This is valid for infinite K according to the DE
#' equations.
#' @description This function calculates the log-likelihood of observing a non-endemic lineage
#' with fixed colonization time. This is valid for infinite K according to the DE
#' equations.
#' @inheritParams default_params_doc
#' @return the loglikelihood
#' @examples
#'
#' # Select a dataset from a DAISIE package
#'
#' data(Galapagos_datalist)
#' datalist <- Galapagos_datalist
#' brts <- datalist[[3]]$branching_times
#' # Define example parameters
#' pars1 <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)
#'
#' # choose the method to solve the system of differential equations
#' log_likelihood <- DAISIE_DE_logpNE_general(brts = brts,
#'                                    stac = 4,
#'                                    pars1 = pars1,
#'                                    methode = "lsodes",
#'                                    reltolint = 1e-16,
#'                                    abstolint = 1e-16,
#'                                    use_rcpp = TRUE)
#' @noRd
DAISIE_DE_logpNE_general_cpp <- function(brts,
                                     pars1,
                                     stac = 4,
                                     methode = "odeint::runge_kutta_cash_karp54",
                                     reltolint,
                                     abstolint,
                                     use_rcpp = FALSE) {

  if (!(stac %in% c(1, 4, 8))) {
    stop("NE only supports stac values of 1, 4 and 8")
  }

  lambda_c <- pars1[[1]]
  mu      <- pars1[[2]]
  gamma   <- pars1[[4]]
  lambda_a <- pars1[[5]]

  res <- .Call("DAISIE_DE_logpNE_general_rcpp",
               brts,
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
