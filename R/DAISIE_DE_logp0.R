#' @name DAISIE_DE_logp0
#' @title Log-likelihood of a mainland species that colonizes the island
#' but leaves no descendants. This is valid for infinite K according to the DE
#' equations.
#' @description
#' Computes the log-likelihood of a colonization event where a mainland species
#' arrives on the island but does not leave any surviving descendants. This is
#' valid for infinite K according to the DE equations.
#' @inheritParams default_params_doc
#' @return the loglikelihood
#' @examples
#' # Example model parameters
#' pars1 <- c(0.2, 0.1, 0.05, 0.02, 0.03)
#'
#' # Compute log-likelihood
#' log_likelihood <- DAISIE_DE_logp0(island_age = 10,
#'                                  pars1 = pars1,
#'                                  methode = "lsodes",
#'                                  reltolint = 1E-12,
#'                                  abstolint = 1E-12)
#' @noRd

DAISIE_DE_logp0 <- function(island_age,
                            pars1,
                            methode,
                            reltolint,
                            abstolint,
                            use_rcpp = FALSE) {
  t0 <- island_age
  tp <- 0

  # Set initial conditions
  initial_conditions0 <- c(DA1 = 1, Dm1 = 0, E = 0)

  # Time sequence for interval [t0, tp]
  time0 <- c(tp, t0)

  solution0 <- DAISIE_DE_solve_branch(interval_func = interval3,
                                      initial_conditions = initial_conditions0,
                                      parameter = pars1,
                                      time = time0,
                                      methode = methode,
                                      atol = abstolint,
                                      rtol = reltolint,
                                      use_rcpp = use_rcpp)

  # Extract log-likelihood
  L0 <- solution0[, "DA1"][[2]]
  logL0b <- log(L0)
  return(logL0b)
}

