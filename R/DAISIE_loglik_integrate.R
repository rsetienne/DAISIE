DAISIE_loglik_integrate <- function(
  pars1,
  pars2,
  brts,
  stac,
  missnumspec,
  CS_version,
  methode,
  abstolint,
  reltolint,
  verbose
) {

  rho <- function(K, K_dist_pars) {
    return(dgamma(x = K, shape = K_dist_pars[2], scale = K_dist_pars[1]))
  }

  DAISIE_loglik_K <- function(K) {
    loglik_K <- exp(DAISIE_loglik(
      pars1 = c(pars1[1:2], K, pars1[4:5]),
      pars2 = pars2,
      brts = brts,
      stac = stac,
      missnumspec = missnumspec,
      methode = methode,
      abstolint = abstolint,
      reltolint = reltolint,
      verbose = FALSE)) * rho(K = K, K_dist_pars = c(pars1[3],CS_version - 1))
    # if CS_version = 2, then exponential distribution
    # if CS_version = 3, then gamma distribution with shape 2
    return(loglik_K)
  }

  int_loglik <- log(integrate(f = DAISIE_loglik_K,
            lower = 0,
            upper = Inf
            ))
  return(int_loglik)
}
