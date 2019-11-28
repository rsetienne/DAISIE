DAISIE_loglik_integrate <- function(
  pars1 = pars1,
  pars2 = pars2,
  brts = brts,
  stac = stac,
  missnumspec = missnumspec,
  methode = methode,
  abstolint = abstolint,
  reltolint = reltolint,
  verbose = FALSE
) {

  rho <- function(K, K_dist_pars) {
    exp(- K / K_dist_pars[1])
    #can add gamma
    #extra parameter
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
      verbose = FALSE)) * rho(K = K, K_dist_pars = pars1[3])
    return(loglik_K)
  }

  int_loglik <- log(integrate(f = DAISIE_loglik_K,
            lower = 0,
            upper = Inf
            ))
  return(int_loglik)
}
