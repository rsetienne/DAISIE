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
  rho <- function(DAISIE_par, DAISIE_par_dist_pars) {
    return(dgamma(x = DAISIE_par, shape = DAISIE_par_dist_pars[2], scale = DAISIE_par_dist_pars[1]))
  }

  DAISIE_loglik_par <- function(DAISIE_par) {
    if (DAISIE_version == 2) {
      DAISIE_version <- 2.3
    }
    pick <- 10 * (CS_version - floor(CS_version))
    pars1new <- pars1
    pars1[pick] <- DAISIE_par
    loglik_DAISIE_par <- exp(DAISIE_loglik(
      pars1 = pars1new,
      pars2 = pars2,
      brts = brts,
      stac = stac,
      missnumspec = missnumspec,
      methode = methode,
      abstolint = abstolint,
      reltolint = reltolint,
      verbose = verbose)) * rho(DAISIE_par = DAISIE_par, DAISIE_par_dist_pars = c(pars1[pick],floor(CS_version - 1)))
    # if CS_version = 2, then exponential distribution
    # if CS_version = 3, then gamma distribution with shape 2
    return(loglik_DAISIE_par)
  }

  int_loglik <- log(integrate(f = DAISIE_loglik_K,
            lower = 0,
            upper = Inf
            ))
  return(int_loglik)
}
