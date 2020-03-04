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
  distr <- 'gamma'
  if(dist == 'lognormal') {
    rho <- function(DAISIE_par, DAISIE_par_dist_pars) {
      return(dlnorm(x = DAISIE_par,meanlog = log(DAISIE_par_dist_pars[1]) - DAISIE_par_dist_pars[2]^2/2, sdlog = DAISIE_par_dist_pars[2]))
    }
    k_or_sigma <- 0.1
  } else {
    rho <- function(DAISIE_par, DAISIE_par_dist_pars) {
      return(dgamma(x = DAISIE_par, shape = DAISIE_par_dist_pars[2], scale = DAISIE_par_dist_pars[1]/DAISIE_par_dist_pars[2]))
    }
    # if CS_version = 2, then exponential distribution
    # if CS_version = 3, then gamma distribution with shape 2
    k_or_sigma <- floor(CS_version - 1)
  }

  if (DAISIE_version == 2) {
    DAISIE_version <- 2.3
  }
  pick <- 10 * (CS_version - floor(CS_version))
  mean_par <- pars1[pick]
  pars1new <- pars1

  DAISIE_loglik_par <- function(DAISIE_par) {
    pars1new[pick] <- DAISIE_par
    loglik_DAISIE_par <- exp(DAISIE_loglik(
      pars1 = pars1new,
      pars2 = pars2,
      brts = brts,
      stac = stac,
      missnumspec = missnumspec,
      methode = methode,
      abstolint = abstolint,
      reltolint = reltolint,
      verbose = verbose)) * rho(DAISIE_par = DAISIE_par, DAISIE_par_dist_pars = c(mean_par,k_or_sigma))
    return(loglik_DAISIE_par)
  }

  int_loglik <- log(integrate(f = DAISIE_loglik_par,
            lower = 0,
            upper = Inf
            ))
  return(int_loglik)
}
