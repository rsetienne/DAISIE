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
  # DaVinci code for CS_version:
  # The sign determines whether the lognormal (negative) or the gamma distribution (positive)
  # should be chosen.
  # The number before the decimal separator determines what DAISIE parameter should be chosen where
  # counting starts at 2. For example, CS_version = 2.x corresponds to the 1st parameter,
  # lambda^c, whereas CS_version = 4.x corresponds to the 3rd parameter, K.
  # Then the number after the decimal separator divided by 10 sets the standard deviation for the
  # lognormal and it sets the shape for the gamma.
  # For example CS_version = -x.34 sets the sd of the lognormal to 3.4. The sd is thus limited to
  # values below 10.
  # A value of CS_version = +4.31 thus indicates a gamma distribution (+) for the 3rd parameter (4 - 1)
  # with a shape parameter of 3.1.
  # For most purposes values of -4.01 (lognormal on K, with sd of 0.1) or +4.2 (exponential) or +4.3
  # (gamma with shape 2) are sufficient.
  if(sgn(CS_version) == -1) {
    distr <- 'lognormal'
    rho <- function(DAISIE_par, DAISIE_par_dist_pars) {
      return(dlnorm(x = DAISIE_par,
                    meanlog = log(DAISIE_par_dist_pars[1]) - (DAISIE_par_dist_pars[2]^2)/2,
                    sdlog = sqrt(log(1 + DAISIE_par_dist_pars[2] / DAISIE_par_dist_pars[1]^2))))
    }
  } else {
    distr <- 'gamma'
    rho <- function(DAISIE_par, DAISIE_par_dist_pars) {
      return(dgamma(x = DAISIE_par,
                    shape = DAISIE_par_dist_pars[2]^2,
                    scale = DAISIE_par_dist_pars[1]/DAISIE_par_dist_pars[2]^2))
    }
    if (CS_version == 2) {
      CS_version <- 2.3
    }
    # if CS_version = x.2, then exponential distribution
    # if CS_version = x.3, then gamma distribution with shape 2
  }

  k_or_sigma <- 10 * (abs(CS_version) - floor(abs(CS_version)))
  pick <- floor(abs(CS_version)) - 1
  mean_par <- pars1[pick]
  pars1new <- pars1

  DAISIE_loglik_integrand <- function(DAISIE_par) {
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
      verbose = verbose)
      ) *
      rho(
        DAISIE_par = DAISIE_par,
        DAISIE_par_dist_pars = c(mean_par,k_or_sigma)
      )
    return(loglik_DAISIE_par)
  }

  integrated_loglik <- log(integrate(f = DAISIE_loglik_integrand,
                                     lower = 0,
                                     upper = Inf
  ))
  return(integrated_loglik)
}
