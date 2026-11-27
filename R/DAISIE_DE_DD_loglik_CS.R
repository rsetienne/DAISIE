DAISIE_DE_DD_loglik <- function(pars1,
                                pars2,
                                brts,
                                missnumspec,
                                stac,
                                methode = "odeint::runge_kutta_cash_karp54",
                                reltolint = 1e-10,
                                abstolint = 1e-16,
                                sampling = 'n',
                                dd_recolonization = FALSE,
                                mu_NE = pars1[2]) {
  # An endemic clade with a maximum colonization time and only one sampled
  # species (and possibly missing species) is an endemic singleton with a
  # maximum colonization time (as in DAISIE_DE_loglik).
  if (stac == 6 && length(brts) == 2) {
    stac <- 5
  }
  # with linear diversity dependence, the clade can not be larger than K
  if (pars2[2] %in% c(1, 11)) {
    S <- (length(brts) - 1) * !(stac %in% c(0, 1, 4, 8)) + (stac %in% c(1, 3, 4, 7, 8))
    if (ceiling(pars1[3]) < S + missnumspec) {
      return(-Inf)
    }
  }
  if (stac == 0) {
    loglikelihood <- DAISIE_DE_DD_logp0(island_age = brts[1],
                                        pars1 = pars1,
                                        pars2 = pars2,
                                        methode = methode,
                                        reltolint = reltolint,
                                        abstolint = abstolint,
                                        dd_recolonization = dd_recolonization,
                                        mu_NE = mu_NE)
  } else if (stac == 1 || stac == 4 || stac == 8) {
    loglikelihood <- DAISIE_DE_DD_logpNE(brts = brts,
                                         missnumspec = missnumspec,
                                         pars1 = pars1,
                                         pars2 = pars2,
                                         stac = stac,
                                         methode = methode,
                                         reltolint = reltolint,
                                         abstolint = abstolint,
                                         dd_recolonization = dd_recolonization,
                                         mu_NE = mu_NE)
  } else if (((stac == 2 || stac == 3 || stac == 5 || stac == 7) && length(brts) == 2) || stac == 9) {
    loglikelihood <- DAISIE_DE_DD_logpES(brts = brts,
                                         missnumspec = missnumspec,
                                         pars1 = pars1,
                                         pars2 = pars2,
                                         stac = stac,
                                         methode = methode,
                                         reltolint = reltolint,
                                         abstolint = abstolint,
                                         sampling = sampling,
                                         dd_recolonization = dd_recolonization,
                                         mu_NE = mu_NE)
  } else if (((stac == 2 || stac == 3 || stac == 5 || stac == 7) && length(brts) > 2) || stac == 6) {
    loglikelihood <- DAISIE_DE_DD_logpEC(brts = brts,
                                         missnumspec = missnumspec,
                                         pars1 = pars1,
                                         pars2 = pars2,
                                         stac = stac,
                                         methode = methode,
                                         reltolint = reltolint,
                                         abstolint = abstolint,
                                         sampling = sampling,
                                         dd_recolonization = dd_recolonization,
                                         mu_NE = mu_NE)
  } else {
    stop("Unknown stac value: ", stac)
  }
  return(loglikelihood)
}

#' @name DAISIE_DE_DD_loglik_CS
#' @title Computes the loglikelihood of the diversity-dependent DAISIE_DE
#' model given data and a set of model parameters
#' @description Diversity-dependent version of \code{DAISIE_DE_loglik_CS}.
#' Each DAISIE_DE variable carries the number of hidden species (species that
#' are not in the data) alive at that time, so the rates can depend on the
#' total number of species of the clade on the island. This gives the same
#' likelihood as \code{DAISIE_loglik_CS} with finite K.
#' @param pars1 Contains the model parameters: \cr \cr
#' \code{pars1[1]} corresponds to lambda^c (cladogenesis rate) \cr
#' \code{pars1[2]} corresponds to mu (extinction rate) \cr
#' \code{pars1[3]} corresponds to K (clade-level carrying capacity) \cr
#' \code{pars1[4]} corresponds to gamma (immigration rate) \cr
#' \code{pars1[5]} corresponds to lambda^a (anagenesis rate) \cr
#' @param pars2 As in \code{DAISIE_loglik_CS}: \cr \cr
#' \code{pars2[1]} largest number of hidden species kept when the clade size
#' is not limited by K (as lx in DAISIE) \cr
#' \code{pars2[2]} ddmodel: 0, 1, 11, 2 or 21, as in DAISIE \cr
#' \code{pars2[3]} cond: 1 to condition on island age and non-extinction of
#' the island biota, 0 for no conditioning \cr
#' \code{pars2[4]} verbose \cr
#' @param datalist Data object, as in \code{DAISIE_DE_loglik_CS}
#' @param methode Method of the ODE-solver: an "odeint::" stepper (default
#' "odeint::runge_kutta_cash_karp54") or a deSolve method such as "lsoda"
#' @param abstolint Absolute tolerance of the integration
#' @param reltolint Relative tolerance of the integration
#' @param sampling 'n' (default): missing_species is the exact number of
#' missing species; 'rho': each species is sampled with probability
#' S / (S + missing_species), as in \code{DAISIE_DE_loglik_CS}
#' @param dd_recolonization FALSE (default) or TRUE: is recolonization
#' diversity dependent? As in \code{create_CS_version}
#' @param mu_NE Extinction rate of the non-endemic species; NULL (default)
#' means equal to mu, as in DAISIE
#' @return The loglikelihood
#' @export DAISIE_DE_DD_loglik_CS
DAISIE_DE_DD_loglik_CS <- function(pars1,
                                   pars2 = c(100, 11, 0, 0),
                                   datalist,
                                   methode = "odeint::runge_kutta_cash_karp54",
                                   abstolint = 1e-16,
                                   reltolint = 1e-10,
                                   sampling = 'n',
                                   dd_recolonization = FALSE,
                                   mu_NE = NULL) {
  if (length(pars1) != 5) {
    stop("DAISIE_DE_DD_loglik_CS: only the model with one set of 5 parameters is implemented.")
  }
  if (is.null(mu_NE)) {
    mu_NE <- pars1[2]
  }
  if (min(pars1) < 0) {
    return(-Inf)
  }
  cond <- pars2[3]
  island_age <- datalist[[1]]$island_age
  logp0 <- DAISIE_DE_DD_logp0(island_age = island_age,
                              pars1 = pars1,
                              pars2 = pars2,
                              methode = methode,
                              reltolint = reltolint,
                              abstolint = abstolint,
                              dd_recolonization = dd_recolonization,
                              mu_NE = mu_NE)
  if (is.null(datalist[[1]]$not_present)) {
    loglik <- (datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2) * logp0
    numimm <- (datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2) + length(datalist) - 1
  } else {
    loglik <- datalist[[1]]$not_present * logp0
    numimm <- datalist[[1]]$not_present + length(datalist) - 1
  }
  logcond <- (cond == 1) * log(1 - exp(numimm * logp0))
  loglik <- loglik - logcond

  vec_loglikelihood <- c()
  for (i in seq_len(length(datalist) - 1) + 1) {
    stac <- datalist[[i]]$stac
    brts <- datalist[[i]]$branching_times
    missnumspec <- datalist[[i]]$missing_species
    loglikelihood <- DAISIE_DE_DD_loglik(pars1 = pars1,
                                         pars2 = pars2,
                                         brts = brts,
                                         missnumspec = missnumspec,
                                         stac = stac,
                                         methode = methode,
                                         reltolint = reltolint,
                                         abstolint = abstolint,
                                         sampling = sampling,
                                         dd_recolonization = dd_recolonization,
                                         mu_NE = mu_NE)
    vec_loglikelihood <- c(vec_loglikelihood, loglikelihood)

    print_parameters_and_loglik(
      pars = c(stac, pars1),
      loglik = loglikelihood,
      verbose = pars2[4],
      parnames = c("lambda^c", "mu", "K", "gamma", "lambda^a"),
      type = 'clade_loglik'
    )
  }
  loglik <- sum(vec_loglikelihood) + loglik
  return(loglik)
}
