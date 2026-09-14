DAISIE_DE_loglik <- function(pars1,
                             brts,
                             missnumspec,
                             stac,
                             methode,
                             reltolint = 1e-15,
                             abstolint = 1e-15,
                             sampling = 'rho') {
  S <- length(brts) - 1
  fac <- S * (log(S) - log(S + missnumspec))
  if( stac == 0) {
    loglikelihood <- DAISIE_DE_logp0(island_age = brts[1],
                                     pars1 = pars1,
                                     reltolint = abstolint,
                                     abstolint = reltolint,
                                     methode = methode)
  } else if (stac == 1 || stac == 4 || stac == 8) {
    loglikelihood <- DAISIE_DE_logpNE(brts = brts,
                                      pars1 = pars1,
                                      stac = stac,
                                      methode = methode,
                                      reltolint = reltolint,
                                      abstolint = abstolint)

  } else if (((stac == 2 || stac == 3 || stac == 5 || stac == 7) && length(brts) == 2) || stac == 9) {
    if(sampling == 'rho' || missnumspec == 0)
      loglikelihood <- DAISIE_DE_logpES(brts = brts,
                                        missnumspec = missnumspec,
                                        stac = stac,
                                        pars1 = pars1,
                                        methode = methode,
                                        reltolint = reltolint,
                                        abstolint = abstolint) + fac
    else
      loglikelihood <- DAISIE_DE_n(DAISIE_DE_function = DAISIE_DE_logpES,
                                   brts = brts,
                                   missnumspec = missnumspec,
                                   stac = stac,
                                   pars1 = pars1,
                                   methode = methode,
                                   reltolint = reltolint,
                                   abstolint = abstolint)
  } else if (((stac == 2 || stac == 3 || stac == 5 || stac == 7) && length(brts) > 2) || stac == 6) {
    if(sampling == 'rho' || missnumspec == 0)
      loglikelihood <- DAISIE_DE_logpEC(brts = brts,
                                        missnumspec = missnumspec,
                                        stac = stac,
                                        pars1 = pars1,
                                        methode = methode,
                                        reltolint = reltolint,
                                        abstolint = abstolint) + fac
    else
      loglikelihood <- DAISIE_DE_n(DAISIE_DE_function = DAISIE_DE_logpEC,
                                   brts = brts,
                                   missnumspec = missnumspec,
                                   stac = stac,
                                   pars1 = pars1,
                                   methode = methode,
                                   reltolint = reltolint,
                                   abstolint = abstolint)
  } else {
    stop("Unknown stac value: ", stac)
  }
  return(loglikelihood)
}

#' @name DAISIE_DE_loglik_CS
#' @title Computes the loglikelihood of the DAISIE_DE model given data and a set
#' of model parameters
#' @description Computes the loglikelihood of the DAISIE_DE model given
#' colonization and branching times for lineages on an island, and a set of model
#' parameters. The output is a loglikelihood value
#' @inheritParams default_params_doc
#' @param pars1 Contains the model parameters: \cr \cr
#' \code{pars1[1]} corresponds to lambda^c (cladogenesis rate) \cr
#' \code{pars1[2]} corresponds to mu (extinction rate of endemic species) \cr
#' \code{pars1[3]} corresponds to mu2 (extinction rate of non-endemic species) \cr
#' \code{pars1[4]} corresponds to gamma (immigration rate) \cr
#' \code{pars1[5]} corresponds to lambda^a (anagenesis rate) \cr
#' @param pars2 Contains the model settings \cr \cr
#' \code{pars2[1]} irrelevant for DAISIE_DE \cr
#' \code{pars2[2]} irrelevant for DAISIE_DE \cr
#' \code{pars2[3]} corresponds to cond = setting of conditioning\cr \cr
#' cond = 0 : conditioning on island age \cr
#' cond = 1 : conditioning on island age and non-extinction of the island biota \cr \cr
#' cond > 1 : conditioning on island age and having at least cond colonizations on the island \cr \cr
#' \code{pars2[4]} sets the level of verbosity. When equal to 0, no output is generated. At higher values
#' (1 or 2) more output will be generated.
#' @param datalist Data object containing information on colonisation and
#' branching times. This object can be generated using the DAISIE_dataprep
#' function, which converts a user-specified data table into a data object, but
#' the object can of course also be entered directly. It is an R list object
#' with the following elements.\cr The first element of the list has two or
#' three components:
#' \cr \cr \code{$island_age} - the island age \cr
#' Then, depending on whether a distinction between types is made, we have:\cr
#' \code{$not_present} - the number of mainland lineages that are not present
#' on the island \cr
#' or:\cr
#' \code{$not_present_type1} - the number of mainland lineages of type 1 that are not present on the island \cr
#' \code{$not_present_type2} - the number of mainland lineages of type 2 that
#' are not present on the island \cr \cr
#' The remaining elements of the list
#' each contains information on a single colonist lineage on the island and has
#' 5 components:\cr \cr
#' \code{$colonist_name} - the name of the species or
#' clade that colonized the island \cr
#' \code{$branching_times} - island age and
#' stem age of the population/species in the case of Non-endemic,
#' Non-endemic_MaxAge and Endemic anagenetic species. For cladogenetic species
#' these should be island age and branching times of the radiation including
#' the stem age of the radiation.\cr
#' \code{$stac} - the status of the colonist \cr \cr
#' - Non_endemic_MaxAge: 1 \cr
#' - Endemic: 2 \cr
#' - Endemic&Non_Endemic: 3 \cr
#' - Non_Endemic: 4 \cr
#' - Endemic_Singleton_MaxAge: 5 \cr
#' - Endemic_Clade_MaxAge: 6 \cr
#' - Endemic&Non_Endemic_Clade_MaxAge: 7 \cr
#' - Non_endemic_MaxAge_MinAge: 8 \cr
#' - Endemic_Singleton_MaxAge_MinAge: 9 \cr
#' \code{$missing_species} - number of island species that were not sampled for
#' particular clade (only applicable for endemic clades) \cr
#' \code{$type1or2} - whether the colonist belongs to type 1 or type 2. Currently
#' not implemented for DAISIE_DE \cr
#' @param methode Method of the ODE-solver. See package deSolve for details.
#' Default is "odeint::runge_kutta_cask_karp54"
#' @param abstolint Absolute tolerance of the integration
#' @param reltolint Relative tolerance of the integration
#' @return The loglikelihood
#' @author Rampal S. Etienne & Bart Haegeman
#' @seealso \code{\link{DAISIE_ML}}
#' @references O.N. Dehayem et al. 2026. Preprint.
#' @export DAISIE_DE_loglik_CS

DAISIE_DE_loglik_CS <- function( pars1,
                                 pars2,
                                 datalist,
                                 methode = "odeint::runge_kutta_cash_karp54",
                                 abstolint = 1e-15,
                                 reltolint = 1e-15,
                                 equal_extinction = TRUE,
                                 sampling = 'n') {

  # Apply equal extinction condition AFTER initializing pars1
  if (equal_extinction) {
    pars1[3] <- pars1[2]
  }
  cond <- pars2[3]
  island_age <- datalist[[1]]$island_age
  parameter <- pars1
  if (length(parameter) == 5) {
    logp0 <- DAISIE_DE_logp0(island_age = island_age,
                             pars1 = pars1,
                             reltolint = reltolint, #was 1E-12
                             abstolint = abstolint, #was 1E-12
                             methode = methode)
    if (is.null(datalist[[1]]$not_present)) {
      loglik <- (datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2) * logp0
      numimm <- (datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2) + length(datalist) - 1
    } else {
      loglik <- datalist[[1]]$not_present * logp0
      numimm <- datalist[[1]]$not_present + length(datalist) - 1
    }
    logcond <- (cond == 1) * log(1 - exp(numimm * logp0))
    for (i in 2:length(datalist)) {
      datalist[[i]]$type1or2 <- 1
    }
  } else {
    numimm <- datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2 + length(datalist) - 1
    numimm_type2 <- length(which(unlist(datalist)[which(names(unlist(datalist)) == "type1or2")] == 2))
    numimm_type1 <- length(datalist) - 1 - numimm_type2
    if (!is.na(parameter[11])) {
      if (parameter[11] < numimm_type2 / numimm | parameter[11] > (1 - numimm_type1 / numimm)) {
        return(-Inf)
      }
      datalist[[1]]$not_present_type2 <- max(0, round(parameter[11] * numimm) - numimm_type2)
      datalist[[1]]$not_present_type1 <- numimm - (length(datalist) - 1) - datalist[[1]]$not_present_type2
    }
    logp0_type1 <- DAISIE_DE_logp0(island_age = island_age,
                                   pars1 = pars1,
                                   methode = methode,
                                   reltolint = reltolint,
                                   abstolint = abstolint)
    logp0_type2 <- DAISIE_DE_logp0(island_age = island_age,
                                   pars1 = pars1,
                                   methode = methode,
                                   reltolint = reltolint,
                                   abstolint = abstolint)
    loglik <- datalist[[1]]$not_present_type1 * logp0_type1 + datalist[[1]]$not_present_type2 * logp0_type2
    logcond <- (cond == 1) * log(1 - exp((datalist[[1]]$not_present_type1 + numimm_type1) * logp0_type1 +
                                         (datalist[[1]]$not_present_type2 + numimm_type2) * logp0_type2))
  }

  loglik <- loglik - logcond
  vec_loglikelihood <- c()

  for (i in 2:length(datalist)) {
    stac <- datalist[[i]]$stac
    brts <- datalist[[i]]$branching_times
    missnumspec <- datalist[[i]]$missing_species
    loglikelihood <- DAISIE_DE_loglik(pars1 = pars1,
                                      brts = brts,
                                      missnumspec = missnumspec,
                                      stac = stac,
                                      methode = methode,
                                      reltolint = reltolint,
                                      abstolint = abstolint,
                                      sampling = sampling)

    vec_loglikelihood <- c(vec_loglikelihood, loglikelihood)

    print_parameters_and_loglik(
      pars = c(stac, parameter),
      loglik = loglikelihood,
      verbose = pars2[4],
      parnames = c("lambda^c", "mu", "K", "gamma", "lambda^a", "prob_init_pres"),
      type = 'clade_loglik'
    )
  }
  loglik <- sum(vec_loglikelihood) + loglik
  return(loglik)
}
