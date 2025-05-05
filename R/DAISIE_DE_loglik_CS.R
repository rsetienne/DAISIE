
DAISIE_DE_loglik_CS <- function( pars1,
                                 pars2,
                                 datalist,
                                 methode = "lsodes",
                                 abstolint = 1e-15,
                                 reltolint = 1e-15,
                                 equal_extinction = TRUE)

{
  # Apply equal extinction condition AFTER initializing pars1
  if (equal_extinction) {
    pars1[3] <- pars1[2]
  }

  cond <- pars2[3]

  if (length(pars1) == 5) {
    logp0 <- DAISIE_DE_logp0(datalist, pars1, methode)
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
    if (!is.na(pars1[11])) {
      if (pars1[11] < numimm_type2 / numimm | pars1[11] > (1 - numimm_type1 / numimm)) {
        return(-Inf)
      }
      datalist[[1]]$not_present_type2 <- max(0, round(pars1[11] * numimm) - numimm_type2)
      datalist[[1]]$not_present_type1 <- numimm - (length(datalist) - 1) - datalist[[1]]$not_present_type2
    }
    logp0_type1 <- DAISIE_DE_logp0(datalist, pars1[1:5], methode)
    logp0_type2 <- DAISIE_DE_logp0(datalist, pars1[6:10], methode)
    loglik <- datalist[[1]]$not_present_type1 * logp0_type1 + datalist[[1]]$not_present_type2 * logp0_type2
    logcond <- (cond == 1) * log(1 - exp((datalist[[1]]$not_present_type1 + numimm_type1) * logp0_type1 +
                                           (datalist[[1]]$not_present_type2 + numimm_type2) * logp0_type2))
  }

  loglik <- loglik - logcond
  vec_loglikelihood <- c()

  for (i in 2:length(datalist)) {

    stac <- datalist[[i]]$stac
    brts = datalist[[i]]$branching_times
    missnumspec = datalist[[i]]$missing_species

    if (stac == 1) {
      loglikelihood <- DAISIE_DE_logpNE_max_age_coltime(datalist,i,pars1,methode,reltolint,abstolint)
    } else if (stac == 2) {
      if (length(brts) == 2)
        loglikelihood <- DAISIE_DE_logpES(datalist,i,pars1,methode,reltolint,abstolint)
      else
        loglikelihood <- DAISIE_DE_logpEC(datalist,i,pars1,methode,reltolint,abstolint)
    } else if (stac == 3) {
      if (length(brts) == 2)
        loglikelihood <- DAISIE_DE_logpES_mainland(datalist,i,pars1,methode,reltolint,abstolint)
      else
        loglikelihood <- DAISIE_DE_logpEC_mainland(datalist,i,pars1,methode,reltolint,abstolint)
    } else if (stac == 4) {
      loglikelihood <- DAISIE_DE_logpNE(datalist,i,pars1,methode,reltolint,abstolint)
    } else if (stac == 5) {
      loglikelihood <- DAISIE_DE_logpES_max_age_coltime(datalist,i,pars1,methode,reltolint,abstolint)
    } else if (stac == 6) {
      loglikelihood <- DAISIE_DE_logpEC_max_age_coltime(datalist,i,pars1,methode,reltolint,abstolint)
    } else if (stac == 7) {
      if (length(brts) == 2)
        loglikelihood <- DAISIE_DE_logpES_max_age_coltime_and_mainland(datalist,i,pars1,methode,reltolint,abstolint)
      else
        loglikelihood <- DAISIE_DE_logpEC_max_age_coltime_and_mainland(datalist,i,pars1,methode,reltolint,abstolint)
    } else if (stac == 8) {
      loglikelihood <- DAISIE_DE_logpNE_max_min_age_coltime(datalist,i,pars1,methode,reltolint,abstolint)
    } else if (stac == 9) {
      loglikelihood <- DAISIE_DE_logpES_max_min_age_coltime(datalist,i,pars1,methode,reltolint,abstolint)
    } else {
      stop("Unknown stac value: ", stac)
    }

    vec_loglikelihood <- c(vec_loglikelihood, loglikelihood)

    print_parameters_and_loglik(
      pars = c(stac, pars1),
      loglik = loglikelihood,
      verbose = pars2[4],
      parnames = c("lambda^c", "mu1", "mu2", "gamma", "lambda^a", "prob_init_pres"),
      type = 'clade_loglik'
    )
  }
  loglik <- sum(vec_loglikelihood) + loglik
  return(loglik)
}
