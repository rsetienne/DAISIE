DAISIE_DE_loglik_CS <- function(
    pars1 = pars1,
    pars2 = pars2,
    datalist = datalist,
    methode = methode,
    rtol,
    atol) {

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
    logcond <- (cond == 1) * log(1 - exp((datalist[[1]]$not_present_type1 + numimm_type1) * logp0_type1 + (datalist[[1]]$not_present_type2 + numimm_type2) * logp0_type2))
  }
  loglik <- loglik - logcond
  vec_likelihood <- c()
  for (i in 2:length(datalist)) {
    likelihood <- switch(datalist[[i]]$stac,
                         `1` = DAISIE_DE_logpNE_max_age_coltime(datalist, i, pars1, methode, rtol, atol),
                         `2` = if (length(datalist[[i]]$branching_times) == 2 && datalist[[i]]$missing_species == 0)
                           DAISIE_DE_logpES(datalist, i, pars1, methode, rtol, atol)
                         else
                           DAISIE_DE_logpEC(datalist, i, pars1, methode, rtol, atol),
                         `3` = if (length(datalist[[i]]$branching_times) == 2 && datalist[[i]]$missing_species == 0)
                           DAISIE_DE_logpES_mainland(datalist, i, pars1, methode, rtol, atol)
                         else
                           DAISIE_DE_logpEC_mainland(datalist, i, pars1, methode, rtol, atol),
                         `4` = DAISIE_DE_logpNE(datalist, i, pars1, methode, rtol, atol),
                         `5` = DAISIE_DE_logpES_max_age_coltime(datalist, i, pars1, methode, rtol, atol),
                         `6` = if (length(datalist[[i]]$branching_times) > 2)
                           DAISIE_DE_logpEC_max_age_coltime(datalist, i, pars1, methode, rtol, atol)
                         else
                           function_Factor_Loglik_EC_max_Age_approximation_2(datalist, i, pars1)
    )
    vec_likelihood <- c(vec_likelihood, likelihood)
    cat(sprintf('Status of colonist: %d, Parameters: %f %f %f %f %f, Loglikelihood: %f\n',
                datalist[[i]]$stac, pars1[1], pars1[2], pars1[3], pars1[4], pars1[5], likelihood))
  }
  cat(sprintf('Parameters: %f %f %f %f %f, Loglikelihood: %f\n', pars1[1], pars1[2], pars1[3], pars1[4], pars1[5], sum(vec_likelihood) + loglik))
  return(sum(vec_likelihood) + loglik)
}
