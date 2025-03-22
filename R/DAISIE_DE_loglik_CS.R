DAISIE_DE_loglik_CS <- function(
    pars1 = pars1,
    pars2 = pars2,
    datalist = datalist,
    methode = methode,
    abstolint = 1E-16,
    reltolint = 1E-16,
    equal_extinction = TRUE)

  {
  pars1 <- pars1

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
    logcond <- (cond == 1) * log(1 - exp((datalist[[1]]$not_present_type1 + numimm_type1) * logp0_type1 + (datalist[[1]]$not_present_type2 + numimm_type2) * logp0_type2))
  }
  loglik <- loglik - logcond
  vec_loglikelihood <- c()
  for (i in 2:length(datalist)) {
    loglikelihood <- switch(datalist[[i]]$stac,
                            `1` = DAISIE_DE_logpNE_max_age_coltime(datalist, i, pars1, methode, reltolint, abstolint),
                            `2` = if (length(datalist[[i]]$branching_times) == 2)
                              DAISIE_DE_logpES(datalist, i, pars1, methode, reltolint, abstolint)
                            else
                              DAISIE_DE_logpEC(datalist, i, pars1, methode, reltolint, abstolint),
                            `3` = if (length(datalist[[i]]$branching_times) == 2)
                              DAISIE_DE_logpES_mainland(datalist, i, pars1, methode, reltolint, abstolint)
                            else
                              DAISIE_DE_logpEC_mainland(datalist, i, pars1, methode, reltolint, abstolint),
                            `4` = DAISIE_DE_logpNE(datalist, i, pars1, methode, reltolint, abstolint),
                            `5` = DAISIE_DE_logpES_max_age_coltime(datalist, i, pars1, methode, reltolint, abstolint),
                            `6` = DAISIE_DE_logpEC_max_age_coltime(datalist, i, pars1, methode, reltolint, abstolint)
    )
    vec_loglikelihood <- c(vec_loglikelihood, loglikelihood)

    DAISIE:::print_parameters_and_loglik(pars = c(datalist[[i]]$stac,pars1),
                                loglik = loglikelihood,
                                verbose = pars2[4],
                                parnames = c("lambda^c", "mu1", "mu2", "gamma", "lambda^a", "prob_init_pres"),
                                type = 'clade_loglik')
  }
  loglik <- sum(vec_loglikelihood) + loglik
  return(loglik)
}










DAISIE_loglik_CS <- DAISIE_loglik_all <- function(
    pars1,
    pars2,
    datalist,
    methode = "lsodes",
    CS_version = 1,
    abstolint = 1E-16,
    reltolint = 1E-10) {
  if (length(pars1) == 14) {
    if (datalist[[1]]$island_age > pars1[11]) {
      stop(
        "The island age in the area parameters is inconsistent with the island
        data."
      )
    }
    peak <- calc_peak(
      total_time = datalist[[1]]$island_age,
      area_pars = create_area_pars(
        max_area = pars1[8],
        current_area = pars1[9],
        proportional_peak_t = pars1[10],
        total_island_age = pars1[11],
        sea_level_amplitude = pars1[12],
        sea_level_frequency = pars1[13],
        island_gradient_angle = pars1[14]
      )
    )
    pars1 <- c(
      pars1,
      island_ontogeny = pars2[5],
      sea_level = pars2[6],
      datalist[[1]]$island_age,
      peak
    )
  }

  pars1 <- as.numeric(pars1)
  cond <- pars2[3]
  if (length(pars1) == 6) {
    endpars1 <- 6
  } else {
    endpars1 <- 5
  }

  if(length(pars1) %in% c(5,6) | !is.na(pars2[5])) {
    if(!is.na(pars2[5]))
    {
      endpars1 <- length(pars1)
    }
    logp0 <- DAISIE_loglik_CS_choice(
      pars1 = pars1,
      pars2 = pars2,
      brts = datalist[[1]]$island_age,
      stac = 0,
      missnumspec = 0,
      methode = methode,
      CS_version = CS_version,
      abstolint = abstolint,
      reltolint = reltolint
    )
    if(logp0 >= 0 & pars1[2]/pars1[1] > 100)
    {
      logp0 <- approximate_logp0(gamma = pars1[4],
                                 mu = pars1[2],
                                 t = datalist[[1]]$island_age)
    }
    if(logp0 >= 0)
    {
      message('Positive values of loglik encountered without possibility for approximation. Setting loglik to -Inf.')
      loglik <- -Inf
      print_parameters_and_loglik(pars = pars,
                                  loglik = loglik,
                                  verbose = pars2[4],
                                  type = 'island_loglik')
      return(loglik)
    }
    if (is.null(datalist[[1]]$not_present)) {
      loglik <- (datalist[[1]]$not_present_type1 +
                   datalist[[1]]$not_present_type2) * logp0
      numimm <- (datalist[[1]]$not_present_type1 +
                   datalist[[1]]$not_present_type2) + length(datalist) - 1
    } else {
      loglik <- datalist[[1]]$not_present * logp0
      numimm <- datalist[[1]]$not_present + length(datalist) - 1
    }
    logcond <- logcondprob(numcolmin = cond,numimm = numimm,logp0 = logp0)
    if (length(datalist) > 1) {
      for (i in 2:length(datalist)) {
        datalist[[i]]$type1or2 <- 1
      }
    }
  } else {
    numimm <- datalist[[1]]$not_present_type1 +
      datalist[[1]]$not_present_type2 + length(datalist) - 1
    numimm_type2 <- length(
      which(unlist(datalist)[which(names(unlist(datalist)) == "type1or2")] == 2)
    )
    numimm_type1 <- length(datalist) - 1 - numimm_type2

    if (is.na(pars1[11]) == FALSE && length(pars1) == 11) {
      if (pars1[11] < numimm_type2 / numimm |
          pars1[11] > (1 - numimm_type1 / numimm)) {
        return(-Inf)
      }
      datalist[[1]]$not_present_type2 <- max(
        0,
        round(pars1[11] * numimm) - numimm_type2
      )
      datalist[[1]]$not_present_type1 <- numimm - (length(datalist) - 1) -
        datalist[[1]]$not_present_type2
    }
    logp0_type1 <- DAISIE_loglik_CS_choice(
      pars1 = pars1[1:5],
      pars2 = pars2,
      brts = datalist[[1]]$island_age,
      stac = 0,
      missnumspec = 0,
      methode = methode,
      CS_version = CS_version,
      abstolint = abstolint,
      reltolint = reltolint
    )
    if(logp0_type1 >= 0 & pars1[2]/pars1[1] > 100)
    {
      logp0_type1 <- approximate_logp0(gamma = pars1[4], mu = pars1[2], t = datalist[[1]]$island_age)
    }
    logp0_type2 <- DAISIE_loglik_CS_choice(
      pars1 = pars1[6:10],
      pars2 = pars2,
      brts = datalist[[1]]$island_age,
      stac = 0,
      missnumspec = 0,
      methode = methode,
      CS_version = CS_version,
      abstolint = abstolint,
      reltolint = reltolint
    )
    if(logp0_type2 >= 0 & pars1[7]/pars1[6] > 100)
    {
      logp0_type2 <- approximate_logp0(gamma = pars1[9], mu = pars1[7], t = datalist[[1]]$island_age)
    }
    loglik <- datalist[[1]]$not_present_type1 * logp0_type1 +
      datalist[[1]]$not_present_type2 * logp0_type2
    #logcond <- (cond == 1) *
    #  log(1 - exp((datalist[[1]]$not_present_type1 + numimm_type1) *
    #                logp0_type1 +
    #                (datalist[[1]]$not_present_type2 + numimm_type2) *
    #                logp0_type2))
    logcond <- logcondprob(numcolmin = cond,
                           numimm = c(datalist[[1]]$not_present_type1 + numimm_type1,datalist[[1]]$not_present_type2 + numimm_type2),
                           logp0 = c(logp0_type1,logp0_type2) )
  }
  loglik <- loglik - logcond

  if(length(datalist) > 1)
  {
    for(i in 2:length(datalist))
    {
      if(datalist[[i]]$type1or2 == 1)
      {
        pars <- pars1[1:endpars1]
      } else {
        pars <- pars1[6:10]
      }
      loglik <- loglik + DAISIE_loglik_CS_choice(
        pars1 = pars,
        pars2 = pars2,
        datalist = datalist[[i]],
        brts = datalist[[i]]$branching_times,
        stac = datalist[[i]]$stac,
        missnumspec = datalist[[i]]$missing_species,
        methode = methode,
        CS_version = CS_version,
        abstolint = abstolint,
        reltolint = reltolint)
    }
  }

  print_parameters_and_loglik(pars = pars,
                              loglik = loglik,
                              verbose = pars2[4],
                              parnames = c("lambda^c", "mu", "K", "gamma", "lambda^a", "prob_init_pres"),
                              type = 'island_loglik')
  return(loglik)
}
