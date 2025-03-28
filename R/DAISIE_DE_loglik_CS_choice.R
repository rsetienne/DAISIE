DAISIE_DE_loglik_CS_choice <- function(
                                       datalist,
                                       brts,
                                       pars1,
                                       missnumspec,
                                       stac,
                                       methode = "lsodes",
                                       abstolint = 1e-16,
                                       reltolint = 1e-16,
                                       equal_extinction = TRUE) {

  equal_extinction <- equal_extinction

  # Apply equal extinction condition AFTER initializing pars1
  if (equal_extinction) {
    pars1[3] <- pars1[2]
  }
  loglik_vector <- numeric()


    if (stac == 1) {
      loglikelihood <- DAISIE_DE_logpNE_max_age_coltime(datalist, brts, pars1, missnumspec, methode, reltolint, abstolint)
    } else if (stac == 2) {
      if (length(brts) == 2)
        loglikelihood <- DAISIE_DE_logpES(datalist, brts, pars1, missnumspec, methode, reltolint, abstolint)
      else
        loglikelihood <- DAISIE_DE_logpEC(datalist, brts, pars1, missnumspec, methode, reltolint, abstolint)
    } else if (stac == 3) {
      if (length(brts) == 2)
        loglikelihood <- DAISIE_DE_logpES_mainland(datalist, brts, pars1, missnumspec, methode, reltolint, abstolint)
      else
        loglikelihood <- DAISIE_DE_logpEC_mainland(datalist, brts, pars1, missnumspec, methode, reltolint, abstolint)
    } else if (stac == 4) {
      loglikelihood <- DAISIE_DE_logpNE(datalist, brts, pars1, missnumspec, methode, reltolint, abstolint)
    } else if (stac == 5) {
      loglikelihood <- DAISIE_DE_logpES_max_age_coltime(datalist, brts, pars1, missnumspec, methode, reltolint, abstolint)
    } else if (stac == 6) {
      loglikelihood <- DAISIE_DE_logpEC_max_age_coltime(datalist, brts, pars1, missnumspec, methode, reltolint, abstolint)
    } else if (stac == 7) {
      if (length(brts) == 2)
        loglikelihood <- DAISIE_DE_logpES_max_age_coltime_and_mainland(datalist, brts, pars1, missnumspec, methode, reltolint, abstolint)
      else
        loglikelihood <- DAISIE_DE_logpEC_max_age_coltime_and_mainland(datalist, brts, pars1, missnumspec, methode, reltolint, abstolint)
    } else if (stac == 8) {
      loglikelihood <- DAISIE_DE_logpNE_max_min_age_coltime(datalist, brts, pars1, missnumspec, methode, reltolint, abstolint)
    } else if (stac == 9) {
      loglikelihood <- DAISIE_DE_logpES_max_min_age_coltime(datalist, brts, pars1, missnumspec, methode, reltolint, abstolint)
    } else {
      stop("Unknown stac value: ", stac)
    }

    # Print in custom format
    cat("Status of colonist:", stac,
        ", Parameters:", paste(format(pars1, digits = 7), collapse = ", "),
        ", Loglikelihood:", format(loglikelihood, digits = 7), "\n")

    loglik_vector <- loglikelihood


  return(loglik_vector)
}

i <- 7
DAISIE_DE_loglik_CS_choice(datalist = datalist,
                                     pars1 = pars1,
                                     brts = datalist[[i]]$branching_times,
                                     stac = datalist[[i]]$stac,
                                     missnumspec = datalist[[i]]$missing_species,
                                     methode = "lsodes",
                                     abstolint = 1e-15,
                                     reltolint = 1e-15,
                                     equal_extinction = TRUE)

#pars2 <- c(300,  11,   0,   2)
#pars11 <- pars1
#pars11[3] <- Inf
DAISIE:::DAISIE_loglik_CS_choice(pars1 = pars11,
                                 pars2 = pars2,
                                 datalist = datalist,
                                 brts = datalist[[i]]$branching_times,
                                 stac = datalist[[i]]$stac,
                                 missnumspec = datalist[[i]]$missing_species,
                                 methode = "lsodes",
                                 abstolint =  1e-16,
                                 reltolint =  1e-16)

