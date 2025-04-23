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

    # Print in custom format
    cat("Status of colonist:", stac,
        ", Parameters:", paste(format(pars1, digits = 7), collapse = ", "),
        ", Loglikelihood:", format(loglikelihood, digits = 7), "\n")

    loglik_vector <- loglikelihood


  return(loglik_vector)
}


