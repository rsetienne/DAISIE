DAISIE_loglik_IW_M1 <- function(
  pars1,
  pars2,
  datalist,
  brts,
  stac,
  missnumspec,
  methode,
  abstolint = 1E-16,
  reltolint = 1E-14,
  verbose
  ) {
  if(missnumspec > 0) {
    stop('This likelihood computation cannot deal with missing species.')
  }
  if(!(stac %in% c(0,2,4)) & is.null(datalist$all_colonisations)) {
    stop('This likelihood computation must have explicit colonization times or none at all.')
  }
  datalist2 <- list()
  datalist2[[1]] <- list(island_age = max(abs(brts)), not_present = as.numeric(is.null(datalist)))
  if(!is.null(datalist)) {
    pars2[3] <- 0
    datalist2[[2]] <- datalist
    loglik <- DAISIE_loglik_IW(
      pars1 = pars1,
      pars2 = pars2,
      datalist = datalist2,
      methode = methode,
      abstolint = abstolint,
      reltolint = reltolint,
      verbose = verbose)
  } else {
    loglik <- DAISIE_loglik(
      pars1 = pars1,
      pars2 = pars2,
      brts = max(abs(brts)),
      stac = 0,
      missnumspec = missnumspec,
      methode = 'ode45',
      abstolint = abstolint,
      reltolint = reltolint,
      verbose = verbose
    )
  }
  return(loglik)
}
