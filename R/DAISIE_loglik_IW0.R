DAISIE_loglik_IW_M1 <- function(
  pars1,
  pars2,
  datalist,
  brts,
  stac,
  missnumspec,
  methode = "ode45",
  abstolint = 1E-16,
  reltolint = 1E-14,
  verbose
  ) {
  if(missnumspec > 0) {
    stop('This likelihood computation cannot deal with missing species.')
  }
  if(!(stac %in% c(0,2,4))) {
    stop('This likelihood computation must have explicit colonization times.')
  }
  datalist2 <- list()
  datalist2[[1]] <- list(island_age = max(abs(brts)), not_present = as.numeric(is.null(datalist)))
  if(!is.null(datalist)) datalist2[[2]] <- datalist
  loglik <- DAISIE_loglik_IW(
    pars1 = pars1,
    pars2 = pars2,
    datalist = datalist2,
    methode = 'odeint::runge_kutta_fehlberg78',
    abstolint = abstolint,
    reltolint = reltolint,
    verbose = verbose
  )

}
