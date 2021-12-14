# Don't document this function. For internal use only.
DAISIE_loglik_all_choosepar4 <- function(trparsopt,
                                         trparsfix,
                                         idparsopt,
                                         idparsfix,
                                         pars2,
                                         datalist,
                                         methode,
                                         CS_version,
                                         abstolint = 1E-16,
                                         reltolint = 1E-10) {
  trpars1 <- rep(0, 6)
  trpars1[idparsopt] <- trparsopt
  if (length(idparsfix) != 0) {
    trpars1[idparsfix] <- trparsfix
  }
  if (max(trpars1) > 1 | min(trpars1) < 0) {
    loglik <- -Inf
  } else {
    pars1 <- trpars1 / (1 - trpars1)
    CS_version$sd <- pars1[6]
    pars1 <- pars1[-6]
    if (min(pars1) < 0) {
      loglik <- -Inf
    } else {
      loglik <- DAISIE::DAISIE_loglik_all(
        pars1 = pars1,
        pars2 = pars2,
        datalist = datalist,
        methode = methode,
        CS_version = CS_version,
        abstolint = abstolint,
        reltolint = reltolint
      )
    }
    if (is.nan(loglik) || is.na(loglik)) {
      cat("There are parameter values used
          which cause numerical problems.\n")
      loglik <- -Inf
    }
  }
  return(loglik)
}

#' Computes MLE for single type species under a clade specific scenario where
#' one parameter may vary over the clades governed by a specific distribution
#'
#' @inheritParams default_params_doc
#'
#' @return The output is a dataframe containing estimated parameters and
#' maximum loglikelihood.
#' \item{lambda_c}{ gives the maximum likelihood
#' estimate of lambda^c, the rate of cladogenesis}
#' \item{mu}{ gives the maximum
#' likelihood estimate of mu, the extinction rate}
#' \item{K}{ gives the maximum
#' likelihood estimate of K, the carrying-capacity}
#' \item{gamma}{ gives the
#' maximum likelihood estimate of gamma, the immigration rate }
#' \item{lambda_a}{ gives the maximum likelihood estimate of lambda^a, the rate
#' of anagenesis}
#' \item{sd}{gives the maximum likelihood estimate of the standard
#' deviation for the parameter which is allowed to vary}
#' \item{loglik}{ gives the maximum loglikelihood}
#' \item{df}{
#' gives the number
#' of estimated parameters, i.e. degrees of freedom}
#' \item{conv}{ gives a
#' message on convergence of optimization; conv = 0 means convergence}
#' @keywords internal
DAISIE_ML4 <- function(
  datalist,
  initparsopt,
  idparsopt,
  parsfix,
  idparsfix,
  res = 100,
  ddmodel = 0,
  cond = 0,
  tol = c(1E-4, 1E-5, 1E-7),
  maxiter = 1000 * round((1.25) ^ length(idparsopt)),
  methode = "odeint::runge_kutta_fehlberg78",
  optimmethod = "subplex",
  CS_version = create_CS_version(model = 2,
                                 relaxed_par = "cladogenesis"),
  verbose = 0,
  tolint = c(1E-16, 1E-10),
  island_ontogeny = NA,
  jitter = 0,
  num_cycles = 1) {

  out2err <- data.frame(
    lambda_c = NA,
    mu = NA,
    K = NA,
    gamma = NA,
    lambda_a = NA,
    sd = NA,
    loglik = NA,
    df = NA,
    conv = NA
  )
  out2err <- invisible(out2err)
  namepars <- c(
    "lambda_c",
    "mu",
    "K",
    "gamma",
    "lambda_a",
    "sd"
  )

  if (length(namepars[idparsopt]) == 0) {
    optstr <- "nothing"
  } else {
    optstr <- namepars[idparsopt]
  }
  cat("You are optimizing", optstr, "\n")
  if (length(namepars[idparsfix]) == 0) {
    fixstr <- "nothing"
  } else {
    fixstr <- namepars[idparsfix]
  }
  cat("You are fixing", fixstr, "\n")
  idpars <- sort(c(idparsopt, idparsfix))
  missnumspec <- unlist(lapply(datalist, function(list) {list$missing_species})) # nolint
  if (sum(missnumspec) > (res - 1)) {
    cat(
      "The number of missing species is too large relative to the
      resolution of the ODE.\n")
    return(out2err)
  }
  if (length(idpars) != 6) {
    cat("You have too few or too many parameters to be optimized or fixed.\n")
    return(out2err)
  }
  if ((prod(idpars == (1:6)) != 1) || # nolint
      (length(initparsopt) != length(idparsopt)) ||
      (length(parsfix) != length(idparsfix))) {
    cat("The parameters to be optimized and/or fixed are incoherent.\n")
    return(out2err)
  }
  cat("Calculating the likelihood for the initial parameters.", "\n")
  utils::flush.console()
  trparsopt <- initparsopt / (1 + initparsopt)
  trparsopt[which(initparsopt == Inf)] <- 1
  trparsfix <- parsfix / (1 + parsfix)
  trparsfix[which(parsfix == Inf)] <- 1
  pars2 <- c(
    res,
    ddmodel,
    cond,
    verbose,
    island_ontogeny
  )

  optimpars <- c(tol, maxiter)
  initloglik <- DAISIE_loglik_all_choosepar4(
    trparsopt = trparsopt,
    trparsfix = trparsfix,
    idparsopt = idparsopt,
    idparsfix = idparsfix,
    pars2 = pars2,
    datalist = datalist,
    methode = methode,
    CS_version = CS_version,
    abstolint = tolint[1],
    reltolint = tolint[2]
  )
  cat(
    "The loglikelihood for the initial parameter values is",
    initloglik,
    "\n"
  )
  if (initloglik == -Inf) {
    cat(
      "The initial parameter values have a likelihood that is equal to 0 or
      below machine precision. Try again with different initial values.\n"
    )
    return(out2err)
  }
  cat("Optimizing the likelihood - this may take a while.", "\n")
  utils::flush.console()
  out <- DDD::optimizer(
    optimmethod = optimmethod,
    optimpars = optimpars,
    fun = DAISIE_loglik_all_choosepar4,
    trparsopt = trparsopt,
    idparsopt = idparsopt,
    trparsfix = trparsfix,
    idparsfix = idparsfix,
    pars2 = pars2,
    datalist = datalist,
    methode = methode,
    CS_version = CS_version,
    abstolint = tolint[1],
    reltolint = tolint[2],
    jitter = jitter,
    num_cycles = num_cycles
  )
  if (out$conv != 0) {
    cat(
      "Optimization has not converged.
      Try again with different initial values.\n")
    out2 <- out2err
    out2$conv <- out$conv
    return(out2)
  }
  MLtrpars <- as.numeric(unlist(out$par))
  MLpars <- MLtrpars / (1 - MLtrpars)
  ML <- as.numeric(unlist(out$fvalues))

  MLpars1 <- rep(0, 6)
  MLpars1[idparsopt] <- MLpars
  if (length(idparsfix) != 0) {
    MLpars1[idparsfix] <- parsfix
  }
  if (MLpars1[3] > 10 ^ 7) {
    MLpars1[3] <- Inf
  }
  out2 <- data.frame(
    lambda_c = MLpars1[1],
    mu = MLpars1[2],
    K = MLpars1[3],
    gamma = MLpars1[4],
    lambda_a = MLpars1[5],
    sd = MLpars1[6],
    loglik = ML,
    df = length(initparsopt),
    conv = unlist(out$conv)
  )
  s1 <- sprintf(
    "Maximum likelihood parameter estimates: lambda_c: %f, mu: %f, K: %f,
      gamma: %f, lambda_a: %f, sd: %f",
    MLpars1[1],
    MLpars1[2],
    MLpars1[3],
    MLpars1[4],
    MLpars1[5],
    MLpars1[6]
  )
  s2 <- sprintf("Maximum loglikelihood: %f", ML)
  cat("\n", s1, "\n", s2, "\n")
  return(invisible(out2))
}
