DAISIE_loglik_IW_choosepar <- function(
  trparsopt,
  trparsfix,
  idparsopt,
  idparsfix,
  pars2,
  M,
  datalist,
  methode = "ode45",
  abstolint = 1E-16,
  reltolint = 1E-14
) {
  trpars1 <- rep(0, 5)
  trpars1[idparsopt] <- trparsopt
  if (length(idparsfix) != 0) {
    trpars1[idparsfix] <- trparsfix
  }
  if (max(trpars1) > 1 | min(trpars1) < 0) {
    loglik <- -Inf
  } else {
    pars1 <- c(trpars1 / (1 - trpars1), M)
    if (min(pars1) < 0) {
      loglik <- -Inf
    } else {
      loglik = DAISIE_loglik_IW(
        pars1 = pars1,
        pars2 = pars2,
        datalist = datalist,
        methode = methode,
        abstolint = abstolint,
        reltolint = reltolint
      )
    }
    if (is.nan(loglik) || is.na(loglik)) {
      warning("There are parameter values used which cause numerical problems.")
      loglik <- -Inf
    }
  }
  return(loglik)
}



#' Maximization of the loglikelihood under the DAISIE model with island-wide
#' diversity-dependence
#'
#' This function computes the maximum likelihood estimates of the parameters of
#' the DAISIE model with island-wide diversity-dependence for data from
#' lineages colonizing an island. It also outputs the corresponding
#' loglikelihood that can be used in model comparisons.
#'
#' The result of sort(c(idparsopt, idparsfix)) should be identical to c(1:5).
#' If not, an error is reported that the input is incoherent. The same happens
#' when the length of initparsopt is different from the length of idparsopt,
#' and the length of parsfix is different from the length of idparsfix.\cr
#'
#' @inheritParams default_params_doc
#'
#' @return The output is a dataframe containing estimated parameters and
#' maximum loglikelihood.  \item{lambda_c}{ gives the maximum likelihood
#' estimate of lambda^c, the rate of cladogenesis} \item{mu}{ gives the maximum
#' likelihood estimate of mu, the extinction rate} \item{K}{ gives the maximum
#' likelihood estimate of K, the carrying-capacity} \item{gamma}{ gives the
#' maximum likelihood estimate of gamma, the immigration rate }
#' \item{lambda_a}{ gives the maximum likelihood estimate of lambda^a, the rate
#' of anagenesis} \item{loglik}{ gives the maximum loglikelihood} \item{df}{
#' gives the number of estimated parameters, i.e. degrees of feedom}
#' \item{conv}{ gives a message on convergence of optimization; conv = 0 means
#' convergence}
#' @author Rampal S. Etienne
#' @seealso \code{\link{DAISIE_loglik_IW}}, \code{\link{DAISIE_ML_CS}}
#' \code{\link{DAISIE_sim_cr}}
#' @references Valente, L.M., A.B. Phillimore and R.S. Etienne (2015).
#' Equilibrium and non-equilibrium dynamics simultaneously operate in the
#' Galapagos islands. Ecology Letters 18: 844-852. <DOI:10.1111/ele.12461>.
#' @keywords models
#' @export DAISIE_ML_IW
DAISIE_ML_IW <- function(
  datalist,
  initparsopt,
  idparsopt,
  parsfix,
  idparsfix,
  res = 100,
  ddmodel = 11,
  cond = 0,
  tol = c(1E-4, 1E-5, 1E-7),
  maxiter = 1000 * round((1.25) ^ length(idparsopt)),
  methode = "ode45",
  optimmethod = "subplex",
  verbose = 0,
  tolint = c(1E-16, 1E-14),
  jitter = 0,
  num_cycles = 1) {

  out2err <- data.frame(lambda_c = NA, mu = NA, K = NA, gamma = NA, lambda_a = NA, loglik = NA, df = NA, conv = NA)
  out2err <- invisible(out2err)
  if (is.null(datalist[[1]]$brts_table)) {
    datalist <- add_brt_table(datalist)
  }
  np <- datalist[[1]]$not_present
  if (is.null(np)) {
    np <- datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2
  }
  np = datalist[[1]]$not_present
  if (is.null(np)) {
    np = datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2
    warning('Number of species not present is misspecified.\n')
    return(invisible(out2err))
  }
  M <- length(datalist) - 1 + np

  idpars <- sort(c(idparsopt, idparsfix))
  if ((prod(idpars == (1:5)) != 1) || (length(initparsopt) != length(idparsopt)) || (length(parsfix) != length(idparsfix))) {
    warning("The parameters to be optimized and/or fixed are incoherent.\n")
    return(out2err)
  }
  if (length(idparsopt) > 11) {
    warning("The number of parameters to be optimized is too high.\n")
    return(out2err)
  }
  namepars <- c("lambda_c", "mu", "K'", "gamma", "lambda_a")

  print_ml_par_settings(
    namepars = namepars,
    idparsopt = idparsopt,
    idparsfix = idparsfix,
    idparsnoshift = NA,
    all_no_shift = NA,
    verbose = verbose
  )

  trparsopt <- initparsopt / (1 + initparsopt)
  trparsopt[which(initparsopt == Inf)] <- 1
  trparsfix <- parsfix / (1 + parsfix)
  trparsfix[which(parsfix == Inf)] <- 1
  pars2 <- c(res, ddmodel, cond, verbose)
  optimpars <- c(tol, maxiter)
  initloglik <- DAISIE_loglik_IW_choosepar(trparsopt = trparsopt, trparsfix = trparsfix, idparsopt = idparsopt, idparsfix = idparsfix, M = M, pars2 = pars2, datalist = datalist, methode = methode, abstolint = tolint[1], reltolint = tolint[2])
  message("The loglikelihood for the initial parameter values is ", initloglik)
  if (initloglik == -Inf) {
    warning("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n")
    return(out2err)
  }
  message("Optimizing the likelihood - this may take a while.")
  out <- DDD::optimizer(
    optimmethod = optimmethod,
    optimpars = optimpars,
    fun = DAISIE_loglik_IW_choosepar,
    trparsopt = trparsopt,
    idparsopt = idparsopt,
    trparsfix = trparsfix,
    idparsfix = idparsfix,
    M = M,
    pars2 = pars2,
    datalist = datalist,
    methode = methode,
    abstolint = tolint[1],
    reltolint = tolint[2],
    jitter = jitter,
    num_cycles = num_cycles
  )
  if (out$conv != 0) {
    warning("Optimization has not converged. Try again with different initial values.")
    out2 <- out2err
    out2$conv <- out$conv
    return(out2)
  }
  MLtrpars <- as.numeric(unlist(out$par))
  MLpars <- MLtrpars / (1 - MLtrpars)
  ML <- as.numeric(unlist(out$fvalues))
  MLpars1 <- rep(0, 5)
  MLpars1[idparsopt] <- MLpars
  if (length(idparsfix) != 0) { MLpars1[idparsfix] <- parsfix }
  if (MLpars1[3] > 10 ^ 7){ MLpars1[3] <- Inf }
  out2 <- data.frame(lambda_c = MLpars1[1], mu = MLpars1[2], K = MLpars1[3], gamma = MLpars1[4], lambda_a = MLpars1[5], loglik = ML, df = length(initparsopt), conv = unlist(out$conv))
  print_parameters_and_loglik(pars = MLpars1[1:5],
                              loglik = ML,
                              verbose = TRUE,
                              type = 'island_ML')
  return(invisible(out2))
}
