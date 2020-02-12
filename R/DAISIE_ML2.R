DAISIE_loglik_all_choosepar2 <- function(
  trparsopt,
  trparsfix,
  idparsopt,
  idparsfix,
  idparsmat,
  pars2,
  datalist,
  methode,
  abstolint = 1E-16,
  reltolint = 1E-10
  ) {
   trpars1 <- 0 * idparsmat
   trpars1[idparsopt] <- trparsopt
   if (length(idparsfix) != 0) {
     for (i in 1:length(idparsfix)) {
       trpars1[which(idparsmat == idparsfix[i])] <- trparsfix[i]
     }
   }
   for (i in 1:length(idparsopt)) {
     trpars1[which(idparsmat == idparsopt[i])] <- trparsopt[i]
   }
   if (max(trpars1) > 1 || min(trpars1) < 0) {
      loglik <- -Inf
   } else {
      pars1 <- trpars1 / (1 - trpars1)
      loglik <- 0
      for (i in 1:length(datalist)) {
        loglik <- loglik + DAISIE_loglik_all(pars1 = pars1[i, ], pars2 = pars2, datalist = datalist[[i]], methode = methode, abstolint = abstolint, reltolint = reltolint)
      }
   }
   if (is.nan(loglik) || is.na(loglik)) {
      cat("There are parameter values used which cause numerical problems.\n")
      loglik <- -Inf
   }
   return(loglik)
}

#' Computes MLE for two sets of species under a clade specific scenario
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
#' of anagenesis} \item{lambda_c2}{ gives the maximum likelihood estimate of
#' lambda^c2, the rate of cladogenesis for the second group of
#' species} \item{mu2}{ gives the maximum likelihood estimate of mu2, the
#' extinction rate for the second group of species} \item{K2}{ gives
#' the maximum likelihood estimate of K2, the carrying-capacity for the
#'  second group of species} \item{gamma2}{ gives the maximum
#' likelihood estimate of gamma2, the immigration rate for the second
#' group of species} \item{lambda_a2}{ gives the maximum likelihood estimate of
#' lambda^a2, the rate of anagenesis for the second group of species}
#' \item{loglik}{ gives the maximum loglikelihood} \item{df}{ gives the number
#' of estimated parameters, i.e. degrees of feedom} \item{conv}{ gives a
#' message on convergence of optimization; conv = 0 means convergence}
#'
DAISIE_ML2 <- function(
  datalist,
  initparsopt,
  idparsopt,
  parsfix,
  idparsfix,
  idparsmat,
  res = 100,
  ddmodel = 0,
  cond = 0,
  island_ontogeny = NA,
  tol = c(1E-4, 1E-5, 1E-7),
  maxiter = 1000 * round((1.25) ^ length(idparsopt)),
  methode = "lsodes",
  optimmethod = "subplex",
  verbose = 0,
  tolint = c(1E-16, 1E-10)
  ) {
# datalist = list of all data: branching times, status of clade, and numnber of missing species
# datalist[[,]][1] = list of branching times (positive, from present to past)
# - max(brts) = age of the island
# - next largest brts = stem age / time of divergence from the mainland
# The interpretation of this depends on stac (see below)
# For stac = 0, this needs to be specified only once.
# For stac = 1, this is the time since divergence from the immigrant's sister on the mainland.
# The immigrant must have immigrated at some point since then.
# For stac = 2 and stac = 3, this is the time since divergence from the mainland.
# The immigrant that established the clade on the island must have immigrated precisely at this point.
# For stac = 3, it must have reimmigrated, but only after the first immigrant had undergone speciation.
# - min(brts) = most recent branching time (only for stac = 2, or stac = 3)
# datalist[[,]][2] = list of status of the clades formed by the immigrant
#  . stac == 0 : immigrant is not present and has not formed an extant clade
# Instead of a list of zeros, here a number must be given with the number of clades having stac = 0
#  . stac == 1 : immigrant is present but has not formed an extant clade
#  . stac == 2 : immigrant is not present but has formed an extant clade
#  . stac == 3 : immigrant is present and has formed an extant clade
#  . stac == 4 : immigrant is present but has not formed an extant clade, and it is known when it immigrated.
#  . stac == 5 : immigrant is not present and has not formed an extent clade, but only an endemic species
# datalist[[,]][3] = list with number of missing species in clades for stac = 2 and stac = 3;
# for stac = 0 and stac = 1, this number equals 0.
# idparsmat = matrix containing the ids of the parameters, linking them to initparsopt and parsfix.
# Per island system we use the following order
# - lac = (initial) cladogenesis rate
# - mu = extinction rate
# - K = maximum number of species possible in the clade
# - gam = (initial) immigration rate
# - laa = (initial) anagenesis rate
# Example: idparsmat = rbind(c(1,2,3,4,5),c(1,2,3,6,7)) has different rates of colonization and anagenesis for the two islands.
# initparsopt, parsfix = values of optimized and fixed model parameters
# idparsopt, idparsfix = ids of optimized and fixed model parameters
# - res = pars2[1] = lx = length of ODE variable x
# - ddmodel = pars2[2] = diversity-dependent model,mode of diversity-dependence
#  . ddmodel == 0 : no diversity-dependence
#  . ddmodel == 1 : linear dependence in speciation rate (anagenesis and cladogenesis)
#  . ddmodel == 11 : linear dependence in speciation rate and immigration rate
#  . ddmodel == 3 : linear dependence in extinction rate
# - cond = conditioning; currently only cond = 0 is possible
#  . cond == 0 : no conditioning
#  . cond == 1 : conditioning on presence on the island

  options(warn = -1)
  out2err <- data.frame(lambda_c = NA, mu = NA, K = NA, gamma = NA, lambda_a = NA, loglik = NA, df = NA, conv = NA)
  out2err <- invisible(out2err)
  numisl <- length(datalist)
  missnumspec <- 0
  for (i in 1:numisl) {
    missnumspec <- missnumspec + sum(unlist(lapply(datalist[[i]], function(list) {list$missing_species})))
  }

  if (missnumspec > (res - 1)) {
    cat("The number of missing species is too large relative to the resolution of the ODE.\n")
    return(out2err)
  }
  if (all((sort(unique(as.vector(idparsmat))) != sort(c(idparsopt, idparsfix)))) ||
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
  pars2 <- c(res, ddmodel, cond, 0, island_ontogeny)
  optimpars <- c(tol, maxiter)
  initloglik <- DAISIE_loglik_all_choosepar2(trparsopt = trparsopt, trparsfix = trparsfix, idparsopt = idparsopt, idparsfix = idparsfix, idparsmat = idparsmat, pars2 = pars2, datalist = datalist, methode, abstolint = tolint[1], reltolint = tolint[2])
  cat("The loglikelihood for the initial parameter values is", initloglik, "\n")
  if (initloglik == -Inf) {
    cat("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n")
    return(out2err)
  }
  cat("Optimizing the likelihood - this may take a while.", "\n")
  utils::flush.console()
  out <- DDD::optimizer(optimmethod = optimmethod, optimpars = optimpars, fun = DAISIE_loglik_all_choosepar2, trparsopt = trparsopt, idparsopt = idparsopt, trparsfix = trparsfix, idparsfix = idparsfix, idparsmat = idparsmat, pars2 = pars2, datalist = datalist, methode = methode, abstolint = tolint[1], reltolint = tolint[2])
  if (out$conv != 0) {
    cat("Optimization has not converged. Try again with different initial values.\n")
    out2 <- out2err
    out2$conv <- out$conv
    return(out2err)
  }
  MLtrpars <- as.numeric(unlist(out$par))
  MLpars <- MLtrpars / (1 - MLtrpars)
  ML <- as.numeric(unlist(out$fvalues))
  MLpars1 <- 0 * idparsmat
  if (length(idparsfix) != 0) {
    for (i in 1:length(idparsfix)) {
      MLpars1[which(idparsmat == idparsfix[i])] <- parsfix[i]
    }
  }
  for (i in 1:length(idparsopt)) {
    MLpars1[which(idparsmat == idparsopt[i])] <- MLpars[i]
  }
  for (i in 1:numisl) {
    if (MLpars1[i, 3] > 10 ^ 7) {
      MLpars1[i, 3] <- Inf
    }
  }
  out2 <- data.frame(lambda_c = MLpars1[, 1], mu = MLpars1[, 2], K = MLpars1[, 3], gamma = MLpars1[, 4], lambda_a = MLpars1[, 5], loglik = ML, df = length(initparsopt), conv = unlist(out$conv))
  s1 <- sprintf("Maximum likelihood parameter estimates: %f", MLpars1)
  s2 <- sprintf("Maximum loglikelihood: %f", ML)
  cat("\n", s1, "\n", s2, "\n")
  return(invisible(out2))
}
