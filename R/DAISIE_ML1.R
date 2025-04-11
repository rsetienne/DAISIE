# Don't document this function. For internal use only.
DAISIE_loglik_all_choosepar <- function(trparsopt,
                                        trparsfix,
                                        idparsopt,
                                        idparsfix,
                                        idparsnoshift,
                                        idparseq,
                                        pars2,
                                        datalist,
                                        methode,
                                        CS_version = 1,
                                        abstolint = 1E-16,
                                        reltolint = 1E-10,
                                        equal_extinction = FALSE) {
  all_no_shift <- 6:10
  non_oceanic_option <- FALSE
  if (max(idparsopt,-Inf) <= 6 &&
      max(idparsfix,-Inf) <= 6 &&
      (6 %in% idparsopt || 6 %in% idparsfix)) {
    idparsnoshift <- 7:11
    all_no_shift <- 7:11
    non_oceanic_option <- TRUE
  }
  if (sum(idparsnoshift %in% (all_no_shift)) != 5) {
    trpars1 <- rep(0, 11)
  } else {
    trpars1 <- rep(0, 6)
    prop_type2_present <- which(idparsfix == 11)
    if (length(prop_type2_present) > 0) {
       trparsfix <- trparsfix[-prop_type2_present]
       idparsfix <- idparsfix[-prop_type2_present]
    }
  }
  trpars1[idparsopt] <- trparsopt
  if (length(idparsfix) != 0) {
    trpars1[idparsfix] <- trparsfix
  }
  if (sum(idparsnoshift %in% all_no_shift) != 5) {
    trpars1[idparsnoshift] <- trpars1[idparsnoshift - 5]
  }
  if (max(trpars1) > 1 | min(trpars1) < 0) {
    loglik <- -Inf
  } else {
    pars1 <- trpars1 / (1 - trpars1)
    if (pars2[6] > 0) {
      pars1 <- DAISIE_eq(datalist, pars1, pars2[-5])
      if (sum(idparsnoshift %in% all_no_shift) != 5) {
        pars1[idparsnoshift] <- pars1[idparsnoshift - 5]
      }
    }
    if (min(pars1) < 0 | (pars1[6] > 1 && non_oceanic_option == TRUE)) {
      loglik <- -Inf
    } else {
      loglik <- DAISIE_loglik_all(
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
      message("There are parameter values used which cause numerical problems.")
      loglik <- -Inf
    }
  }
  return(loglik)
}

#' Computes MLE for single type species under a clade specific scenario
#'
#' @inheritParams default_params_doc
#'
#' @keywords internal
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
#' \item{loglik}{ gives the maximum loglikelihood}
#' \item{df}{
#' gives the number
#' of estimated parameters, i.e. degrees of feedom}
#' \item{conv}{ gives a
#' message on convergence of optimization; conv = 0 means convergence}
DAISIE_ML1 <- function(
  datalist,
  initparsopt,
  idparsopt,
  parsfix,
  idparsfix,
  idparsnoshift = 6:10,
  res = 100,
  ddmodel = 0,
  cond = 0,
  eqmodel = 0,
  x_E = 0.95,
  x_I = 0.98,
  tol = c(1E-4, 1E-5, 1E-7),
  maxiter = 1000 * round((1.25) ^ length(idparsopt)),
  methode = "lsodes",
  optimmethod = "subplex",
  CS_version = 1,
  verbose = 0,
  tolint = c(1E-16, 1E-10),
  island_ontogeny = NA,
  jitter = 0,
  num_cycles = 1,
  equal_extinction = FALSE) {
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
  # initparsopt, parsfix = optimized and fixed model parameters
  # - pars1[1] = lac = (initial) cladogenesis rate
  # - pars1[2] = mu = extinction rate
  # - pars1[3] = K = maximum number of species possible in the clade
  # - pars1[4] = gam = (initial) immigration rate
  # - pars1[5] = laa = (initial) anagenesis rate
  # - pars1[6]...pars1[10] = same as pars1[1]...pars1[5], but for a second type of immigrant
  # - pars1[11] = proportion of type 2 immigrants in species pool
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
  # - eqmodel = equilibrium model
  #  . eqmodel = 0 : no equilibrium is assumed
  #  . eqmodel = 1 : equilibrium is assumed on deterministic equation for total number of species
  #  . eqmodel = 2 : equilibrium is assumed on total number of species using deterministic equation for endemics and immigrants
  #  . eqmodel = 3 : equilibrium is assumed on endemics using deterministic equation for endemics and immigrants
  #  . eqmodel = 4 : equilibrium is assumed on immigrants using deterministic equation for endemics and immigrants
  #  . eqmodel = 5 : equilibrium is assumed on endemics and immigrants using deterministic equation for endemics and immigrants

  if(!is.list(CS_version)) CS_version <- as.list(CS_version)
  function_to_optimize <- CS_version$function_to_optimize
  if(is.null(function_to_optimize)) function_to_optimize <- 'DAISIE'
  if(function_to_optimize == 'DAISIE_DE') {
    DAISIE_loglik_all_choosepar_fun <- DAISIE_DE_loglik_all_choosepar
  } else
  {
    DAISIE_loglik_all_choosepar_fun <- DAISIE_loglik_all_choosepar
  }

  out2err <- data.frame(
    lambda_c = NA,
    mu = NA,
    K = NA,
    gamma = NA,
    lambda_a = NA,
    loglik = NA,
    df = NA,
    conv = NA
  )
  out2err <- invisible(out2err)
  idparseq <- c()
  if (eqmodel == 1 | eqmodel == 3 | eqmodel == 13) {
    idparseq <- 2
  }
  if (eqmodel == 2 | eqmodel == 4) {
    idparseq <-  4
  }
  if (eqmodel == 5 | eqmodel == 15) {
    idparseq <- c(2, 4)
  }
  namepars <- c(
    "lambda_c",
    "mu",
    "K",
    "gamma",
    "lambda_a",
    "lambda_c2",
    "mu2",
    "K2",
    "gamma2",
    "lambda_a2",
    "prop_type2"
  )
  all_no_shift <- 6:10
  max_idpars <- 11

  if (max(idparsopt, -Inf) <= 6 &&
      max(idparsfix, -Inf) <= 6 &&
      (6 %in% idparsopt || 6 %in% idparsfix)) {
    max_idpars <- 12
    idparsnoshift <- 7:11
    all_no_shift <- 7:11
    namepars <- c(
      "lambda_c",
      "mu",
      "K",
      "gamma",
      "lambda_a",
      "prob_init_pres",
      "lambda_c2",
      "mu2",
      "K2",
      "gamma2",
      "lambda_a2",
      "prop_type2"
    )
    nc <- NA
    names(nc) <- "prob_init_pres"
    out2err <- add_column_to_dataframe(df = out2err,
                                       position = 'lambda_a',
                                       column_to_insert = nc)
  }

  print_ml_par_settings(
    namepars = namepars,
    idparsopt = idparsopt,
    idparsfix = idparsfix,
    idparsnoshift = idparsnoshift,
    all_no_shift = all_no_shift,
    verbose = verbose
  )

  idpars <- sort(c(idparsopt, idparsfix, idparsnoshift, idparseq))

  missnumspec <- unlist(lapply(datalist, function(list) {list$missing_species})) # nolint
  if (max(missnumspec) > (res - 1)) {
    warning(
      "The number of missing species is too large relative to the
       resolution of the ODE.")
    return(out2err)
  }

  if (max(missnumspec) > res/10 && max(missnumspec) <= (res - 1)) {
    warning(
      "The number of missing species is quite high relative to the
        resolution of the ODE.")
  }

  if ((length(idpars) != max(idpars))) {
    warning("The parameters to be optimized and/or fixed are incoherent.")
    return(out2err)
  }

  if ((!all(idpars == 1:max(idpars))) || # nolint
      (length(initparsopt) != length(idparsopt)) ||
      (length(parsfix) != length(idparsfix))) {
    warning("The parameters to be optimized and/or fixed are incoherent.")
    return(out2err)
  }
  if (length(idparseq) == 0) {
  } else {
    if (ddmodel == 3) {
      warning("Equilibrium optimization is not implemented for ddmodel = 3")
    } else {
      message(
        "You are assuming equilibrium. Extinction and/or immigration will
          be considered a function of the other parameters, the species
          pool size, the number of endemics,
          and/or the number of non-endemics"
      )
    }
  }
  trparsopt <- initparsopt / (1 + initparsopt)
  trparsopt[which(initparsopt == Inf)] <- 1
  trparsfix <- parsfix / (1 + parsfix)
  trparsfix[which(parsfix == Inf)] <- 1
  pars2 <- c(
    res,
    ddmodel,
    cond,
    verbose,
    island_ontogeny,
    eqmodel,
    tol,
    maxiter,
    x_E,
    x_I
  )

  optimpars <- c(tol, maxiter)
  initloglik <- DAISIE_loglik_all_choosepar_fun(
    trparsopt = trparsopt,
    trparsfix = trparsfix,
    idparsopt = idparsopt,
    idparsfix = idparsfix,
    idparsnoshift = idparsnoshift,
    idparseq = idparseq,
    pars2 = pars2,
    datalist = datalist,
    methode = methode,
    CS_version = CS_version,
    abstolint = tolint[1],
    reltolint = tolint[2],
    equal_extinction = equal_extinction
  )

  print_init_ll(initloglik = initloglik, verbose = verbose)

  if (initloglik == -Inf) {
    warning(
      "The initial parameter values have a likelihood that is equal to 0 or
       below machine precision. Try again with different initial values."
    )
    return(out2err)
  }

  out <- DDD::optimizer(
    optimmethod = optimmethod,
    optimpars = optimpars,
    fun = DAISIE_loglik_all_choosepar_fun,
    trparsopt = trparsopt,
    idparsopt = idparsopt,
    trparsfix = trparsfix,
    idparsfix = idparsfix,
    idparsnoshift = idparsnoshift,
    idparseq = idparseq,
    pars2 = pars2,
    datalist = datalist,
    methode = methode,
    CS_version = CS_version,
    abstolint = tolint[1],
    reltolint = tolint[2],
    jitter = jitter,
    num_cycles = num_cycles,
    equal_extinction = equal_extinction
  )
  if (out$conv != 0) {
    warning(
      "Optimization has not converged.
       Try again with different initial values.")
    out2 <- out2err
    out2$conv <- out$conv
    return(out2)
  }
  MLtrpars <- as.numeric(unlist(out$par))
  MLpars <- MLtrpars / (1 - MLtrpars)
  ML <- as.numeric(unlist(out$fvalues))

  if (sum(idparsnoshift %in% (all_no_shift)) != 5) {
    MLpars1 <- rep(0, 11)
  } else {
    MLpars1 <- rep(0, 6)
  }
  MLpars1[idparsopt] <- MLpars
  if (length(idparsfix) != 0) {
    MLpars1[idparsfix] <- parsfix
  }

  if (eqmodel > 0) {
    MLpars1 <- DAISIE_eq(datalist, MLpars1, pars2[-5])
  }

  if (MLpars1[3] > 10 ^ 7) {
    MLpars1[3] <- Inf
  }

  if (sum(idparsnoshift %in% (all_no_shift)) != 5) {
    if (length(idparsnoshift) != 0) {
      MLpars1[idparsnoshift] <- MLpars1[idparsnoshift - 5]
    }
    if (MLpars1[8] > 10 ^ 7) {
      MLpars1[8] <- Inf
    }
    out2 <- data.frame(
      lambda_c = MLpars1[1],
      mu = MLpars1[2],
      K = MLpars1[3],
      gamma = MLpars1[4],
      lambda_a = MLpars1[5],
      lambda_c2 = MLpars1[6],
      mu2 = MLpars1[7],
      K2 = MLpars1[8],
      gamma2 = MLpars1[9],
      lambda_a2 = MLpars1[10],
      prop_type2 = MLpars1[11],
      loglik = ML,
      df = length(initparsopt),
      conv = unlist(out$conv)
    )
    pars_to_print <- MLpars1[1:11]
    parnames <- c('lambda^c','mu','K','gamma','lambda^a','lambda^c2','mu2','K2','gamma2','lambda^a2','prop_type2')
  } else if (all(all_no_shift == 7:11)) {
    out2 <- data.frame(
      lambda_c = MLpars1[1],
      mu = MLpars1[2],
      K = MLpars1[3],
      gamma = MLpars1[4],
      lambda_a = MLpars1[5],
      prob_init_pres = MLpars1[6],
      loglik = ML,
      df = length(initparsopt),
      conv = unlist(out$conv)
    )
    pars_to_print <- MLpars1[1:6]
    parnames <- c("lambda^c", "mu", "K", "gamma", "lambda^a", "prob_init_pres")
  } else {
    out2 <- data.frame(
      lambda_c = MLpars1[1],
      mu = MLpars1[2],
      K = MLpars1[3],
      gamma = MLpars1[4],
      lambda_a = MLpars1[5],
      loglik = ML,
      df = length(initparsopt),
      conv = unlist(out$conv)
    )
    pars_to_print <- MLpars1[1:5]
    parnames <- c('lambda^c','mu','K','gamma','lambda^a')
  }
  if(function_to_optimize != 'DAISIE') {
    parnames[which(parnames == 'mu')] <- 'mu_E'
    parnames[which(parnames == 'K')] <- 'mu_NE'
    parnames[which(parnames == 'mu2')] <- 'mu2_E'
    parnames[which(parnames == 'K2')] <- 'mu2_NE'
  }
  print_parameters_and_loglik(pars = pars_to_print,
                              loglik = ML,
                              verbose = verbose,
                              parnames = parnames,
                              type = 'island_ML')
  if (eqmodel > 0) {
    M <- calcMN(datalist, MLpars1)
    ExpEIN <- DAISIE_ExpEIN(datalist[[1]]$island_age, MLpars1, M) # nolint start
    message(
      paste0("The expected number of endemics, non-endemics, and the total at ",
        "these parameters is: "),
      paste(ExpEIN[[1]], ExpEIN[[2]], ExpEIN[[3]])
    ) # nolint end
  }
  if(function_to_optimize != 'DAISIE') {
    names(out2[2:3]) <- c('mu_E','mu_NE')
  }
  return(invisible(out2))
}
