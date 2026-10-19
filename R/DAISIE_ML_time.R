# Don't document this function. For internal use only.
# Transforms the parameters from (0, 1) to (0, Inf), adds the fixed area
# parameters and computes the loglikelihood of the time-dependent CS model.
DAISIE_loglik_all_choosepar_time <- function(trparsopt,
                                             trparsfix,
                                             idparsopt,
                                             idparsfix,
                                             pars2,
                                             area_pars,
                                             datalist,
                                             methode,
                                             CS_version = list(model = 1, function_to_optimize = 'DAISIE', sampling = 'n'),
                                             abstolint = 1E-16,
                                             reltolint = 1E-10) {
  trpars1 <- rep(0, 7)
  trpars1[idparsopt] <- trparsopt
  if (length(idparsfix) != 0) {
    trpars1[idparsfix] <- trparsfix
  }
  if (max(trpars1) > 1 || min(trpars1) < 0) {
    return(-Inf)
  }
  pars1 <- trpars1 / (1 - trpars1)
  loglik <- DAISIE_loglik_all(
    pars1 = c(pars1, area_pars),
    pars2 = pars2,
    datalist = datalist,
    methode = methode,
    CS_version = CS_version,
    abstolint = abstolint,
    reltolint = reltolint
  )
  if (is.nan(loglik) || is.na(loglik)) {
    message("There are parameter values used which cause numerical problems.")
    loglik <- -Inf
  }
  return(loglik)
}

#' Maximum likelihood estimation for the time-dependent CS model
#'
#' Estimates parameters of the clade-specific (CS) model with island area
#' that changes through time (island ontogeny and/or sea level), for a single
#' type of species. The likelihood is that of \code{\link{DAISIE_loglik_all}}
#' with 14 model parameters. This function is called by
#' \code{\link{DAISIE_ML}} when \code{island_ontogeny} is not \code{NA}.
#'
#' The model parameters that can be estimated or fixed are, in this order:
#' 1 = cladogenesis rate, 2 = extinction rate, 3 = carrying capacity,
#' 4 = immigration rate, 5 = anagenesis rate, 6 = \code{d}, the exponent of
#' the area dependence of cladogenesis, 7 = \code{x}, the exponent of the area
#' dependence of extinction. All seven must be either in \code{idparsopt} or in
#' \code{idparsfix}. The parameters are assumed to be non-negative. Diversity
#' dependence is linear in cladogenesis and immigration (\code{ddmodel = 11}),
#' which is what the time-dependent likelihood implements.
#'
#' The area parameters in \code{area_pars} are not estimated.
#'
#' @inheritParams default_params_doc
#' @param initparsopt initial values of the parameters to be optimized
#' @param idparsopt ids (1 to 7) of the parameters to be optimized
#' @param parsfix values of the parameters that are fixed
#' @param idparsfix ids (1 to 7) of the parameters that are fixed
#' @param area_pars numeric vector of length 7 (or the list returned by
#' \code{\link{create_area_pars}}): maximum area, current area,
#' proportion of the island's lifetime at which area peaks, total island age,
#' sea level amplitude, sea level frequency and island gradient angle (see
#' \code{\link{create_area_pars}}). The total island age must exceed the age of
#' the island in the data when \code{island_ontogeny = 1}.
#' @param island_ontogeny 1 or "beta" for beta-function island ontogeny,
#' 0 or "const" for constant ontogeny.
#' @param sea_level 1 or "sine" for sinusoidal sea level, 0 or "const" for
#' constant sea level.
#' @param res resolution of the ODE (number of state-vector elements per
#' block)
#' @param cond conditioning: 0 for none, 1 for conditioning on presence of
#' at least one species
#' @param tol tolerances of the optimization
#' @param maxiter maximum number of iterations
#' @param methode ODE solver; an odeint stepper such as
#' \code{"odeint::runge_kutta_dopri5"} or a deSolve method such as
#' \code{"lsodes"}
#' @param optimmethod optimization method, see \code{\link[DDD]{optimizer}}
#' @param tolint tolerances of the ODE solver (absolute, relative)
#' @param jitter amount of jitter on the initial parameters
#' @param num_cycles number of cycles of the optimization
#' @param verbose verbosity
#' @return A data frame with the estimated parameters (\code{lambda_c},
#' \code{mu}, \code{K}, \code{gamma}, \code{lambda_a}, \code{d}, \code{x}),
#' the maximum loglikelihood \code{loglik}, the number of estimated parameters
#' \code{df} and the convergence code \code{conv} (0 means convergence).
#' @keywords internal
DAISIE_ML_time <- function(datalist,
                           initparsopt,
                           idparsopt,
                           parsfix,
                           idparsfix,
                           area_pars,
                           island_ontogeny = 1,
                           sea_level = 0,
                           res = 100,
                           cond = 0,
                           tol = c(1E-4, 1E-5, 1E-7),
                           maxiter = 1000 * round((1.25) ^ length(idparsopt)),
                           methode = "odeint::runge_kutta_dopri5",
                           optimmethod = "simplex",
                           CS_version = list(model = 1, function_to_optimize = 'DAISIE', sampling = 'n'),
                           verbose = 0,
                           tolint = c(1E-16, 1E-10),
                           jitter = 0,
                           num_cycles = 1) {
  if (identical(island_ontogeny, "beta")) island_ontogeny <- 1
  if (identical(island_ontogeny, "const")) island_ontogeny <- 0
  if (identical(sea_level, "sine")) sea_level <- 1
  if (identical(sea_level, "const")) sea_level <- 0
  if (!isTRUE(island_ontogeny %in% c(0, 1))) {
    stop("island_ontogeny must be 0 ('const') or 1 ('beta').")
  }
  if (!isTRUE(sea_level %in% c(0, 1))) {
    stop("sea_level must be 0 ('const') or 1 ('sine').")
  }
  if (is.list(area_pars)) {
    # as created by create_area_pars()
    area_pars <- unlist(area_pars[c(
      "max_area", "current_area", "proportional_peak_t", "total_island_age",
      "sea_level_amplitude", "sea_level_frequency", "island_gradient_angle"
    )])
  }
  if (length(area_pars) != 7 || !is.numeric(area_pars)) {
    stop("area_pars must be a numeric vector of length 7.")
  }
  if (!is.null(CS_version) && length(CS_version) > 1 &&
      isTRUE(CS_version$function_to_optimize == "DAISIE_DE")) {
    stop("Time dependent estimation is not available for DAISIE_DE.")
  }
  idpars <- sort(c(idparsopt, idparsfix))
  if (length(idpars) != 7 || !all(idpars == 1:7) ||
      length(initparsopt) != length(idparsopt) ||
      length(parsfix) != length(idparsfix)) {
    stop(
      "For the time-dependent model all parameters 1 to 7 (lambda^c, mu, K, ",
      "gamma, lambda^a, d, x) must be specified exactly once in idparsopt or ",
      "idparsfix, with matching initial values or fixed values."
    )
  }
  if (min(c(initparsopt, parsfix)) < 0) {
    stop("The parameters must be non-negative.")
  }
  if (length(datalist[[1]]$island_age) != 1 ||
      (island_ontogeny == 1 && datalist[[1]]$island_age >= area_pars[4])) {
    stop(
      "The island age in the data must be smaller than the total island age ",
      "in area_pars (area_pars[4])."
    )
  }
  if (jitter > 0) {
    initparsopt <- initparsopt * exp(stats::runif(length(initparsopt), -jitter, jitter))
  }
  pars2 <- c(res, 11, cond, verbose, island_ontogeny, sea_level)
  trparsopt <- initparsopt / (1 + initparsopt)
  trparsopt[which(initparsopt == Inf)] <- 1
  trparsfix <- parsfix / (1 + parsfix)
  trparsfix[which(parsfix == Inf)] <- 1
  optimpars <- c(tol, maxiter)

  out2err <- data.frame(
    lambda_c = NA, mu = NA, K = NA, gamma = NA, lambda_a = NA, d = NA, x = NA,
    loglik = NA, df = NA, conv = NA
  )
  initloglik <- DAISIE_loglik_all_choosepar_time(
    trparsopt = trparsopt, trparsfix = trparsfix, idparsopt = idparsopt,
    idparsfix = idparsfix, pars2 = pars2, area_pars = area_pars,
    datalist = datalist, methode = methode, CS_version = CS_version,
    abstolint = tolint[1], reltolint = tolint[2]
  )
  print_init_ll(initloglik = initloglik, verbose = verbose)
  if (initloglik == -Inf) {
    warning(
      "The initial parameter values have a likelihood that is equal to 0 or ",
      "below machine precision. Try again with different initial values."
    )
    return(out2err)
  }
  out <- DDD::optimizer(
    optimmethod = optimmethod,
    optimpars = optimpars,
    fun = DAISIE_loglik_all_choosepar_time,
    trparsopt = trparsopt,
    idparsopt = idparsopt,
    trparsfix = trparsfix,
    idparsfix = idparsfix,
    pars2 = pars2,
    area_pars = area_pars,
    datalist = datalist,
    methode = methode,
    CS_version = CS_version,
    abstolint = tolint[1],
    reltolint = tolint[2],
    jitter = jitter,
    num_cycles = num_cycles
  )
  if (out$conv != 0) {
    warning(
      "Optimization has not converged. ",
      "Try again with different initial values."
    )
    out2err$conv <- out$conv
    return(out2err)
  }
  MLtrpars <- as.numeric(unlist(out$par))
  MLpars <- MLtrpars / (1 - MLtrpars)
  ML <- as.numeric(unlist(out$fvalues))
  MLpars1 <- rep(0, 7)
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
    d = MLpars1[6],
    x = MLpars1[7],
    loglik = ML,
    df = length(initparsopt),
    conv = unlist(out$conv)
  )
  print_parameters_and_loglik(
    pars = MLpars1,
    loglik = ML,
    verbose = verbose,
    parnames = c("lambda^c", "mu", "K", "gamma", "lambda^a", "d", "x"),
    type = "island_ML"
  )
  return(out2)
}
