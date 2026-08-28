
#' @export
# Don't document this function. For internal use only.
DAISIE_loglik_all_choosepar_time <- function(
    trparsopt,
    trparsfix,
    idparsopt,
    idparsfix,
    idparsnoshift = NULL,
    idparseq = NULL,
    pars2,
    datalist,
    methode,
    CS_version = list(
      model = 1,
      function_to_optimize = "DAISIE",
      sampling = "n"
    ),
    abstolint = 1E-16,
    reltolint = 1E-10,
    equal_extinction = TRUE
) {

  # Time-dependent pars1 has 14 parameters:
  # 1  lambda_c
  # 2  mu
  # 3  K
  # 4  gamma
  # 5  lambda_a
  # 6  d
  # 7  x
  # 8  max_area
  # 9  current_area
  # 10 proportional_peak_t
  # 11 total_island_age
  # 12 sea_level_amplitude
  # 13 sea_level_frequency
  # 14 island_gradient_angle

  trpars1 <- rep(0, 14)

  trpars1[idparsopt] <- trparsopt

  if (length(idparsfix) != 0) {
    trpars1[idparsfix] <- trparsfix
  }

  if (max(trpars1) > 1 || min(trpars1) < 0) {
    return(-Inf)
  }

  pars1 <- trpars1 / (1 - trpars1)
  pars1[which(trpars1 == 1)] <- Inf

  # Basic parameter checks
  if (
    any(is.na(pars1)) ||
    any(is.nan(pars1)) ||
    any(pars1[1:14] < 0)
  ) {
    return(-Inf)
  }

  # Island age used for the likelihood
  total_time <- datalist[[1]]$island_age

  if (is.null(total_time) || is.na(total_time)) {
    total_time <- max(unlist(lapply(
      datalist,
      function(x) max(abs(as.numeric(x$branching_times)), na.rm = TRUE)
    )), na.rm = TRUE)
  }

  # The island age used in the data should not exceed the total island age
  # in the area parameters.
  if (total_time > pars1[11]) {
    return(-Inf)
  }

  area_pars <- create_area_pars(
    max_area = pars1[8],
    current_area = pars1[9],
    proportional_peak_t = pars1[10],
    total_island_age = pars1[11],
    sea_level_amplitude = pars1[12],
    sea_level_frequency = pars1[13],
    island_gradient_angle = pars1[14]
  )

  peak <- tryCatch(
    calc_peak(
      total_time = total_time,
      area_pars = area_pars
    ),
    error = function(e) NA
  )

  if (is.na(peak) || is.nan(peak)) {
    return(-Inf)
  }

  island_ontogeny <- 0
  sea_level <- 0

  if (length(pars2) >= 5 && !is.na(pars2[5])) {
    island_ontogeny <- pars2[5]
  }

  if (length(pars2) >= 6 && !is.na(pars2[6])) {
    sea_level <- pars2[6]
  }

  pars1_time <- c(
    pars1,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    total_time = total_time,
    peak = peak
  )

  loglik <- DAISIE_loglik_all(
    pars1 = pars1_time,
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


#' @export
DAISIE_ML5 <- function(
    datalist,
    initparsopt,
    idparsopt,
    parsfix,
    idparsfix,
    res = 100,
    ddmodel = 0,
    cond = 0,
    eqmodel = 0,
    x_E = 0.95,
    x_I = 0.98,
    tol = c(1E-4, 1E-5, 1E-7),
    maxiter = 1000 * round((1.25) ^ length(idparsopt)),
    methode = "lsodes",
    optimmethod = "simplex",
    CS_version = list(
      model = 1,
      function_to_optimize = "DAISIE",
      sampling = "n"
    ),
    verbose = 0,
    tolint = c(1E-16, 1E-10),
    island_ontogeny = NA,
    sea_level = NA,
    jitter = 0,
    num_cycles = 1,
    equal_extinction = TRUE
) {

  if (!is.list(CS_version)) {
    CS_version <- as.list(CS_version)
  }

  function_to_optimize <- CS_version$function_to_optimize
  if (is.null(function_to_optimize)) {
    function_to_optimize <- "DAISIE"
  }

  if (function_to_optimize == "DAISIE_DE") {
    stop("DAISIE_DE is not implemented in this time-dependent ML1 wrapper.")
  }

  if (eqmodel != 0) {
    warning("Equilibrium optimization is not implemented for the time-dependent model.")
    return(invisible(data.frame(
      lambda_c = NA,
      mu = NA,
      K = NA,
      gamma = NA,
      lambda_a = NA,
      d = NA,
      x = NA,
      max_area = NA,
      current_area = NA,
      proportional_peak_t = NA,
      total_island_age = NA,
      sea_level_amplitude = NA,
      sea_level_frequency = NA,
      island_gradient_angle = NA,
      loglik = NA,
      df = NA,
      conv = NA
    )))
  }

  out2err <- invisible(data.frame(
    lambda_c = NA,
    mu = NA,
    K = NA,
    gamma = NA,
    lambda_a = NA,
    d = NA,
    x = NA,
    max_area = NA,
    current_area = NA,
    proportional_peak_t = NA,
    total_island_age = NA,
    sea_level_amplitude = NA,
    sea_level_frequency = NA,
    island_gradient_angle = NA,
    loglik = NA,
    df = NA,
    conv = NA
  ))

  namepars <- c(
    "lambda_c",
    "mu",
    "K",
    "gamma",
    "lambda_a",
    "d",
    "x",
    "max_area",
    "current_area",
    "proportional_peak_t",
    "total_island_age",
    "sea_level_amplitude",
    "sea_level_frequency",
    "island_gradient_angle"
  )

  max_idpars <- 14

  idparsnoshift <- integer(0)
  idparseq <- integer(0)
  all_no_shift <- integer(0)

  print_ml_par_settings(
    namepars = namepars,
    idparsopt = idparsopt,
    idparsfix = idparsfix,
    idparsnoshift = idparsnoshift,
    all_no_shift = all_no_shift,
    verbose = verbose
  )

  idpars <- sort(c(idparsopt, idparsfix))

  if (max(idpars, -Inf) > max_idpars) {
    warning("For DAISIE_ML1_time, parameter IDs must be between 1 and 14.")
    return(out2err)
  }

  missnumspec <- unlist(lapply(datalist, function(list) {
    list$missing_species
  }))

  if (max(missnumspec) > (res - 1)) {
    warning(
      "The number of missing species is too large relative to the resolution of the ODE."
    )
    return(out2err)
  }

  if (max(missnumspec) > res / 10 && max(missnumspec) <= (res - 1)) {
    warning(
      "The number of missing species is quite high relative to the resolution of the ODE."
    )
  }

  if (length(idpars) != max(idpars)) {
    warning("The parameters to be optimized and/or fixed are incoherent.")
    return(out2err)
  }

  if (
    !all(idpars == 1:max(idpars)) ||
    length(initparsopt) != length(idparsopt) ||
    length(parsfix) != length(idparsfix)
  ) {
    warning("The parameters to be optimized and/or fixed are incoherent.")
    return(out2err)
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
    sea_level,
    eqmodel,
    tol,
    maxiter,
    x_E,
    x_I
  )

  optimpars <- c(tol, maxiter)

  initloglik <- DAISIE_loglik_all_choosepar_time(
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
      "The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values."
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
      "Optimization has not converged. Try again with different initial values."
    )
    out2 <- out2err
    out2$conv <- out$conv
    return(out2)
  }

  MLtrpars <- as.numeric(unlist(out$par))
  MLpars <- MLtrpars / (1 - MLtrpars)
  MLpars[which(MLtrpars == 1)] <- Inf

  ML <- as.numeric(unlist(out$fvalues))

  MLpars1 <- rep(NA, 14)

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
    max_area = MLpars1[8],
    current_area = MLpars1[9],
    proportional_peak_t = MLpars1[10],
    total_island_age = MLpars1[11],
    sea_level_amplitude = MLpars1[12],
    sea_level_frequency = MLpars1[13],
    island_gradient_angle = MLpars1[14],
    loglik = ML,
    df = length(initparsopt),
    conv = unlist(out$conv)
  )

  pars_to_print <- MLpars1[1:14]

  parnames <- c(
    "lambda^c",
    "mu",
    "K",
    "gamma",
    "lambda^a",
    "d",
    "x",
    "max_area",
    "current_area",
    "proportional_peak_t",
    "total_island_age",
    "sea_level_amplitude",
    "sea_level_frequency",
    "island_gradient_angle"
  )

  print_parameters_and_loglik(
    pars = pars_to_print,
    loglik = ML,
    verbose = verbose,
    parnames = parnames,
    type = "island_ML_time"
  )

  return(invisible(out2))
}

