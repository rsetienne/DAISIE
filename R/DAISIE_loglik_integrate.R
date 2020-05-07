DAISIE_loglik_integrate <- function(
  pars1,
  pars2,
  brts,
  stac,
  missnumspec,
  CS_version,
  methode,
  abstolint,
  reltolint,
  verbose
) {
  # DaVinci code for CS_version:
  # The sign determines whether the lognormal (negative) or the gamma
  # distribution (positive) is chosen.
  # The number before the decimal separator determines what DAISIE parameter
  # is chosen where counting starts at 2. For example, CS_version = 2.x
  # corresponds to the 1st parameter, lambda^c, whereas CS_version = 4.x
  # corresponds to the 3rd parameter, K.
  # The number after the decimal separator divided by 10 sets the standard
  # deviation for the lognormal and the gamma.
  # For example CS_version = -x.34 sets the sd of the lognormal to 3.4.
  # The sd is thus limited to values below 10.
  # A value of CS_version = +4.31 thus indicates a gamma distribution (+) for
  # the 3rd parameter (4 - 1) with a standard deviation of 3.1.
  if (sign(CS_version) == -1) {
    rho <- function(DAISIE_par, DAISIE_dist_pars) {
      sigma_squared <- log(1 + (DAISIE_dist_pars[2] / DAISIE_dist_pars[1])^2)
      return(stats::dlnorm(x = DAISIE_par,
                    meanlog = log(DAISIE_dist_pars[1]) - sigma_squared / 2,
                    sdlog = sqrt(sigma_squared)))
    }
  } else {
    rho <- function(DAISIE_par, DAISIE_dist_pars) {
      return(stats::dgamma(x = DAISIE_par,
                           shape = DAISIE_dist_pars[1]^2 / DAISIE_dist_pars[2]^2,
                           scale = DAISIE_dist_pars[2]^2 / DAISIE_dist_pars[1]))
    }
  }

  sd_par <- 10 * (abs(CS_version) - floor(abs(CS_version)))
  pick <- floor(abs(CS_version)) - 1
  mean_par <- pars1[pick]

  DAISIE_loglik_integrand <- function(DAISIE_par,
                                      pars1,
                                      pars2,
                                      brts,
                                      stac,
                                      missnumspec,
                                      methode,
                                      abstolint,
                                      reltolint,
                                      verbose,
                                      pick,
                                      mean_par,
                                      sd_par) {
    pars1[pick] <- DAISIE_par
    loglik_DAISIE_par <- exp(DAISIE_loglik(
      pars1 = pars1,
      pars2 = pars2,
      brts = brts,
      stac = stac,
      missnumspec = missnumspec,
      methode = methode,
      abstolint = abstolint,
      reltolint = reltolint,
      verbose = verbose)
      ) *
      rho(
        DAISIE_par = DAISIE_par,
        DAISIE_dist_pars = c(mean_par, sd_par)
      )
    return(loglik_DAISIE_par)
  }

  DAISIE_loglik_integrand_vectorized <- function(DAISIE_par_vec,
                                                 pars1,
                                                 pars2,
                                                 brts,
                                                 stac,
                                                 missnumspec,
                                                 methode,
                                                 abstolint,
                                                 reltolint,
                                                 verbose,
                                                 pick,
                                                 mean_par,
                                                 sd_par) {
    cpus <- 1
    cpus <- min(cpus, parallel::detectCores())
    if (cpus > 1) {
      if (.Platform$OS.type == "windows") {
        cl <- parallel::makeCluster(cpus - 1)
        doParallel::registerDoParallel(cl)
        on.exit(parallel::stopCluster(cl))
      } else {
        doMC::registerDoMC(cpus - 1)
      }
      X <- NULL; rm(X)
      loglik_vec <- rep(NA,length(DAISIE_par_vec))
      loglik_vec <- foreach::foreach(
        X = DAISIE_par_vec,
        .combine = c,
        .export = c("DAISIE_loglik_integrand","rho"),
        .packages = c('DAISIE','foreach','deSolve','doParallel'),
        .verbose = FALSE) %dopar%
        DAISIE_loglik_integrand(DAISIE_par = X,
                                pars1 = pars1,
                                pars2 = pars2,
                                brts = brts,
                                stac = stac,
                                missnumspec = missnumspec,
                                methode = methode,
                                abstolint = abstolint,
                                reltolint = reltolint,
                                verbose = verbose,
                                pick = pick,
                                mean_par = mean_par,
                                sd_par = sd_par)
      parallel::stopCluster(cl)
    } else {
      loglik_vec <- rep(NA,length(DAISIE_par_vec))
      for(i in 1:length(DAISIE_par_vec)) {
        loglik_vec[i] <- DAISIE_loglik_integrand(DAISIE_par = DAISIE_par_vec[i],
                                                 pars1 = pars1,
                                                 pars2 = pars2,
                                                 brts = brts,
                                                 stac = stac,
                                                 missnumspec = missnumspec,
                                                 methode = methode,
                                                 abstolint = abstolint,
                                                 reltolint = reltolint,
                                                 verbose = verbose,
                                                 pick = pick,
                                                 mean_par = mean_par,
                                                 sd_par = sd_par)
      }
    }
    return(loglik_vec)
  }

  integrated_loglik <- log(stats::integrate(
    f = DAISIE_loglik_integrand_vectorized,
    lower = 0,
    upper = Inf,
    pars1 = pars1,
    pars2 = pars2,
    brts = brts,
    stac = stac,
    missnumspec = missnumspec,
    methode = methode,
    abstolint = abstolint,
    reltolint = reltolint,
    verbose = verbose,
    pick = pick,
    mean_par = mean_par,
    sd_par = sd_par)$value)

  return(integrated_loglik)
}

#' Vectorize a function and allow multithreading
#'
#' @inheritParams base::Vectorize
#' @inheritParams parallel::mcmapply
#' @param mc.cores Number of cores to use
#' @return A function with the same arguments as FUN wrapping
#' a call to \code{mcmapply}
mcVectorize <- function(FUN,
                        vectorize.args = arg.names,
                        SIMPLIFY = TRUE,
                        USE.NAMES = TRUE,
                        mc.preschedule = TRUE,
                        mc.set.seed = TRUE,
                        mc.silent = FALSE,
                        mc.cores = getOption('mc.cores',2L),
                        mc.cleanup = TRUE)
{
  arg.names <- as.list(formals(FUN))
  arg.names[["..."]] <- NULL
  arg.names <- names(arg.names)
  vectorize.args <- as.character(vectorize.args)
  if (!length(vectorize.args))
    return(FUN)
  if (!all(vectorize.args %in% arg.names))
    stop("must specify names of formal arguments for 'vectorize'")
  collisions <- arg.names %in% c("FUN", "SIMPLIFY",
                                 "USE.NAMES", "vectorize.args")
  if (any(collisions))
    stop(sQuote("FUN"), " may not have argument(s) named ",
         paste(sQuote(arg.names[collisions]), collapse = ", "))
  FUNV <- function() {
    args <- lapply(as.list(match.call())[-1L], eval, parent.frame())
    names <- if (is.null(names(args)))
      character(length(args))
    else names(args)
    dovec <- names %in% vectorize.args
    if(Sys.info()['sysname'] == 'Windows') {
      if(mc.cores > 1) {
        future::plan(multiprocess)
        do.call(future.apply::future_mapply, c(FUN = FUN,
                                               args[dovec],
                                               MoreArgs = list(args[!dovec]),
                                               SIMPLIFY = SIMPLIFY,
                                               USE.NAMES = USE.NAMES,
                                               future.scheduling = Inf))
      } else {
        do.call(mapply, c(FUN = FUN,
                          args[dovec],
                          MoreArgs = list(args[!dovec]),
                          SIMPLIFY = SIMPLIFY,
                          USE.NAMES = USE.NAMES))      }
    } else {
      do.call(parallel::mcmapply, c(FUN = FUN,
                        args[dovec],
                        MoreArgs = list(args[!dovec]),
                        SIMPLIFY = SIMPLIFY,
                        USE.NAMES = USE.NAMES,
                        mc.preschedule = mc.preschedule,
                        mc.set.seed = mc.set.seed,
                        mc.silent = mc.silent,
                        mc.cores = mc.cores,
                        mc.cleanup = mc.cleanup))
    }
  }
  formals(FUNV) <- formals(FUN)
  FUNV
}
