#' Integrates the loglikelihood of a single clade across a parameter weighted by a
#' given distribution
#'
#' @inheritParams DAISIE_loglik_CS
#' @param CS_version a list with the following elements:
#' \itemize{
#'   \item{choice: the choice of loglikelihood - this should be 2}
#'   \item{pick_parameter: the parameter to integrate over. This is \code{'lambda^c'},
#'    \code{'mu'},\code{'K'},\code{'gamma'} or \code{'lambda^a'}}
#'   \item{distribution: distribution to weigh the likelihood by; either \code{'lognormal'}
#'   or \code{'gamma'}}
#'   \item{sd_par: standard deviation of the distribution}
#'   \item{number_of_cores: number_of_cores to use in the integration}
#'   }
#' @return A loglikelihood value
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
  verbose) {
  testit::assert(is.list(CS_version))
  cpus <- CS_version$number_of_cores
  if (CS_version$distribution == 'lognormal') {
  ##if (sign(CS_version) == -1) {
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

  sd_par <- CS_version$sd_par
  pick <- which(c('lambda^c','mu','K','gamma','lambda^a') == CS_version$pick_parameter)
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
                                                 sd_par,
                                                 cpus) {
    if(cpus > 1) {
      X <- NULL; rm(X)
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
    sd_par = sd_par,
    cpus = cpus)$value)

    return(integrated_loglik)
}

initiate_cluster <- function(cpus = 1, cl = NULL)
{
  cpus <- min(cpus,parallel::detectCores())
  if(!is.null(cl)) try(parallel::stopCluster(cl))
  cl <- parallel::makeCluster(cpus)
  doParallel::registerDoParallel(cl, cores = cpus)
  #on.exit(parallel::stopCluster(cl))
  return(cl)
}
