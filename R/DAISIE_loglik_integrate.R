#' Integrates the loglikelihood of a single clade across a parameter weighted by a
#' given distribution
#'
#' @inheritParams DAISIE_loglik_CS
#' @param CS_version a list with the following elements:
#' \itemize{
#'   \item{model: the CS model to run, options are \code{1} for single rate
#'   DAISIE model, \code{2} for multi-rate DAISIE, or \code{0} for IW test
#'   model}
#'   \item{pick_parameter: the parameter to relax (integrate over). Options are
#' \code{"cladogenesis"}, \code{"extinction"}, \code{"carrying_capacity"},
#' \code{"immigration"}, or \code{"anagenesis"}}
#'   \item{distribution: the distribution to weigh the likelihood, either
#' \code{"lognormal"} or \code{"gamma"}}
#'   \item{sd: standard deviation of the distribution}
#'   \item{num_cores: number of cores to use in the integration}
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
  if (CS_version$distribution == "lognormal") {
    rho <- function(DAISIE_par, DAISIE_dist_pars) {
      sigma_squared <- log(1 + (DAISIE_dist_pars[2] / DAISIE_dist_pars[1])^2)
      return(stats::dlnorm(x = DAISIE_par,
                    meanlog = log(DAISIE_dist_pars[1]) - sigma_squared / 2,
                    sdlog = sqrt(sigma_squared)))
    }
  }
  if (CS_version$distribution == "gamma") {
    rho <- function(DAISIE_par, DAISIE_dist_pars) {
      return(stats::dgamma(x = DAISIE_par,
                           shape = DAISIE_dist_pars[1]^2 / DAISIE_dist_pars[2]^2,
                           scale = DAISIE_dist_pars[2]^2 / DAISIE_dist_pars[1]))
    }
  }

  sd <- CS_version$sd
  pick <- which(c("cladogenesis",
                  "extinction",
                  "carrying_capacity",
                  "immigration",
                  "anagenesis") == CS_version$pick_parameter)
  mean <- pars1[pick]

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
                                      mean,
                                      sd) {
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
        DAISIE_dist_pars = c(mean, sd)
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
                                                 mean,
                                                 sd,
                                                 cpus) {
    if(cpus > 1) {
      X <- NULL; rm(X)
      loglik_vec <- foreach::foreach(
        X = DAISIE_par_vec,
        .combine = c,
        .export = c("DAISIE_loglik_integrand","rho"),
        .packages = c("DAISIE", "foreach", "deSolve", "doParallel"),
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
                                mean = mean,
                                sd = sd)
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
                                                 mean = mean,
                                                 sd = sd)
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
    mean = mean,
    sd = sd,
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
