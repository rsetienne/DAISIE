#' Integrates the loglikelihood of a single clade across a parameter weighted
#' by a given distribution
#'
#' @inheritParams default_params_doc
#'
#' @keywords internal
#'
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
  par_sd <- CS_version$par_sd
  par_upper_bound <- CS_version$par_upper_bound
  pick <- which(c("cladogenesis",
                  "extinction",
                  "carrying_capacity",
                  "immigration",
                  "anagenesis") == CS_version$relaxed_par)
  par_mean <- pars1[pick]

  if (is.infinite(par_mean) || is.infinite(par_sd)) {
    stop("The relaxed parameter mean or standard deviation is infinite")
  }

  if(is.null(CS_version$integration_method)) CS_version$integration_method <- 'standard'
  if(CS_version$integration_method == 'MC') {
    if(is.null(CS_version$seed)) CS_version$seed <- 42
    if(is.null(CS_version$sample_size)) CS_version$sample_size <- 1000
    set.seed(CS_version$seed)
    gamma_pars <- transform_gamma_pars(
      par_mean = par_mean,
      par_sd = par_sd)
    DAISIE_par <- rgamma(n = CS_version$sample_size,shape = gamma_pars$shape, scale = gamma_pars$scale)
    DAISIE_loglik_MC_Vectorized <- Vectorize(DAISIE_loglik_MC, vectorize.args = "DAISIE_par")
    integrated_loglik <- -log(CS_version$sample_size) +
      log(sum(DAISIE_loglik_MC_Vectorized(DAISIE_par = DAISIE_par,
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
                                          par_mean = par_mean,
                                          par_sd = par_sd)))
  } else {
    integrated_loglik <- integral_peak(
      logfun = Vectorize(DAISIE_loglik_integrand,
                         vectorize.args = "DAISIE_par"),
      xx = sort(c(seq(-20, min(20,log(par_upper_bound)), 2),
                  seq(log(par_mean) - 1, log(par_mean) + 1),
                  log((par_mean + 10 * par_sd) / par_mean))),
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
      par_mean = par_mean,
      par_sd = par_sd,
      par_upper_bound = par_upper_bound) -
      cum_rho(par_upper_bound = par_upper_bound,
              DAISIE_dist_pars = list(par_mean = par_mean,
                                      par_sd = par_sd)
      )
  }
  return(integrated_loglik)
}


#' function to calculate the log likelihood for the relaxed rate model.
#'
#' @inheritParams default_params_doc
#'
#' @keywords internal
#'
#' @return A numeric
DAISIE_loglik_MC <- function(DAISIE_par,
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
                             par_mean,
                             par_sd) {
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
    verbose = verbose))
  return(loglik_DAISIE_par)
}


#' Integrand to be integrated to calculate the log likelihood for the relaxed
#' rate model.
#'
#' @inheritParams default_params_doc
#'
#' @keywords internal
#'
#' @return A numeric
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
                                    par_mean,
                                    par_sd) {
  pars1[pick] <- DAISIE_par
  loglik_DAISIE_par <- DAISIE_loglik(
    pars1 = pars1,
    pars2 = pars2,
    brts = brts,
    stac = stac,
    missnumspec = missnumspec,
    methode = methode,
    abstolint = abstolint,
    reltolint = reltolint,
    verbose = verbose) +
    rho(
      DAISIE_par = DAISIE_par,
      DAISIE_dist_pars = list(par_mean = par_mean,
                              par_sd = par_sd)
    )
  return(loglik_DAISIE_par)
}

#' Gamma distribution density parameterised with mean and standard deviation
#'
#' @inheritParams default_params_doc
#'
#' @return Numeric
#' @keywords internal
rho <- function(DAISIE_par, DAISIE_dist_pars) {

  gamma_pars <- transform_gamma_pars(
    par_mean = DAISIE_dist_pars$par_mean,
    par_sd = DAISIE_dist_pars$par_sd)

  gamma_den <- stats::dgamma(
    x = DAISIE_par,
    shape = gamma_pars$shape,
    scale = gamma_pars$scale,
    log = TRUE)

  return(gamma_den)
}

#' Cumulative Gamma distribution parameterised with mean and standard deviation
#'
#' @inheritParams default_params_doc
#'
#' @return Numeric
#' @keywords internal
cum_rho <- function(par_upper_bound, DAISIE_dist_pars) {

  gamma_pars <- transform_gamma_pars(
    par_mean = DAISIE_dist_pars$par_mean,
    par_sd = DAISIE_dist_pars$par_sd)

  gamma_cum_prob <- stats::pgamma(
    q = par_upper_bound,
    shape = gamma_pars$shape,
    scale = gamma_pars$scale,
    log.p = TRUE)

  return(gamma_cum_prob)
}

#' @title Computes integral of a very peaked function, modified from the
#' SADISA package
#' @description computes the logarithm of the integral of exp(logfun) from 0
#' to Inf under the following assumptions:
#' \itemize{
#'  \item{"exp(logfun)" has a single, sharply peaked maximum}
#'  \item{"exp(logfun)" is increasing to the left of the peak and
#'  decreasing to the right of the peak}
#'  \item{"exp(logfun)" can be zero or positive at zero}
#'  \item{"exp(logfun)" tends to zero at infinity}
#' }
#' @param logfun the logarithm of the function to integrate
#' @param xx the initial set of points on which to evaluate the function
#' @param xcutoff when the maximum has been found among the xx, this parameter
#' sets the width of the interval to find the maximum in
#' @param ymaxthreshold sets the deviation allowed in finding the maximum
#' among the xx
#' @return the result of the integration
#' @references Haegeman, B. & R.S. Etienne (2017). A general sampling formula
#' for community structure data. Methods in Ecology & Evolution 8: 1506-1519.
#' https://doi.org/10.1111/2041-210X.12807
#' @keywords internal
integral_peak <- function(logfun,
                          xx = seq(-20, 20, 2),
                          xcutoff = 2,
                          ymaxthreshold = 1E-12,
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
                          par_mean,
                          par_sd,
                          par_upper_bound) {
  fun <- function(x) {
    exp(logfun(x,
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
               par_mean,
               par_sd))
  }
  # determine integrand peak
  yy <- xx + logfun(exp(xx),
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
                    par_mean,
                    par_sd)
  yy[which(is.na(yy) | is.nan(yy))] <- -Inf
  yymax <- max(yy)
  if (yymax == -Inf) {
    logQ <- -Inf
    return(logQ)
  }

  iimax <- which(yy >= (yymax - ymaxthreshold))
  xlft <- xx[iimax[1]] - xcutoff
  xrgt <- xx[iimax[length(iimax)]] + xcutoff
  optfun <- function(x) {
    x + logfun(exp(x),
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
               par_mean,
               par_sd)
  }
  optres <- stats::optimize(
    f = optfun,
    interval = c(xlft, xrgt),
    maximum = TRUE,
    tol = 1e-10)
  xmax <- optres$maximum

  # compute integral
  gamma_pars <- transform_gamma_pars(
    par_mean = par_mean,
    par_sd = par_sd)
  if (gamma_pars$shape < 1) {
    lower <- min(exp(xmax), 1E-3)
    pars1f <- pars1
    pars1f[pick] <- lower / 2
    Q0 <- exp(DAISIE_loglik(
      pars1 = pars1f,
      pars2 = pars2,
      brts = brts,
      stac = stac,
      missnumspec = missnumspec,
      methode = methode,
      abstolint = abstolint,
      reltolint = reltolint,
      verbose = verbose)) *
      pracma::gammainc(lower/gamma_pars$scale,gamma_pars$shape)["reginc"]
    Q0 <- as.numeric(Q0)
  } else {
    lower <- 0
    Q0 <- 0
  }
  Q1 <- stats::integrate(f = fun,
                         lower = lower,
                         upper = min(par_upper_bound,exp(xmax)),
                         subdivisions = 1000,
                         rel.tol = 1e-10,
                         abs.tol = 1e-10)$value
  if(exp(xmax) < par_upper_bound) {
    Q2 <- stats::integrate(f = fun,
                           lower = exp(xmax),
                           upper = par_upper_bound,
                           subdivisions = 1000,
                           rel.tol = 1e-10,
                           abs.tol = 1e-10)$value
  } else {
    Q2 <- 0
  }
  logQ <- log(Q0 + Q1 + Q2)
  return(logQ)
}

#' Transforms mean and standard deviation to shape and scale gamma parameters
#'
#' @param par_mean mean of the relaxed parameter
#' @param par_sd standard deviation of the relaxed parameter
#'
#' @keywords internal
#'
#' @return list to shape and scale parameters
transform_gamma_pars <- function(par_mean, par_sd) {
  shape <- par_mean^2 / par_sd^2
  scale <- par_sd^2 / par_mean
  return(list(shape = shape,
              scale = scale))
}
