#' Integrates the loglikelihood of a single clade across a parameter weighted by a
#' given distribution
#'
#' @inheritParams DAISIE_loglik_CS
#' @param brts Numeric vector of branching times
#' @param stac Numeric of Endemicity status
#' @param missnumspec Number of missing species
#' @param verbose Logical determining whether output is printed to the concole
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
    return(log(loglik_DAISIE_par))
  }

  integrated_loglik <- integral_peak(
    logfun = Vectorize(DAISIE_loglik_integrand,
                       vectorize.args = "DAISIE_par"),
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
    sd = sd
  )
  return(integrated_loglik)
}

#' @title Computes integral of a very peaked function, modified from the SADISA package
#' @description   # computes the logarithm of the integral of exp(logfun) from 0 to Inf under the following assumptions:
# . exp(logfun) has a single, sharply peaked maximum
# . exp(logfun) is increasing to the left of the peak and decreasing to the right of the peak
# . exp(logfun) can be zero or positive at zero
# . exp(logfun) tends to zero at infinity
#' @param logfun the logarithm of the function to integrate
#' @param xx the initial set of points on which to evaluate the function
#' @param xcutoff when the maximum has been found among the xx, this parameter sets the width of the interval to find the maximum in
#' @param ymaxthreshold sets the deviation allowed in finding the maximum among the xx
#' @param ... any arguments of the function to optimize
#' @return the result of the integration
#' @references Haegeman, B. & R.S. Etienne (2017). A general sampling formula for community structure data. Methods in Ecology & Evolution. In press.
#' @export

integral_peak <- function(logfun,
                          xx = seq(-20,20,2),
                          xcutoff = 2,
                          ymaxthreshold = 1E-12,
                          ...) {
  fun <- function(x) exp(logfun(x, ...))
  #logQ <- log(stats::integrate(f = fun, lower = 0, upper = Inf, rel.tol = 1e-10, abs.tol = 1e-10)$value)

  # 1/ determine integrand peak
  yy <- xx + logfun(exp(xx), ...)
  yy[which(is.na(yy) | is.nan(yy))] <- -Inf
  yymax <- max(yy)
  if (yymax == -Inf) {
    logQ <- -Inf
    return(logQ)
  }

  iimax <- which(yy >= (yymax - ymaxthreshold))
  xlft <- xx[iimax[1]] - xcutoff
  xrgt <- xx[iimax[length(iimax)]] + xcutoff
  optfun <- function(x) x + logfun(exp(x), ...)
  optres <- stats::optimize(f = optfun, interval = c(xlft,xrgt), maximum = TRUE, tol = 1e-10)
  xmax <- optres$maximum
  #ymax <- optres$objective

  # 3/ compute integral
  logQ <- log(stats::integrate(f = fun, lower = 0, upper = exp(xmax), rel.tol = 1e-10, abs.tol = 1e-10)$value +
          stats::integrate(f = fun, lower = exp(xmax), upper = Inf, rel.tol = 1e-10, abs.tol = 1e-10)$value)

  #intfun <- function(x) exp((x + logfun(exp(x), ...)) - ymax)
  #corrfact <- stats::integrate(f = intfun, lower = -Inf, upper = xmax, rel.tol = 1e-10, abs.tol = 1e-10)$value +
  #            stats::integrate(f = intfun, lower = xmax, upper = Inf, rel.tol = 1e-10, abs.tol = 1e-10)$value
  #logQ <- ymax + log(corrfact)

  return(logQ)
}
