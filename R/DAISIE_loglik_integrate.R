#' Integrates the loglikelihood of a single clade across a parameter weighted by a
#' given distribution
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
  sd <- CS_version$sd
  pick <- which(c("cladogenesis",
                  "extinction",
                  "carrying_capacity",
                  "immigration",
                  "anagenesis") == CS_version$relaxed_par)
  mean <- pars1[pick]

  integrated_loglik <- integral_peak(
    logfun = Vectorize(DAISIE_loglik_integrand,
                       vectorize.args = "DAISIE_par"),
    xx = sort(c(seq(-20,20,2),
                seq(log(mean) - 1,log(mean) + 1),
                log((mean + 10 * sd)/mean))),
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

#' Gamma distribution density parameterised with mean and standard deviation
#'
#' @inheritParams default_params_doc
#'
#' @return Numeric
#' @keywords internal
rho <- function(DAISIE_par, DAISIE_dist_pars) {
  return(stats::dgamma(x = DAISIE_par,
                       shape = DAISIE_dist_pars[1]^2 / DAISIE_dist_pars[2]^2,
                       scale = DAISIE_dist_pars[2]^2 / DAISIE_dist_pars[1],
                       log = TRUE))
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
                                    mean,
                                    sd) {
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
      DAISIE_dist_pars = c(mean, sd)
    )
  return(loglik_DAISIE_par)
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
#' @keywords internal
integral_peak <- function(logfun,
                          xx = seq(-20,20,2),
                          xcutoff = 2,
                          ymaxthreshold = 1E-12,
                          ...) {
  fun <- function(x) exp(logfun(x, ...))
  #logQ <- log(stats::integrate(f = fun, lower = 0, upper = Inf, rel.tol = 1e-10, abs.tol = 1e-10)$value)

  # 1 determine integrand peak
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
  optres <- stats::optimize(f = optfun,
                            interval = c(xlft,xrgt),
                            maximum = TRUE,
                            tol = 1e-10)
  xmax <- optres$maximum
  #ymax <- optres$objective

  # 2 compute integral
  gamma_pars <- transform_gamma_pars(...)
  if(gamma_pars$shape < 1) {
    lower <- min(exp(xmax),1E-3)
    Q0 <- fun(exp(lower/2))/stats::dgamma(x = lower/2,
                                        shape = gamma_pars$shape,
                                        scale = gamma_pars$scale,
                                        log = FALSE) *
      pracma::gammainc(lower/gamma_pars$scale,gamma_pars$shape)['reginc']
  } else {
    lower <- 0
    Q0 <- 0
  }
  Q1 <- stats::integrate(f = fun,
                         lower = lower,
                         upper = exp(xmax),
                         subdivisions = 1000,
                         rel.tol = 1e-10,
                         abs.tol = 1e-10)
  Q2 <- stats::integrate(f = fun,
                         lower = exp(xmax),
                         upper = Inf,
                         subdivisions = 1000,
                         rel.tol = 1e-10,
                         abs.tol = 1e-10)
  Q1 <- Q1$value
  Q2 <- Q2$value
  logQ <- log(Q0 + Q1 + Q2)

  #intfun <- function(x) exp((x + logfun(exp(x), ...)) - ymax)
  #corrfact <- stats::integrate(f = intfun, lower = -Inf, upper = xmax, rel.tol = 1e-10, abs.tol = 1e-10)$value +
  #            stats::integrate(f = intfun, lower = xmax, upper = Inf, rel.tol = 1e-10, abs.tol = 1e-10)$value
  #logQ <- ymax + log(corrfact)

  return(logQ)
}

transform_gamma_pars <- function(mean, sd, ...) {
  shape <- mean^2 / sd^2
  scale <- sd^2 / mean
  return(list(shape = shape,
              scale = scale))
}
