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
    return(loglik_DAISIE_par)
  }

  integrated_loglik <- log(stats::integrate(
    f = Vectorize(DAISIE_loglik_integrand,
                  vectorize.args = "DAISIE_par"),
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
    sd = sd)$value)

    return(integrated_loglik)
}
