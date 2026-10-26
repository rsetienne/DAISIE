#' Missing-species likelihood with adaptive integration (hcubature)
#'
#' Same inputs and same output as TRAISIE_n, but the Cauchy integral is
#' computed with cubature::hcubature on a circle whose radius is chosen at the
#' saddle point, instead of with an evenly spaced grid of angles.
#' This is the method of the commented "loglikelihood1" code in DAISIE_DE_n,
#' written for TRAISIE. It is meant for comparing the two methods.
#'
#' @param rho_likelihood function of the sampling fractions (one per observed
#'   state) that returns the log-likelihood of the clade.
#' @param sampled number of sampled tips in each observed state.
#' @param missing number of missing species in each observed state.
#' @param tolerance relative accuracy asked from hcubature.
#' @param max_evaluations largest number of likelihood evaluations allowed.
#' @return log(P_n), like TRAISIE_n.
#' @keywords internal
TRAISIE_n_hcubature <- function(
    rho_likelihood,
    sampled,
    missing,
    tolerance = 1e-5,
    max_evaluations = 20000L
) {

  if (!requireNamespace("cubature", quietly = TRUE)) {
    stop("Package 'cubature' is needed for n_integration = 'hcubature'.")
  }

  num_observed_states <- length(sampled)

  # No missing species: rho-sampling likelihood.
  if (all(missing == 0)) {
    return(rho_likelihood(rep(1, num_observed_states)))
  }

  # Only the trait states with missing species are integrated over
  states_with_missing_species <- which(missing > 0)
  n <- missing[states_with_missing_species]
  number_of_dimensions <- length(n)

  # Log of the generating function at the complex point(s) z_part, given only
  # for the states with missing species (the other states get z = 0).
  # It is the same function as in TRAISIE_n:
  #   log g(z) = loglik(rho = 1 - z) - sum(sampled * log(1 - z))
  log_g <- function(z_part) {
    z <- complex(num_observed_states)
    z[states_with_missing_species] <- z_part
    rho_likelihood(1 - z) - sum(sampled * log(1 - z))
  }

  # ------------------------------------------------------------
  # Radius: saddle point on the positive real line.
  # The saddle point is where  log g(r) - sum(n * log(r))  is smallest.
  # (find_saddle_point_radii in DAISIE_DE_n looks for the zero of its
  #  derivative with rootSolve::multiroot. That often jumps outside the unit
  #  circle and stops. This function is convex in log(r), so we minimise it
  #  directly, with r kept between 0.01 and 0.95, which finds the same point.)
  # ------------------------------------------------------------
  log_g_real <- function(r) Re(log_g(complex(real = r)))

  saddle_objective <- function(log_r) {
    log_g_real(exp(log_r)) - sum(n * log_r)
  }

  # start from the radius used by TRAISIE_n
  start_radius <- pmin(pmax(n / (sampled[states_with_missing_species] + n), 0.1), 0.95)

  radius <- tryCatch({
    fit <- stats::optim(par = log(start_radius),
                        fn = saddle_objective,
                        method = "L-BFGS-B",
                        lower = rep(log(0.01), number_of_dimensions),
                        upper = rep(log(0.95), number_of_dimensions))
    exp(fit$par)
  }, error = function(e) start_radius)

  # ------------------------------------------------------------
  # Scale factor: on a circle, |g| is largest at the real point (all
  # coefficients of g are >= 0). Dividing by that value keeps the numbers
  # around 1, and it is the same constant for every evaluation.
  # ------------------------------------------------------------
  log_scale <- log_g_real(radius) - sum(n * log(radius))

  number_of_evaluations <- 0

  # Integrand of the Cauchy formula. hcubature gives a matrix of angles,
  # one column per point (vectorInterface = TRUE).
  integrand <- function(angles) {
    angles <- matrix(angles, nrow = number_of_dimensions)
    values <- numeric(ncol(angles))
    for (k in seq_len(ncol(angles))) {
      z_part <- radius * exp(1i * angles[, k])
      exponent <- log_g(z_part) - sum(n * (log(radius) + 1i * angles[, k]))
      values[k] <- Re(exp(exponent - log_scale))
    }
    number_of_evaluations <<- number_of_evaluations + ncol(angles)
    matrix(values / (2 * pi) ^ number_of_dimensions, nrow = 1)
  }

  # The real part of the integrand is the same at angles t and -t
  # (complex conjugate points, all rates are real numbers). So we only
  # integrate the first angle over [0, pi] and multiply by 2, like the
  # mirror points in TRAISIE_n.
  result <- cubature::hcubature(
    f = integrand,
    lowerLimit = rep(0, number_of_dimensions),
    upperLimit = c(pi, rep(2 * pi, number_of_dimensions - 1)),
    vectorInterface = TRUE,
    tol = tolerance,
    maxEval = max_evaluations
  )

  integral <- 2 * result$integral

  if (!is.finite(integral) || integral <= 0) {
    stop("hcubature gave a non-positive value for the Cauchy integral.")
  }
  if (result$returnCode != 0) {
    warning("hcubature did not reach the requested accuracy.")
  }

  log_P_n <- log_scale + log(integral) -
    sum(lchoose(sampled + missing, missing))

  # extra information for comparisons
  attr(log_P_n, "evaluations") <- number_of_evaluations
  attr(log_P_n, "radius") <- radius
  attr(log_P_n, "relative_error") <- result$error / result$integral

  log_P_n
}
