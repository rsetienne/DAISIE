TRAISIE_n <- function(
    loglik_rho,
    sampled,
    missing,
    radius = 0.5,
    N = 16L,
    tolerance = 1e-5,
    max_evaluations = 20000L
) {
  # sampled and missing:
  #   scalars -> total counts, using one shared sampling fraction
  #   vectors -> counts for each observed trait state
  #
  # loglik_rho(rho, ...) must return:
  #   a scalar log-likelihood, or a list containing $loglik.
  #



  dimension <- length(sampled)

  if (!length(radius) %in% c(1L, dimension)) {
    stop("Supply one radius or one per observed trait state.")
  }

  radius <- rep_len(radius, dimension)

  if (any(!is.finite(radius) | radius <= 0 | radius >= 1)) {
    stop("Every radius must lie strictly between 0 and 1.")
  }

  evaluate_loglik <- function(rho) {
    value <- loglik_rho(rho)

    value
  }

  # No missing species: rho-sampling likelihood.
  if (all(missing == 0)) {
    return(evaluate_loglik(rep(1, dimension)))
  }

  # Select trait states with missing species
  states_with_missing_species <- which(missing > 0)

  # Choose the number of integration points per trait state
  grid_size <- max(N, 2 * (max(missing) + 1))

  contour_loglik <- function(grid_size, radii) {
    number_of_points <- grid_size ^ length(states_with_missing_species)

    if (number_of_points > max_evaluations) {
      stop(
        "Contour grid exceeds max_evaluations. ",
        "Trait-specific sampling becomes expensive with many states."
      )
    }

    grid <- as.matrix(
      expand.grid(
        rep(list(seq_len(grid_size)), length(states_with_missing_species))
      )
    )

    angles <- 2 * pi * (grid - 0.5) / grid_size
    log_terms <- complex(nrow(grid))

    for (i in seq_len(nrow(grid))) {
      z <- complex(dimension)

      z[states_with_missing_species] <- radii[states_with_missing_species] *
        exp(1i * angles[i, ])

      # f(z) = P_(1-z) / prod((1-z)^sampled)
      log_terms[i] <-
        evaluate_loglik(1 - z) -
        sum(sampled * log(1 - z)) -
        1i * sum(missing[states_with_missing_species] * angles[i, ])
    }

    # Rescale complex exponentials to avoid overflow
    log_scale <- max(Re(log_terms))

    # Ensure the log-scale is finite before exponentiation
    if (!is.finite(log_scale))
      stop("Invalid log-scale in Cauchy integration.")

    terms <- exp(log_terms - log_scale)
    coefficient <- mean(terms)

    if (
      Re(coefficient) <= 0 ||
      abs(Im(coefficient)) > tolerance * Re(coefficient) ||
      abs(coefficient) <
      100 * .Machine$double.eps * mean(Mod(terms))
    ) {
      stop(
        "Coefficient extraction is numerically unstable. ",
        "Adjust radius or improve likelihood accuracy."
      )
    }

    log_scale +
      log(Re(coefficient)) -
      sum(missing * log(radii)) -
      sum(lchoose(sampled + missing, missing))
  }

  previous <- contour_loglik(grid_size, radius)

  # Check convergence by doubling the integration grid.
  repeat {
    grid_size <- 2L * grid_size
    current <- contour_loglik(grid_size, radius)

    if (abs(expm1(previous - current)) <= tolerance) {
      break
    }

    previous <- current
  }

  # An independent radius check helps detect numerical problems.
  alternative <- contour_loglik(grid_size, radius * exp(-1 / pmax(missing, 1)))

  if (abs(expm1(current - alternative)) > tolerance) {
    stop(
      "Radius check failed. Adjust radius, N, or ODE accuracy."
    )
  }

  # Return log(P_n), not P_n.
  current
}
