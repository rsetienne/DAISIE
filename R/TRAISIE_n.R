TRAISIE_n <- function(
    rho_likelihood,
    sampled,
    missing,
    radius = 0.5,
    N = 16L,
    tolerance = 1e-5,
    max_evaluations = 20000L
) {


  num_observed_states <- length(sampled)

  if (!length(radius) %in% c(1L, num_observed_states)) {
    stop("Supply one radius or one per observed trait state.")
  }

  radius <- rep_len(radius, num_observed_states)

  if (any(!is.finite(radius) | radius <= 0 | radius >= 1)) {
    stop("Every radius must lie strictly between 0 and 1.")
  }

  compute_loglik_rho_sampling <- function(rho) {
    value <- rho_likelihood(rho)

    value
  }

  # No missing species: rho-sampling likelihood.
  if (all(missing == 0)) {
    return(compute_loglik_rho_sampling(rep(1, num_observed_states)))
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

    # Mirror point of each grid point: angle -> 2 * pi - angle.
    # At the mirror point z becomes its complex conjugate. All rates are real
    # numbers, so the likelihood there is the complex conjugate of the
    # likelihood at the original point. We only compute one point of each
    # pair and copy the other one, which halves the work.
    # (grid row number of a combination g = 1 + sum((g - 1) * grid_size^(k - 1)))
    place_value <- grid_size ^ (seq_len(ncol(grid)) - 1)
    mirror_row <- 1 + as.vector((grid_size - grid) %*% place_value)

    # Evaluate the Cauchy integrand at each combination of angles
    for (i in seq_len(nrow(grid))) {

      if (mirror_row[i] < i) {
        # mirror point already computed: copy its complex conjugate
        log_terms[i] <- Conj(log_terms[mirror_row[i]])
        next
      }

      z <- complex(num_observed_states)

      # Construct complex sampling variables for states with missing species
      z[states_with_missing_species] <-
        radii[states_with_missing_species] * exp(1i * angles[i, ])

      # Log generating function and Fourier coefficient-extraction factor
      log_terms[i] <-
        compute_loglik_rho_sampling(1 - z) -
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
  current <- previous

  # Check convergence by doubling the integration grid, then check the
  # result with a second radius.
  # Only done when one trait state has missing species: with two or more
  # states, these checks multiply the computation time, which makes it far
  # too slow. In that case, check the result once at the end of the
  # optimization instead.
 # if (length(states_with_missing_species) == 1) {
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
        "Radius check failed. Adjust radius, N, or ODE accuracy. ",
        "current = ", current, ", alternative = ", alternative
      )
    }

  #}



  # Return log(P_n), not P_n.
  current
}

