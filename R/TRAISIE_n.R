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
  # For missing > 0, it must accept complex sampling fractions.

  valid_counts <- function(x) {
    is.numeric(x) &&
      length(x) > 0L &&
      all(is.finite(x)) &&
      all(x >= 0 & x == floor(x))
  }

  if (!valid_counts(sampled) ||
      !valid_counts(missing) ||
      length(sampled) != length(missing)) {
    stop("sampled and missing must be equally sized vectors of counts.")
  }

  if (length(N) != 1L || !is.finite(N) ||
      N < 2 || N != floor(N)) {
    stop("N must be an integer >= 2.")
  }

  if (length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance <= 0) {
    stop("tolerance must be positive.")
  }

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

    if (is.list(value)) {
      value <- value$loglik
    }

    if (!(is.numeric(value) || is.complex(value)) ||
        length(value) != 1L ||
        is.na(value) ||
        !is.finite(Im(value)) ||
        Re(value) == Inf) {
      stop("loglik_rho must return one valid log-likelihood.")
    }

    value
  }

  # No missing species: complete-sampling likelihood.
  if (all(missing == 0)) {
    return(evaluate_loglik(rep(1, dimension)))
  }

  # Dimensions with zero missing species are evaluated at x = 0.
  active <- which(missing > 0)

  grid_size <- max(
    as.integer(N),
    2L * (max(missing) + 1L)
  )

  contour_loglik <- function(grid_size, radii) {
    number_of_points <- grid_size ^ length(active)

    if (number_of_points > max_evaluations) {
      stop(
        "Contour grid exceeds max_evaluations. ",
        "Trait-specific sampling becomes expensive with many states."
      )
    }

    grid <- as.matrix(
      expand.grid(
        rep(list(seq_len(grid_size)), length(active))
      )
    )

    angles <- 2 * pi * (grid - 0.5) / grid_size
    log_terms <- complex(nrow(grid))

    for (i in seq_len(nrow(grid))) {
      z <- complex(dimension)

      z[active] <- radii[active] *
        exp(1i * angles[i, ])

      # f(z) = P_(1-z) / prod((1-z)^sampled)
      log_terms[i] <-
        evaluate_loglik(1 - z) -
        sum(sampled * log(1 - z)) -
        1i * sum(missing[active] * angles[i, ])
    }

    # Scale before exponentiating to reduce underflow.
    shift <- max(Re(log_terms))

    if (!is.finite(shift)) {
      stop("All contour likelihoods are zero; cannot extract coefficient.")
    }

    terms <- exp(log_terms - shift)
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

    shift +
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
  alternative <- contour_loglik(grid_size, 0.8 * radius)

  if (abs(expm1(current - alternative)) > tolerance) {
    stop(
      "Radius check failed. Adjust radius, N, or ODE accuracy."
    )
  }

  # Return log(P_n), not P_n.
  current
}
