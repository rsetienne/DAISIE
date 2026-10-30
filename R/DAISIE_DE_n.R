#' @name DAISIE_DE_n
#' @title Likelihood of an endemic clade under n-sampling
#' @description Computes the log-likelihood of an endemic clade with
#'   \code{S = length(brts) - 1} sampled species and exactly
#'   \code{missnumspec = m} missing species (n-sampling) from the D-E
#'   equations. Let \code{L(z)} be the likelihood computed with initial
#'   condition \code{E(0) = z} (i.e. rho-sampling with \code{rho = 1 - z},
#'   without the factor \code{rho^S}). Then the n-sampling likelihood is
#'   \code{[z^m] L(z) / choose(S + m, S)}. Four methods are available to
#'   extract the coefficient \code{[z^m] L(z)}:
#'   \itemize{
#'     \item \code{"power_series"} (default): every state variable of the D-E
#'       equations is propagated as a power series in \code{z} truncated
#'       after degree \code{m}. The coefficients obey a closed system of
#'       ordinary differential equations (Cauchy products of the polynomial
#'       right-hand sides), so the coefficient is obtained exactly (up to the
#'       integration tolerance) from a single real-valued solve.
#'     \item \code{"cauchy"}: the Cauchy contour integral
#'       \code{(2 pi i)^-1 oint L(z) z^(-m-1) dz} over a circle through the
#'       saddle point, evaluated by adaptive quadrature with complex-valued
#'       solves of the D-E equations.
#'     \item \code{"central"} and \code{"forward"}: the \code{m}-th
#'       derivative of \code{L} at zero by numerical differentiation of
#'       \code{log L} and the complete Bell polynomials.
#'   }
#'   If the chosen method fails (an error or a non-finite result), the next
#'   method in the order power series, Cauchy integral, central
#'   differentiation, forward differentiation is tried, with a message.
#' @inheritParams default_params_doc
#' @param DAISIE_DE_function The function computing the log-likelihood of the
#'   clade for given (possibly complex) \code{missnumspec}:
#'   \code{DAISIE_DE_logpES} or \code{DAISIE_DE_logpEC}. For
#'   \code{n_method = "power_series"} it is called with
#'   \code{power_series = TRUE}.
#' @param N Not used.
#' @param n_method Method for n-sampling: \code{"power_series"} (default),
#'   \code{"cauchy"}, \code{"central"} or \code{"forward"}.
#' @return The log-likelihood.
#' @keywords internal
#' @noRd
DAISIE_DE_n <- function(DAISIE_DE_function,
                        brts,
                        missnumspec,
                        stac,
                        pars1,
                        methode,
                        reltolint,
                        abstolint,
                        N = 1000,
                        n_method = "power_series") {

  n_methods <- c("power_series", "cauchy", "central", "forward")
  n_method <- DAISIE_DE_n_method(sampling = "n", n_method = n_method)

  S <- length(brts) - 1          # given S the observed number of species

  f <- function(x) {
    rho2n <- -S + S/(1 - x)
    DAISIE_DE_function <- Vectorize(DAISIE_DE_function, "missnumspec")
    lik <- exp(DAISIE_DE_function(brts = brts,
                                  missnumspec = rho2n,
                                  stac = stac,
                                  pars1 = pars1,
                                  methode = methode,
                                  reltolint = reltolint,
                                  abstolint = abstolint))
    return(lik)
  }

  log_f <- function(x) {
    rho2n <- -S + S/(1 - x)
    DAISIE_DE_function <- Vectorize(DAISIE_DE_function, "missnumspec")
    loglik <- DAISIE_DE_function(brts = brts,
                                 missnumspec = rho2n,
                                 stac = stac,
                                 pars1 = pars1,
                                 methode = methode,
                                 reltolint = reltolint,
                                 abstolint = abstolint)
    return(loglik)
  }

  bell_polynomials_up_to_n <- function(n, g_derivs) {
    B <- numeric(n + 1)
    B[1] <- 1  # B_0
    for (m in 1:n) {
      tmp <- numeric(0)
      for (k in 1:m) {
        tmp <- c(tmp, choose(m - 1, k - 1) * B[m - k + 1] * g_derivs[k])
      }
      o <- order(abs(tmp))
      B[m + 1] <- sum(tmp[o])
    }
    return(B)  # B[1] = B_0, ..., B[n+1] = B_n
  }

  nth_derivative_from_log <- function(n, log_f, numderiv_method = 'central') {
    g_derivs <- rep(0,n)
    for(i in 1:n) {
      g_derivs[i] <- suppressWarnings(pracma::fderiv(log_f, x = 0, n = i, method = numderiv_method))
    }
    B <- bell_polynomials_up_to_n(n, g_derivs)
    return(exp(log_f(0)) * B[n + 1])
  }

  find_saddle_point_radius <- function(log_f, n, lower_r = 0.001, upper_r = 0.999) {
    # Exponent g(r) = log_f(r) - n * log(r)
    # We want to find where g'(r) = 0, which means d(log_f)/dr - n/r = 0
    saddle_equation <- function(r) {
      h <- 1e-6 # Step size for central finite difference
      # Numerical derivative of your log_f function at r
      d_log_f_dr <- (log_f(r + h) - log_f(r - h)) / (2 * h)
      # The condition for the saddle point
      return(d_log_f_dr - (n / r))
    }
    # 2. Use uniroot to solve for where the equation equals 0
    if(sign(saddle_equation(lower_r)) == sign(saddle_equation(upper_r))) {
      result <- 1
    } else
    {
      result <- stats::uniroot(saddle_equation, interval = c(lower_r, upper_r), tol = 1e-8)$root
    }
    return(result)
  }

  find_saddle_point_radii <- function(log_f, n) {
    nd <- length(n)
    # Loss function: we want the gradient of the exponent to equal zero
    saddle_loss <- function(r) {
      # Numerical derivatives of log_f at r
      d_dx <- rep(0,nd)
      diffs <- rep(0,nd)
      for(i in 1:length(n)) {
        h <- rep(0,n)
        h[i] <- 1e-6
        d_dx[i] <- (log_f(r + h) - log_f(r - h)) / (2 * h[i])
        diffs[i] <- d_dx[i] - n[i]/r[i]
      }
      return(diffs)
    }
    # Solve the system of non-linear equations
    sol <- rootSolve::multiroot(f = saddle_loss, start = rep(0.5, nd))
    return(sol$root) # Returns a vector r
  }

  integrand <- function(t, log_f, n, r = 1) {
    nd <- length(n)
    t <- matrix(t, nrow = nd)
    z <- 0 * t
    # 1. Map the real parameter t to a complex circle of radius r
    for(i in 1:nd) z[i,] <- r[i] * exp(1i * t[i,])
    # 2. Evaluate log_f(z) over the vector of z positions
    # 3. Combine the terms in the complex exponent
    complex_exponent <- log_f(z)
    for(i in 1:nd) complex_exponent <- complex_exponent - n[i] * (log(r[i]) + 1i * t[i,])
    # 4. Use the Log-Sum-Exp offset trick to prevent numeric underflow/overflow
    # Shift the exponent relative to its maximum real value before exponentiating
    offset <- max(Re(complex_exponent))
    scaled_fz_dz <- exp(complex_exponent - offset)
    # 5. Bring back the scale factor, divide by 2*pi, and extract the Real part
    result <- Re( (scaled_fz_dz * exp(offset)) / (2 * pi)^nd )
    return(result)
  }

  #loglikelihood <- log(nth_derivative_from_log(n = missnumspec, f_val = log_f) + lfactorial(S) - lfactorial(S + missnumspec)

  check <- function(x) {
    if (!is.finite(x)) stop("non-finite result")
    x
  }

  run_method <- function(method) {
    switch(method,
      power_series = check(
        DAISIE_DE_function(brts = brts,
                           missnumspec = missnumspec,
                           stac = stac,
                           pars1 = pars1,
                           methode = methode,
                           reltolint = reltolint,
                           abstolint = abstolint,
                           power_series = TRUE) -
          lchoose(S + missnumspec, S)),
      cauchy = check(
        log(stats::integrate(integrand, lower = 0, upper = 2 * pi,
                             rel.tol = max(reltolint, 1e-10),
                             abs.tol = abstolint,
                             log_f = log_f, n = missnumspec,
                             r = find_saddle_point_radius(log_f = log_f,
                                                          n = missnumspec))$value) -
          lchoose(S + missnumspec, S)),
      central = check(
        log(nth_derivative_from_log(n = missnumspec, log_f = log_f,
                                    numderiv_method = "central")) +
          lfactorial(S) - lfactorial(S + missnumspec)),
      forward = check(
        log(nth_derivative_from_log(n = missnumspec, log_f = log_f,
                                    numderiv_method = "forward")) +
          lfactorial(S) - lfactorial(S + missnumspec))
    )
  }
  failure_message <- c(
    power_series = "Power series failed (%s); switching to Cauchy integral ...",
    cauchy = "Cauchy integral failed (%s); switching to differentiation ...",
    central = "Central differentiation failed (%s); switching to forward method.")
  final_message <- c(
    power_series = paste("Power series, Cauchy integral, central and forward",
                         "differentiation all failed."),
    cauchy = "Cauchy integral, central and forward differentiation all failed.",
    central = "Central and forward differentiation both failed.",
    forward = "Forward differentiation failed.")

  # Try the chosen method and, if it fails, the next ones in the list.
  chain <- n_methods[match(n_method, n_methods):length(n_methods)]
  try_chain <- function(k) {
    tryCatch(
      run_method(chain[k]),
      error = function(e) {
        if (k == length(chain)) {
          stop(final_message[[n_method]], " Last error: ",
               conditionMessage(e), call. = FALSE)
        }
        message(sprintf(failure_message[[chain[k]]], conditionMessage(e)))
        try_chain(k + 1)
      }
    )
  }
  loglikelihood <- try_chain(1)
  # loglikelihood1 <- tryCatch({
  #   log(cubature::hcubature(f = integrand,
  #                      lowerLimit = rep(0, length(missnumspec)),
  #                      upperLimit = rep(2 * pi, length(missnumspec)),
  #                      vectorInterface = TRUE,
  #                      tol = 1e-5,
  #                      log_f = log_f,
  #                      n = missnumspec,
  #                      r = find_saddle_point_radii(log_f = log_f, n = missnumspec),
  #                      tol = 1e-8)$integral) - lchoose(S + missnumspec, S)
  # }, error = function(e) {
  #   message("Multi-D Cauchy integral failed; switching to differentiation ...")
  #   fallback <- tryCatch({
  #     log(nth_derivative_from_log(n = missnumspec, f_val = f(0))) + lfactorial(S) - lfactorial(S + missnumspec)
  #   }, error = function(e2) {
  #     stop("Both integration and differentiation failed: ", e2$message)
  #   })
  #})

  return(loglikelihood)
}


#' @title Method for n-sampling in DAISIE_DE
#' @description Returns the method used for clades with missing species:
#'   \code{NA} when \code{sampling = "rho"} (the method is then not used),
#'   and otherwise \code{n_method}, checked against the available methods,
#'   with \code{NULL} or \code{NA} meaning the default
#'   \code{"power_series"}.
#' @param sampling \code{"n"} or \code{"rho"}.
#' @param n_method \code{"power_series"}, \code{"cauchy"}, \code{"central"},
#'   \code{"forward"}, \code{NA} or \code{NULL}.
#' @return A character string, or \code{NA_character_} for rho-sampling.
#' @keywords internal
#' @noRd
DAISIE_DE_n_method <- function(sampling, n_method = NULL) {
  if (identical(sampling, "rho")) {
    return(NA_character_)
  }
  if (is.null(n_method) || (length(n_method) == 1 && is.na(n_method))) {
    return("power_series")
  }
  return(match.arg(n_method, c("power_series", "cauchy", "central", "forward")))
}
