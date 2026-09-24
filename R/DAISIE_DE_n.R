DAISIE_DE_n <- function(DAISIE_DE_function,
                        brts,
                        missnumspec,
                        stac,
                        pars1,
                        methode,
                        reltolint,
                        abstolint,
                        N = 1000) {

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

  nth_derivative_from_log <- function(n, f_val) {
    g_derivs <- rep(0,n)
    for(i in 1:n) {
      g_derivs[i] <- suppressWarnings(pracma::fderiv(f_val, x = 0, n = i))
    }
    B <- bell_polynomials_up_to_n(n, g_derivs)
    return(f_val * B[n + 1])
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

  loglikelihood <- tryCatch({
    log(integrate(integrand,
                  lower = 0,
                  upper = 2 * pi,
                  rel.tol = reltolint,
                  abs.tol = abstolint,
                  log_f = log_f,
                  n = missnumspec,
                  r = find_saddle_point_radius(log_f = log_f, n = missnumspec))$value) - lchoose(S + missnumspec, S)
  }, error = function(e) {
    message("Cauchy integral failed; switching to differentiation ...")
    fallback <- tryCatch({
      log(nth_derivative_from_log(n = missnumspec, f_val = f(0))) + lfactorial(S) - lfactorial(S + missnumspec)
    }, error = function(e2) {
      stop("Both integration and differentiation failed: ", e2$message)
    })
  })

  # loglikelihood1 <- tryCatch({
  #   log(cubature::adaptIntegrate(integrand,
  #                      lowerLimit = rep(0, length(missnumspec)),
  #                      upperLimit = rep(2 * pi, length(missnumspec)),
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

  #loglikelihood <- log(nth_derivative_from_log(n = missnumspec, f_val = f(0)) + lfactorial(S) - lfactorial(S + missnumspec)
  return(loglikelihood)
}

