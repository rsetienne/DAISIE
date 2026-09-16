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

  cheb_nodes01 <- function(N){
    k <- 1:N
    (1 + cos((2*k-1)*pi/(2*N)))/2
  }

  cheb_coeff <- function(f,N) {
    x <- cheb_nodes01(N)
    t <- 2*x - 1
    fx <- f(x)

    c <- numeric(N)

    for(j in 0:(N-1)){
      Tj <- cos(j*acos(t))
      c[j+1] <- (2/N)*sum(fx*Tj)
    }

    c[1] <- c[1]/2
    c
  }

  spectral_filter <- function(c) {
    N <- length(c)-1
    alpha <- 36
    p <- 10

    for(k in 0:N){
      sigma <- exp(-alpha*(k/N)^p)
      c[k + 1] <- c[k + 1]*sigma
    }

    c
  }

  cheb_derivative <- function(c) {
    N <- length(c) - 1
    d <- numeric(N)

    d[N] <- 2 * N * c[N + 1]
    d[N - 1] <- 2 * (N - 1) * c[N]

    if(N > 2){
      for(k in (N - 2) : 1){
        d[k] <- d[k + 2] + 2 * k * c[k + 1]
      }
    }

    d[1] <- d[1]/2
    d
  }

  cheb_eval <- function(c,t) {
    N <- length(c) - 1
    b1 <- 0
    b2 <- 0

    for(j in N:1){
      b0 <- 2 * t * b1 - b2 + c[j + 1]
      b2 <- b1
      b1 <- b0
    }
    return(t * b1 - b2 + c[1])
  }

  scaled_coeff <- function(f,n,S,N = 60) {
    #c <- cheb_coeff(f,N)
    c <- pracma::chebCoeff(f, 0, 1, N); c[1] <- c[1]/2
    for(i in 1:n){
      c <- spectral_filter(c)
      c <- cheb_derivative(c)
    }
    loglik <- log(cheb_eval(c,-1)) + n * log(2) + lfactorial(S) - lfactorial(S + n)
    return(loglik)
  }

  nth_deriv_richardson <- function(f, x0, n, h = 1e-5) {
    binom_coef <- choose(n, 0:n)
    deriv <- 0
    for (k in 0:n) {
      term <- (-1)^(n - k) * binom_coef[k + 1] * f(x0 + k * h)
      deriv <- deriv + term
    }
    return(deriv / h^n)
  }

  # nth_deriv_numDeriv <- function(f, x0, n) {
  #   current_f <- f
  #   for (i in 1:n) {
  #     current_f <- function(x) numDeriv::grad(current_f, x)[1]  # Compute next derivative
  #   }
  #   return(current_f(x0))
  # }

  # nth_deriv_TensorFlow <- function(f, x0, n) {
  #   tensorflow::tf$enable_eager_execution()
  #   x <- tensorflow::tf$constant(x0, dtype = tf$float32)
  #   current_y <- f(x)
  #   for (i in 1:n) {
  #     with(tensorflow::tf$GradientTape() %as% tensorflow::tape, {
  #       current_y_val <- current_y
  #     })
  #     current_y <- tensorflow::tape$gradient(current_y_val, x)
  #   }
  # }

  # nth_deriv_madness <- function(f, x0, n) {
  #   if (n == 0) {
  #     return(f(x0))
  #   }
  #   if (n == 1) {
  #     x_mad <- madness::madness(x0)
  #     y_mad <- f(x_mad)
  #     return(madness::dvdx(y_mad))
  #   }
  #   lower_order_wrapper <- function(val) {
  #     nth_derivative_madness(f, val, n - 1)
  #   }
  #   return(madness::numderiv(f = lower_order_wrapper, val = x0))
  # }

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

  nth_derivative_from_log <- function(n, f_val, g_derivs) {
    B <- bell_polynomials_up_to_n(n, g_derivs)
    return(f_val * B[n + 1])
  }

  integrand <- function(t) {
    z <- exp(1i * t)
    fz <- 1/(2*pi*1i) * f(z)/(z^(missnumspec + 1))
    dz_dt <- 1i * z
    result <- Re(fz * dz_dt)
    return(result)
  }

  integrand2 <- function(t) {
    # 1. Map the real parameter t to the complex unit circle
    z <- exp(1i * t)

    # 2. Safely evaluate log_f(z).
    # Using sapply ensures it works even if R's 'integrate' passes a vector of t values.
    ln_fz <- sapply(z, log_f)

    # 3. Combine the terms in the complex exponent: ln(f(z)) - i * n * t
    complex_exponent <- ln_fz - 1i * missnumspec * t

    # 4. Use the Log-Sum-Exp offset trick to prevent numeric underflow/overflow
    # Shift the exponent relative to its maximum real value before exponentiating
    offset <- max(Re(complex_exponent))
    scaled_fz_dz <- exp(complex_exponent - offset)

    # 5. Bring back the scale factor, divide by 2*pi, and extract the Real part
    result <- Re( (scaled_fz_dz * exp(offset)) / (2 * pi) )

    return(result)
  }

  integrate_on_unit_circle0 <- function(f, n, r = 1.0, N = 5000) {
    # 1. Discretize around a circle of radius 'r'
    theta <- seq(0, 2 * pi, length.out = N + 1)[1:N]
    z <- r * exp(1i * theta)

    # 2. Compute log(integrand)
    h <- f(z) - 1i * n * theta - n * log(r)

    # 3. Pull out the maximum real part to stabilize exponentiation
    R_max <- max(Re(h))

    # 4. Exponentiate the safe, scaled-down values
    stabilized_terms <- exp(h - R_max)

    # 5. Average the terms across the grid points
    mean_stabilized <- mean(stabilized_terms)

    # 6. Safety check for zero coefficients
    if (Mod(mean_stabilized) < 1e-12) {
      return(-Inf)
    }

    # 7. Reconstruct the final result in the log domain
    log_c_n <- R_max + log(mean_stabilized)

    return(Re(log_c_n))
  }

  integrate_on_unit_circle1 <- function(f, n, N = 8192) {
    r_opt <- 1.0

    # 1. Discretize around the optimal circle
    # We make N very large for high n to avoid phase aliasing
    theta <- seq(0, 2 * pi, length.out = N + 1)[1:N]
    z <- r_opt * exp(1i * theta)

    # 2. Evaluate your completely RAW, unmodified function
    h_raw <- f(z)

    # 3. AUTOMATIC BRANCH-CUT UNWRAPPING
    # We extract the imaginary part and force it to be continuous
    img_part <- Im(h_raw)

    # 4. Detect jumps greater than pi and smooth them out sequentially
    for (i in 2:length(img_part)) {
      diff_phase <- img_part[i] - img_part[i-1]
      if (diff_phase > pi) {
        img_part[i:length(img_part)] <- img_part[i:length(img_part)] - 2 * pi
      } else if (diff_phase < -pi) {
        img_part[i:length(img_part)] <- img_part[i:length(img_part)] + 2 * pi
      }
    }

    # 5. Re-combine the smooth imaginary part with the real part
    h_smooth <- Re(h_raw) + 1i * img_part

    # 6. Apply Cauchy's structural scaling for the n-th coefficient
    h_total <- h_smooth - 1i * n * theta - n * log(r_opt)

    # 7. Extract max real part for scaling stabilization
    R_max <- max(Re(h_total))
    stabilized_terms <- exp(h_total - R_max)
    mean_stabilized <- mean(stabilized_terms)

    if (Mod(mean_stabilized) < 1e-12) {
      return(-Inf)
    }

    # 8. Reconstruct final log-coefficient
    log_c_n <- R_max + log(mean_stabilized)

    return(Re(log_c_n))
  }

  integrate_on_unit_circle2 <- function(f, n, r = 1.0, N = 20) {

    # 1. We redefine the problem in terms of theta on a circle of radius r
    # z = r * cos(theta) + i * r * sin(theta)
    # The term exp(-i * n * theta) is split into cos(n*theta) - i * sin(n*theta)

    # We first do a quick evaluations to find a scaling factor to prevent overflow
    test_theta <- seq(0, 2 * pi, length.out = N)
    test_z <- r * exp(1i * test_theta)
    R_max <- max(Re(f(test_z))) - n * log(r)

    # 2. Define the scaled real component of your log-function f(z)
    # Note: f only takes 'z'
    integrand_complex <- function(theta) {
      z <- r * exp(1i * theta)
      # Scale down by R_max to protect floating-point capacity
      return(exp(f(z) - n * log(r) - R_max))
    }

    # 3. Split the total complex product into four real-valued integrations:
    # Integral = (Real_f + i*Imag_f) * (cos(n*theta) - i*sin(n*theta))
    # Real Part = Real_f * cos(n*theta) + Imag_f * sin(n*theta)
    # Imag Part = Imag_f * cos(n*theta) - Real_f * sin(n*theta)

    term1 <- function(t) Re(integrand_complex(t)) * cos(n * t)
    term2 <- function(t) Im(integrand_complex(t)) * sin(n * t)
    term3 <- function(t) Im(integrand_complex(t)) * cos(n * t)
    term4 <- function(t) Re(integrand_complex(t)) * sin(n * t)

    # 4. Integrate using R's adaptive internal engine
    # We divide the result by 2*pi at the same time
    I_real <- (integrate(term1, 0, 2*pi)$value + integrate(term2, 0, 2*pi, abs.tol = abstolint, rel.tol = reltolint)$value) / (2 * pi)
    I_imag <- (integrate(term3, 0, 2*pi)$value - integrate(term4, 0, 2*pi, abs.tol = abstolint, rel.tol = reltolint)$value) / (2 * pi)

    integral_scaled <- I_real + 1i * I_imag

    # 5. Safety check for zero coefficients
    if (Mod(integral_scaled) < 1e-14) {
      return(-Inf)
    }

    # 6. Reconstruct the final log-coefficient by adding the scale back
    log_c_n <- R_max + log(integral_scaled)

    return(log_magnitude = Re(log_c_n))
  }

  lderiv <- rep(0,missnumspec)
  for(i in 1:missnumspec) {
    lderiv[i] <- suppressWarnings(pracma::fderiv(log_f, x = 0, n = i))
    #lderiv[i] <- suppressWarnings(nth_deriv_richardson(f = log_f, x0 = 0, n = i)) #simply inaccurate
    #lderiv[i] <- suppressWarnings(nth_deriv_numDeriv(f = log_f, x0 = 0, n = i))   #doesn't work because it gives stack overflow
    #lderiv[i] <- suppressWarnings(nth_deriv_TensorFlow(f = log_f, x0 = 0, n = i)) #doesn't work because it requires tensorflow to be installed
    #lderiv[i] <- suppressWarnings(nth_deriv_madness(f = log_f, x0 = 0, n = i))    #doesn't work, errors
    #lderiv[i] <- suppressWarnings(calculus::derivative(f = log_f, var = c(x = 0), order = i)) #doesn't work as it gives NaN
    #lderiv[i] <- suppressWarnings(pnd::GenD(FUN = log_f, x = 0, deriv.order = i)) #gives deviations
  }

  #loglikelihood <- integrate_on_unit_circle0(f = log_f, n = missnumspec, N = 8092) - lchoose(S + missnumspec, S)
  #loglikelihood <- integrate_on_unit_circle1(f = log_f, n = missnumspec, N = 8092) - lchoose(S + missnumspec, S)
  #loglikelihood <- integrate_on_unit_circle2(f = log_f, n = missnumspec, N = 8092) - lchoose(S + missnumspec, S)
  #loglikelihood <- log(integrate(integrand, lower = 0, upper = 2 * pi, abs.tol = abstolint, rel.tol = reltolint)$value) - lchoose(S + missnumspec, S)
  loglikelihood <- log(integrate(integrand2, lower = 0, upper = 2 * pi, abs.tol = abstolint, rel.tol = reltolint)$value) - lchoose(S + missnumspec, S)
  #loglikelihood <- log(nth_derivative_from_log(n = missnumspec, f_val = f(0), g_derivs = lderiv)) + lfactorial(S) - lfactorial(S + missnumspec)
  #loglikelihood <- log(pracma::fderiv(f, x = 0, n = missnumspec)) + lfactorial(S) - lfactorial(S + missnumspec)
  #loglikelihood <- log(calculus::derivative(f, var = c(x = 0), order = missnumspec)) + lfactorial(S) - lfactorial(S + missnumspec)
  #loglikelihood <- scaled_coeff(f, n = missnumspec, S = S,N = N)

  return(loglikelihood)
}
