DAISIE_DE_n <- function(DAISIE_DE_function,
                        brts,
                        missnumspec,
                        stac,
                        pars1,
                        methode,
                        reltolint,
                        abstolint,
                        N_cheb = 120) {

  S <- length(brts) - 1          # given S the observed number of species

  f <- function(x) {
    rho2n <- -S + S/(1 - x)
    DAISIE_DE_function <- Vectorize(DAISIE_DE_function, "missnumspec")
    lik <- exp(-S * log(1 - x) + DAISIE_DE_function(brts = brts,
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
    loglik <- -S * log(1 - x) + DAISIE_DE_function(brts = brts,
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

  bell_polynomials_up_to_n <- function(n, g_derivs) {
    B <- numeric(n + 1)
    B[1] <- 1  # B_0

    for (m in 1:n) {
      sum_val <- 0
      for (k in 1:m) {
        sum_val <- sum_val + choose(m - 1, k - 1) * B[m - k + 1] * g_derivs[k]
      }
      B[m + 1] <- sum_val
    }

    return(B)  # B[1] = B_0, ..., B[n+1] = B_n
  }

  nth_derivative_from_log <- function(n, f_val, g_derivs) {
    B <- bell_polynomials_up_to_n(n, g_derivs)
    return(f_val * B[n + 1])
  }

  lderiv <- rep(0,missnumspec)
  for(i in 1:missnumspec) {
    lderiv[i] <- suppressWarnings(pracma::fderiv(log_f, x = 0, n = i))
  }

  loglikelihood <- log(nth_derivative_from_log(n = missnumspec, f_val = f(0), g_derivs = lderiv)) + lfactorial(S) - lfactorial(S + missnumspec)
  #loglikelihood <- log(pracma::fderiv(f, x = 0, n = missnumspec)) + lfactorial(S) - lfactorial(S + missnumspec)
  #loglikelihood <- log(calculus::derivative(f, var = c(x = 0), order = missnumspec)) + lfactorial(S) - lfactorial(S + missnumspec)
  #loglikelihood <- scaled_coeff(f, n = missnumspec, S = S, N = N_cheb)

  return(loglikelihood)
}
