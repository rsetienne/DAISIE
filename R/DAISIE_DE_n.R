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

  loglikelihood <- scaled_coeff(f, n = missnumspec, S = S, N = N_cheb)
  return(loglikelihood)
}
