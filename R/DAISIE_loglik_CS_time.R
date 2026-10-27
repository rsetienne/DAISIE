#' Computes island_area, but takes vector as argument (needed by )
#' @param timeval current time of simulation
#' @param area_pars A vector similar to list produced by create_area_pars
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the
#'   peak area is achieved}
#'   \item{[3]: total island age}
#' }
#' @param island_ontogeny a string describing the type of island ontogeny. Can be \code{NULL},
#' \code{"beta"} for a beta function describing area through time.
#' @param sea_level a numeric describing the type of sea level.
#' @family rate calculations
#' @author Pedro Neves
#' @keywords internal
#'
#' @references Valente, Luis M., Rampal S. Etienne, and Albert B. Phillimore.
#' "The effects of island ontogeny on species diversity and phylogeny."
#' Proceedings of the Royal Society of London B: Biological Sciences 281.1784 (2014): 20133227.
island_area_vector <- function(timeval,
                               area_pars,
                               island_ontogeny,
                               sea_level,
                               total_time,
                               peak) {
  # Constant: no island ontogeny and no sea level. (With sea level but without
  # island ontogeny the area does change through time, see island_area().)
  if (is.na(island_ontogeny) ||
      (island_ontogeny == 0 && (is.na(sea_level) || sea_level == 0))) {
    if (area_pars[1] != 1 || is.null(area_pars[1])) {
      warning("Constant ontogeny requires a maximum area of 1.")
    }
    return(1)
  } else { # Ontogeny
    area_pars <- create_area_pars(max_area = area_pars[1],
                                  current_area = area_pars[2],
                                  proportional_peak_t = area_pars[3],
                                  total_island_age = area_pars[4],
                                  sea_level_amplitude = area_pars[5],
                                  sea_level_frequency = area_pars[6],
                                  island_gradient_angle = area_pars[7])
    area <- island_area(
      timeval = timeval,
      area_pars = area_pars,
      island_ontogeny = island_ontogeny,
      sea_level = sea_level,
      total_time = total_time,
      peak = peak
    )
    # floor: an area that underflows to exactly 0 gives K * A = Inf * 0 = NaN
    # and A^-x = Inf (rates wrongly set to 0 / extcutoff)
    return(max(area, 1e-300))
  }
}

DAISIE_loglik_rhs_time <- function(t, x, parsvec) {
  lac0 <- parsvec[1]
  mu0 <- parsvec[2]
  K0 <- parsvec[3]
  gam0 <- parsvec[4]
  laa0 <- parsvec[5]
  d <- parsvec[6]
  x_hyperpar <- parsvec[7]
  area_pars <- parsvec[8:14]
  island_ontogeny <- parsvec[15]
  sea_level <- parsvec[16]
  total_time <- parsvec[17]
  peak <- parsvec[18]
  kk <- parsvec[19]
  ddep <- parsvec[20]

  lx <- (length(x) - 1)/2
  lnn <- lx + 4 + 2 * kk
  nn <- -2:(lx + 2 * kk + 1)
  nn <- pmax(rep(0, lnn), nn)

  area <- island_area_vector(
    timeval = total_time - abs(t),
    area_pars = area_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    total_time = total_time,
    peak = peak
  )

  lacvec <- get_clado_rate_per_capita(
    lac = lac0,
    d = d,
    A = area,
    K = K0,
    num_spec = nn
  )
  muvec <- rep(1, lnn) * get_ext_rate_per_capita(
    mu = mu0,
    x = x_hyperpar,
    A = area,
    extcutoff = 1000000
  )
  gamvec <- get_immig_rate_per_capita(
    gam = gam0,
    A = area,
    num_spec = nn,
    K = K0
  )
  laavec <- laa0 * rep(1, lnn)

  xx1 <- c(0, 0, x[1:lx], 0)
  xx2 <- c(0, 0, x[(lx + 1):(2 * lx)], 0)
  xx3 <- x[2 * lx + 1]
  nil2lx <- 3:(lx + 2)
  il1 <- nil2lx + kk - 1
  il2 <- nil2lx + kk + 1
  il3 <- nil2lx + kk
  il4 <- nil2lx + kk - 2
  in1 <- nil2lx + 2 * kk - 1
  in2 <- nil2lx + 1
  in3 <- nil2lx + kk
  ix1 <- nil2lx - 1
  ix2 <- nil2lx + 1
  ix3 <- nil2lx
  ix4 <- nil2lx - 2

  dx1 = laavec[il1 + 1] * xx2[ix1] +
    lacvec[il4 + 1] * xx2[ix4] +
    muvec[il2 + 1] * xx2[ix3] +
    lacvec[il1] * nn[in1] * xx1[ix1] +
    muvec[il2] * nn[in2] * xx1[ix2] +
    -(muvec[il3] + lacvec[il3]) * nn[in3] * xx1[ix3] +
    -gamvec[il3] * xx1[ix3]
  # The next two lines are relicts because the k = 1 case is dealt with by rhs2
  # dx1[1] = dx1[1] + laavec[il3[1]] * xx3 * (kk == 1)
  # dx1[2] = dx1[2] + 2 * lacvec[il3[1]] * xx3 * (kk == 1)

  dx2 = gamvec[il3] * xx1[ix3] +
    lacvec[il1 + 1] * nn[in1] * xx2[ix1] +
    muvec[il2 + 1] * nn[in2] * xx2[ix2] +
    -(muvec[il3 + 1] + lacvec[il3 + 1]) * nn[in3 + 1] * xx2[ix3] +
    -laavec[il3 + 1] * xx2[ix3]

  # The next line is not relevant as xx3 is always 0
  # dx3 = -(laavec[il3[1]] + lacvec[il3[1]] + gamvec[il3[1]] + muvec[il3[1]]) * xx3
  # Still need to specify dx3
  dx3 <- 0

  return(list(c(dx1, dx2, dx3)))
}

DAISIE_loglik_rhs_time1 <- function(t, x, parsvec) {
  lac0 <- parsvec[1]
  mu0 <- parsvec[2]
  K0 <- parsvec[3]
  gam0 <- parsvec[4]
  laa0 <- parsvec[5]
  d <- parsvec[6]
  x_hyperpar <- parsvec[7]
  area_pars <- parsvec[8:14]
  island_ontogeny <- parsvec[15]
  sea_level <- parsvec[16]
  total_time <- parsvec[17]
  peak <- parsvec[18]
  kk <- parsvec[19]
  ddep <- parsvec[20]

  lx <- (length(x))/4
  lnn <- lx + 4 + 2 * kk
  nn <- -2:(lx + 2 * kk + 1)
  nn <- pmax(rep(0, lnn), nn)

  area <- island_area_vector(
    timeval = total_time - abs(t),
    area_pars = area_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    total_time = total_time,
    peak = peak
  )

  lacvec <- get_clado_rate_per_capita(
    lac = lac0,
    d = d,
    A = area,
    K = K0,
    num_spec = nn
  )

  muvec <- rep(1, lnn) * get_ext_rate_per_capita(
    mu = mu0,
    x = x_hyperpar,
    A = area,
    extcutoff = 1000000
  )
  gamvec <- get_immig_rate_per_capita(
    gam = gam0,
    A = area,
    num_spec = nn,
    K = K0
  )

  laavec <- laa0 * rep(1, lnn)

  xx1 <- c(0, 0, x[1:lx], 0) # Q^0_n
  xx2 <- c(0, 0, x[(lx + 1):(2 * lx)], 0) # Q^{M,0}_n
  xx3 <- c(0, 0, x[(2 * lx + 1):(3 * lx)], 0) # Q^0_{M,n}
  xx4 <- c(0, 0, x[(3 * lx + 1):(4 * lx)], 0) # Q^{M,0}_{M,n}

  nil2lx <- 3:(lx + 2)

  il1 <- nil2lx + kk - 1
  il2 <- nil2lx + kk + 1
  il3 <- nil2lx + kk
  il4 <- nil2lx + kk - 2

  in1 <- nil2lx + 2 * kk - 1
  in2 <- nil2lx + 1
  in3 <- nil2lx + kk
  in4 <- nil2lx - 1

  ix1 <- nil2lx - 1
  ix2 <- nil2lx + 1
  ix3 <- nil2lx
  ix4 <- nil2lx - 2

  # The equations are those of DAISIE_loglik_rhs1 (time-independent version);
  # only the rate vectors differ, as these are calculated at time t above.
  # Keep the two functions in sync.
  dx1 <- lacvec[il1] * nn[in1] * xx1[ix1] +
    laavec[il1 + 1] * xx2[ix1] +
    lacvec[il4 + 1] * xx2[ix4] +
    muvec[il2] * nn[in2] * xx1[ix2] +
    muvec[il3 + 1] * xx2[ix3] +
    -(muvec[il3] + lacvec[il3]) * nn[in3] * xx1[ix3] +
    -gamvec[il3] * xx1[ix3]

  dx2 <- gamvec[il3] * xx1[ix3] +
    gamvec[il3] * xx3[ix3] +
    gamvec[il3 + 1] * xx4[ix3] +
    lacvec[il1 + 1] * nn[in1] * xx2[ix1] +
    muvec[il2 + 1] * nn[in2] * xx2[ix2] +
    -(muvec[il3 + 1] + lacvec[il3 + 1]) * nn[in3] * xx2[ix3] +
    -(laavec[il3 + 1] + lacvec[il3 + 1] + muvec[il3 + 1]) * xx2[ix3]

  dx3 <- lacvec[il1] * nn[in1] * xx3[ix1] +
    laavec[il1 + 1] * xx4[ix1] +
    lacvec[il4 + 1] * xx4[ix4] +
    muvec[il2] * nn[in2] * xx3[ix2] +
    muvec[il3 + 1] * xx4[ix3] +
    -(lacvec[il3] + muvec[il3]) * nn[in3] * xx3[ix3] +
    -gamvec[il3] * xx3[ix3]

  dx4 <- lacvec[il1 + 1] * nn[in1] * xx4[ix1] +
    muvec[il2 + 1] * nn[in2] * xx4[ix2] +
    -(lacvec[il3 + 1] + muvec[il3 + 1]) * nn[in3 + 1] * xx4[ix3] +
    -laavec[il3 + 1] * xx4[ix3] +
    -gamvec[il3 + 1] * xx4[ix3]

  return(list(c(dx1, dx2, dx3, dx4)))
}

DAISIE_loglik_rhs_time2 <- function(t, x, parsvec) {
  lac0 <- parsvec[1]
  mu0 <- parsvec[2]
  K0 <- parsvec[3]
  gam0 <- parsvec[4]
  laa0 <- parsvec[5]
  d <- parsvec[6]
  x_hyperpar <- parsvec[7]
  area_pars <- parsvec[8:14]
  island_ontogeny <- parsvec[15]
  sea_level <- parsvec[16]
  total_time <- parsvec[17]
  peak <- parsvec[18]
  kk <- parsvec[19]
  ddep <- parsvec[20]

  lx <- (length(x))/3
  lnn <- lx + 4 + 2 * kk
  nn <- -2:(lx + 2 * kk + 1)
  nn <- pmax(rep(0, lnn), nn)

  area <- island_area_vector(
    timeval = total_time - abs(t),
    area_pars = area_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    total_time = total_time,
    peak = peak
  )

  lacvec <- get_clado_rate_per_capita(
    lac = lac0,
    d = d,
    A = area,
    K = K0,
    num_spec = nn
  )

  muvec <- rep(1, lnn) * get_ext_rate_per_capita(
    mu = mu0,
    x = x_hyperpar,
    A = area,
    extcutoff = 1000000
  )
  gamvec <- get_immig_rate_per_capita(
    gam = gam0,
    A = area,
    num_spec = nn,
    K = K0
  )

  laavec <- laa0 * rep(1, lnn)

  xx1 <- c(0, 0, x[1:lx], 0)
  xx2 <- c(0, 0, x[(lx + 1):(2 * lx)], 0)
  xx3 <- c(0, 0, x[(2 * lx + 1):(3 * lx)], 0)

  nil2lx <- 3:(lx + 2)

  il1 <- nil2lx + kk - 1
  il2 <- nil2lx + kk + 1
  il3 <- nil2lx + kk
  il4 <- nil2lx + kk - 2

  in1 <- nil2lx + 2 * kk - 1
  in2 <- nil2lx + 1
  in3 <- nil2lx + kk
  in4 <- nil2lx - 1

  ix1 <- nil2lx - 1
  ix2 <- nil2lx + 1
  ix3 <- nil2lx
  ix4 <- nil2lx - 2

  # inflow:
  # anagenesis in colonist when k = 1: Q_M,n -> Q^1_n; n+k species present
  # cladogenesis in colonist when k = 1: Q_M,n-1 -> Q^1_n;
  # n+k-1 species present; rate twice
  # anagenesis of reimmigrant: Q^M,k_n-1 -> Q^k,n; n+k-1+1 species present
  # cladogenesis of reimmigrant: Q^M,k_n-2 -> Q^k,n;
  # n+k-2+1 species present; rate once
  # extinction of reimmigrant: Q^M,k_n -> Q^k,n; n+k+1 species present
  # cladogenesis in one of the n+k-1 species: Q^k_n-1 -> Q^k_n;
  # n+k-1 species present; rate twice for k species
  # extinction in one of the n+1 species: Q^k_n+1 -> Q^k_n; n+k+1 species
  # present
  # outflow:
  # all events with n+k species present
  dx1 = (laavec[il3] * xx3[ix3] + 2 * lacvec[il1] * xx3[ix1]) * (kk == 1) +
    laavec[il1 + 1] * xx2[ix1] +
    lacvec[il4 + 1] * xx2[ix4] +
    muvec[il2 + 1] * xx2[ix3] +
    lacvec[il1] * nn[in1] * xx1[ix1] +
    muvec[il2] * nn[in2] * xx1[ix2] +
    -(muvec[il3] + lacvec[il3]) * nn[in3] * xx1[ix3] +
    -gamvec[il3] * xx1[ix3]

  # inflow:
  # immigration when there are n+k species: Q^k,n -> Q^M,k_n;
  # n+k species present
  # cladogenesis in n+k-1 species: Q^M,k_n-1 -> Q^M,k_n;
  # n+k-1+1 species present; rate twice for k species
  # extinction in n+1 species: Q^M,k_n+1 -> Q^M,k_n; n+k+1+1 species present
  # outflow:
  # all events with n+k+1 species present
  dx2 <- gamvec[il3] * xx1[ix3] +
    lacvec[il1 + 1] * nn[in1] * xx2[ix1] +
    muvec[il2 + 1] * nn[in2] * xx2[ix2] +
    -(muvec[il3 + 1] + lacvec[il3 + 1]) * nn[in3 + 1] * xx2[ix3] +
    -laavec[il3 + 1] * xx2[ix3]

  # only when k = 1
  # inflow:
  # cladogenesis in one of the n-1 species: Q_M,n-1 -> Q_M,n;
  # n+k-1 species present; rate once
  # extinction in one of the n+1 species: Q_M,n+1 -> Q_M,n;
  # n+k+1 species present
  # outflow:
  # all events with n+k species present
  dx3 <- lacvec[il1] * nn[in4] * xx3[ix1] + muvec[il2] * nn[in2] * xx3[ix2] +
    -(lacvec[il3] + muvec[il3]) * nn[in3] * xx3[ix3] +
    -(laavec[il3] + gamvec[il3]) * xx3[ix3]

  return(list(c(dx1, dx2, dx3)))
}

DAISIE_loglik_rhs_time3 <- function(t, x, parsvec) {
  lac0 <- parsvec[1]
  mu0 <- parsvec[2]
  K0 <- parsvec[3]
  gam0 <- parsvec[4]
  laa0 <- parsvec[5]
  d <- parsvec[6]
  x_hyperpar <- parsvec[7]
  area_pars <- parsvec[8:14]
  island_ontogeny <- parsvec[15]
  sea_level <- parsvec[16]
  total_time <- parsvec[17]
  peak <- parsvec[18]
  # parsvec[19] (kk) and parsvec[20] (ddep) are not used: the rhs3 system
  # always has k = 0 and ddep is not used in the time-dependent model

  # state: x1 (lx), x2 (lx x lx), x3 (lx x lx) with lx1 = lx2 = lx
  lx <- round(-1 / 4 + 1 / 4 * sqrt(1 + 8 * length(x)))

  area <- island_area_vector(
    timeval = total_time - abs(t),
    area_pars = area_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    total_time = total_time,
    peak = peak
  )

  # per-capita rates as a function of the total number of species m = n1 + n2
  # (m = 0, ..., 2 * lx + 1), see DAISIE_loglik_rhs_time1()
  mm <- 0:(2 * lx + 1)
  lacv <- get_clado_rate_per_capita(
    lac = lac0, d = d, A = area, K = K0, num_spec = mm
  )
  gamv <- get_immig_rate_per_capita(
    gam = gam0, A = area, num_spec = mm, K = K0
  )
  mu <- get_ext_rate_per_capita(
    mu = mu0, x = x_hyperpar, A = area, extcutoff = 1000000
  )
  laa <- laa0

  i0 <- 0:(lx - 1)
  nmat <- outer(i0, i0, "+")
  at <- function(v, shift) {
    matrix(v[pmax(nmat + shift, 0) + 1], lx, lx)
  }
  L0 <- at(lacv, 0)
  Lm <- at(lacv, -1)
  Lp <- at(lacv, 1)
  G0 <- at(gamv, 0)
  im1 <- matrix(pmax(0, i0 - 1), lx, lx)               # max(0, n1 - 1)
  jm1 <- matrix(pmax(0, i0 - 1), lx, lx, byrow = TRUE) # max(0, n2 - 1)
  ip1 <- matrix(i0 + 1, lx, lx)                        # n1 + 1
  jp1 <- matrix(i0 + 1, lx, lx, byrow = TRUE)          # n2 + 1
  ones <- matrix(1, lx, lx)

  # same coefficients as DAISIE_loglik_rhs_precomp2(), but with the
  # time-dependent per-capita rates (which include the diversity dependence)
  cp <- list(
    lx1 = lx,
    lx2 = lx,
    a1 = lacv[i0 + 1] * pmax(0, i0 - 1),
    a2 = mu * (i0 + 1),
    a3 = (lacv[i0 + 2] + mu) * i0 + gamv[i0 + 2] + lacv[i0 + 2] + laa + mu,
    b1 = laa * ones,
    b2 = Lm,
    b3 = mu * ones,
    b4 = mu * rep(1, lx),
    b5 = laa * rep(1, lx),
    b6 = lacv[pmax(0, i0 - 1) + 1],
    b7 = Lm * im1,
    b8 = Lm * jm1,
    b9 = mu * ip1,
    b10 = mu * jp1,
    b11 = (L0 + mu) * nmat + G0,
    c1 = G0,
    c2 = gamv[i0 + 2],
    c3 = L0 * im1,
    c4 = L0 * jm1,
    c5 = mu * ip1,
    c6 = mu * jp1,
    c7 = (Lp + mu) * nmat + laa + Lp + mu
  )
  # the equations are those of DAISIE_loglik_rhs3()
  return(DAISIE_loglik_rhs3(t, x, cp))
}

DAISIE_integrate_time <- function(initprobs,
                                  tvec,
                                  rhs_func,
                                  pars,
                                  rtol,
                                  atol,
                                  method) {
  function_as_text <- as.character(body(rhs_func)[2])
  if (startsWith(method, "odeint")) {
    # C++ (boost::odeint) analog of the deSolve calls below
    do_fun <- grepl(pattern = "rhs <- 0", x = function_as_text)
    do_fun_1 <- grepl(pattern = "rhs <- 1", x = function_as_text)
    do_fun_2 <- grepl(pattern = "rhs <- 2", x = function_as_text)
    do_fun_3 <- grepl(pattern = "rhs <- 3", x = function_as_text)
    runmod <- if (do_fun) {
      "daisie_runmod"
    } else if (do_fun_1) {
      "daisie_runmod1"
    } else if (do_fun_2) {
      "daisie_runmod2"
    } else if (do_fun_3) {
      "daisie_runmod3"
    } else {
      stop(
        "The integrand function is written incorrectly. ",
        "Value of 'function_as_text':", function_as_text
      )
    }
    return(DAISIE_ode_cs_time(
      initprobs = initprobs,
      tvec = tvec,
      pars = pars,
      atol = atol,
      rtol = rtol,
      methode = method,
      runmod = runmod
    ))
  }
  do_fun <- grepl(pattern = "rhs <- 0",x = function_as_text)
  do_fun_1 <- grepl(pattern = "rhs <- 1",x = function_as_text)
  do_fun_2 <- grepl(pattern = "rhs <- 2",x = function_as_text)
  do_fun_3 <- grepl(pattern = "rhs <- 3",x = function_as_text)
  if (do_fun) {
    y <- deSolve::ode(
      initprobs,
      tvec,
      func = DAISIE_loglik_rhs_time,
      pars,
      atol = atol,
      rtol = rtol,
      method = method
    )
  } else if (do_fun_1) {
    y <- deSolve::ode(
      initprobs,
      tvec,
      func = DAISIE_loglik_rhs_time1,
      pars,
      atol = atol,
      rtol = rtol,
      method = method
    )
  } else if (do_fun_2) {
    y <- deSolve::ode(
      initprobs,
      tvec,
      func = DAISIE_loglik_rhs_time2,
      pars,
      atol = atol,
      rtol = rtol,
      method = method
    )
  } else if (do_fun_3) {
    y <- deSolve::ode(
      initprobs,
      tvec,
      func = DAISIE_loglik_rhs_time3,
      pars,
      atol = atol,
      rtol = rtol,
      method = method
    )
  } else {
    stop(
      "The integrand function is written incorrectly. ",
      "Value of 'function_as_text':", function_as_text
    )
  }
  y <- y[-1,-1]
  return(y)
}

#' Integrate the time-dependent CS likelihood ODE with boost::odeint
#'
#' C++ analog of \code{DAISIE_ode_cs()} for the time-dependent model
#' (island ontogeny and/or sea level). The right-hand sides are C++ versions of
#' \code{DAISIE_loglik_rhs_time()}, \code{DAISIE_loglik_rhs_time1()} and
#' \code{DAISIE_loglik_rhs_time2()} and \code{DAISIE_loglik_rhs_time3()}.
#'
#' @param initprobs initial state vector
#' @param tvec two time points: start and end of the integration
#' @param pars the 18 time-dependent model parameters (see
#' \code{DAISIE_loglik_rhs_time()}), followed by \code{kk} and \code{ddep}
#' @param atol absolute tolerance
#' @param rtol relative tolerance
#' @param methode an odeint stepper, e.g. \code{"odeint::runge_kutta_dopri5"}
#' @param runmod one of \code{"daisie_runmod"}, \code{"daisie_runmod1"},
#' \code{"daisie_runmod2"}, \code{"daisie_runmod3"}
#' @return the state vector at the end of \code{tvec}
#' @keywords internal
DAISIE_ode_cs_time <- function(initprobs,
                               tvec,
                               pars,
                               atol,
                               rtol,
                               methode,
                               runmod = "daisie_runmod") {
  N <- length(initprobs)
  lx <- switch(
    runmod,
    daisie_runmod = (N - 1) / 2,
    daisie_runmod1 = N / 4,
    daisie_runmod2 = N / 3,
    daisie_runmod3 = -1 / 4 + 1 / 4 * sqrt(1 + 8 * N),
    stop("DAISIE_ode_cs_time: unknown runmod '", runmod, "'")
  )
  kk <- pars[19]
  if (runmod == "daisie_runmod3") {
    # k = 0 in this system; the C++ rate vectors must cover 2 * lx + 3 species
    kk <- lx
  }
  .Call(
    "daisie_odeint_cs_time", runmod, initprobs, tvec, lx, kk, pars[1:18],
    methode, atol, rtol
  )
}
