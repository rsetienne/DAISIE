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
                               totaltime,
                               peak) {
  # Constant
  if (island_ontogeny == 0 || is.na(island_ontogeny)) {
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
      totaltime = totaltime,
      peak = peak
    )
    return(area)
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
  totaltime <- parsvec[17]
  peak <- parsvec[18]
  kk <- parsvec[19]
  ddep <- parsvec[20]

  lx <- (length(x) - 1)/2
  lnn <- lx + 4 + 2 * kk
  nn <- -2:(lx + 2 * kk + 1)
  nn <- pmax(rep(0, lnn), nn)

  area <- island_area_vector(
    timeval = abs(t),
    area_pars = area_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    totaltime = totaltime,
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
  dx1 <- laavec[il1 + 1] * xx2[ix1] + lacvec[il4 + 1] * xx2[ix4] +
    muvec[il2 + 1] * xx2[ix3] + lacvec[il1] * nn[in1] * xx1[ix1] +
    muvec[il2] * nn[in2] * xx1[ix2] + -(muvec[il3] + lacvec[il3]) *
    nn[in3] * xx1[ix3] + -gamvec[il3] * xx1[ix3]
  dx1[1] <- dx1[1] + laavec[il3[1]] * xx3 * (kk == 1)
  dx1[2] <- dx1[2] + 2 * lacvec[il3[1]] * xx3 * (kk == 1)
  dx2 <- gamvec[il3] * xx1[ix3] +
    lacvec[il1 + 1] * nn[in1] * xx2[ix1] + muvec[il2 + 1] * nn[in2] * xx2[ix2] +
    -(muvec[il3 + 1] + lacvec[il3 + 1]) * nn[in3 + 1] * xx2[ix3] +
    -laavec[il3 + 1] * xx2[ix3]
  dx3 <- -(laavec[il3[1]] + lacvec[il3[1]] + gamvec[il3[1]] +
             muvec[il3[1]]) * xx3
  return(list(c(dx1, dx2, dx3)))
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
  totaltime <- parsvec[17]
  peak <- parsvec[18]
  kk <- parsvec[19]
  ddep <- parsvec[20]

  lx <- (length(x))/3
  lnn <- lx + 4 + 2 * kk
  nn <- -2:(lx + 2 * kk + 1)
  nn <- pmax(rep(0, lnn), nn)

  area <- island_area_vector(
    timeval = abs(t),
    area_pars = area_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    totaltime = totaltime,
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
    x = x,
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
  # cladogenesis in colonist when k = 1: Q_M,n-1 -> Q^1_n; n+k-1 species
  # present; rate twice
  # anagenesis of reimmigrant: Q^M,k_n-1 -> Q^k,n; n+k-1+1 species present
  # cladogenesis of reimmigrant: Q^M,k_n-2 -> Q^k,n; n+k-2+1 species present;
  # rate once
  # extinction of reimmigrant: Q^M,k_n -> Q^k,n; n+k+1 species present
  # cladogenesis in one of the n+k-1 species: Q^k_n-1 -> Q^k_n;
  # n+k-1 species present; rate twice for k species
  # extinction in one of the n+1 species: Q^k_n+1 -> Q^k_n;
  # n+k+1 species present
  # outflow:
  # all events with n+k species present
  dx1 <-
    (laavec[il3] * xx3[ix3] + 2 * lacvec[il1] * xx3[ix1]) * (kk == 1) +
    laavec[il1 + 1] * xx2[ix1] + lacvec[il4 + 1] * xx2[ix4] +
    muvec[il2 + 1] * xx2[ix3] + lacvec[il1] * nn[in1] * xx1[ix1] +
    muvec[il2] * nn[in2] * xx1[ix2] +
    -(muvec[il3] + lacvec[il3]) * nn[in3] * xx1[ix3] - gamvec[il3] * xx1[ix3]
  # inflow:
  # immigration when there are n+k species: Q^k,n -> Q^M,k_n; n+k species present
  # cladogenesis in n+k-1 species: Q^M,k_n-1 -> Q^M,k_n; n+k-1+1 species
  # present; rate twice for k species
  # extinction in n+1 species: Q^M,k_n+1 -> Q^M,k_n; n+k+1+1 species present
  # outflow:
  # all events with n+k+1 species present
  dx2 <- gamvec[il3] * xx1[ix3] +
    lacvec[il1 + 1] * nn[in1] * xx2[ix1] + muvec[il2 + 1] * nn[in2] * xx2[ix2] +
    -(muvec[il3 + 1] + lacvec[il3 + 1]) * nn[in3 + 1] * xx2[ix3] +
    -laavec[il3 + 1] * xx2[ix3]
  # only when k = 1
  # inflow:
  # cladogenesis in one of the n-1 species: Q_M,n-1 -> Q_M,n; n+k-1 species present; rate once
  # extinction in one of the n+1 species: Q_M,n+1 -> Q_M,n; n+k+1 species present
  # outflow:
  # all events with n+k species present
  dx3 <- lacvec[il1] * nn[in4] * xx3[ix1] + muvec[il2] * nn[in2] * xx3[ix2] +
    -(lacvec[il3] + muvec[il3]) * nn[in3] * xx3[ix3] +
    -(laavec[il3] + gamvec[il3]) * xx3[ix3]
  return(list(c(dx1, dx2, dx3)))
}


DAISIE_integrate_time <- function(initprobs,
                                  tvec,
                                  rhs_func,
                                  pars,
                                  rtol,
                                  atol,
                                  method) {
  function_as_text <- as.character(body(rhs_func)[3])
  do_fun_1 <- grepl(pattern = "lx <- \\(length\\(x\\) - 1\\)/2", x = function_as_text)
  do_fun_2 <- grepl(pattern = "lx <- \\(length\\(x\\)\\)/3", x = function_as_text)
  if (do_fun_1) {
    y <- deSolve::ode(
      initprobs,
      tvec,
      func = DAISIE_loglik_rhs_time,
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
  } else {
    stop(
      "The integrand function is written incorrectly. ",
      "Value of 'function_as_text':", function_as_text
    )
  }
  y <- y[-1,-1]
  return(y)
}
