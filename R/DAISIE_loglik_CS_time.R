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
      total_time = total_time,
      peak = peak
    )
    return(area)
  }
}

DAISIE_loglik_rhs_time <- function(t, x, parsvec) {
  rhs <- 0
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

  area <- DAISIE:::island_area_vector(
    timeval = abs(t),
    area_pars = area_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    total_time = total_time,
    peak = peak
  )

  lacvec <- DAISIE:::get_clado_rate_per_capita(
    lac = lac0,
    d = d,
    A = area,
    K = K0,
    num_spec = nn
  )
  muvec <- rep(1, lnn) * DAISIE:::get_ext_rate_per_capita(
    mu = mu0,
    x = x_hyperpar,
    A = area,
    extcutoff = 1000000
  )
  gamvec <- DAISIE:::get_immig_rate_per_capita(
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
  # The next two lines are relicts because the k = 1 case is dealth with by rhs2
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
  rhs <- 1

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

  lx <- length(x) / 4
  lnn <- lx + 4 + 2 * kk

  nn <- -2:(lx + 2 * kk + 1)
  nn <- pmax(rep(0, lnn), nn)



  area <- island_area_vector(
    timeval = abs(t),
    area_pars = area_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    total_time = total_time,
    peak = peak
  )

  lacvec <- DAISIE:::get_clado_rate_per_capita(
    lac = lac0,
    d = d,
    A = area,
    K = K0,
    num_spec = nn
  )

  muvec <- rep(1, lnn) * DAISIE:::get_ext_rate_per_capita(
    mu = mu0,
    x = x_hyperpar,
    A = area,
    extcutoff = 1000000
  )
  gamvec <- DAISIE:::get_immig_rate_per_capita(
    gam = gam0,
    A = area,
    num_spec = nn,
    K = K0
  )

  laavec <- laa0 * rep(1, lnn)

  xx1 <- c(0,0,x[1:lx],0) #Q^0_n
  xx2 <- c(0,0,x[(lx + 1):(2 * lx)],0) #Q^{M,0}_n
  xx3 <- c(0,0,x[(2 * lx + 1):(3 * lx)],0) #Q^0_{M,n}
  xx4 <- c(0,0,x[(3 * lx + 1):(4 * lx)],0) #Q^{M,0}_{M,n}

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

  return(list(c(dx1,dx2,dx3,dx4)))
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
    timeval = abs(t),
    area_pars = area_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    total_time = total_time,
    peak = peak
  )

  lacvec <- DAISIE:::get_clado_rate_per_capita(
    lac = lac0,
    d = d,
    A = area,
    K = K0,
    num_spec = nn
  )

  muvec <- rep(1, lnn) * DAISIE:::get_ext_rate_per_capita(
    mu = mu0,
    x = x_hyperpar,
    A = area,
    extcutoff = 1000000
  )
  gamvec <- DAISIE:::get_immig_rate_per_capita(
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

nndivdep_CS <- function(lx1, lx2, K, k) {
  nn1 <- c(0, 0, 0:(lx1 + 1))
  nn2 <- c(0, 0, 0:(lx2 + 1))
  lnn1 <- length(nn1)
  lnn2 <- length(nn2)
  nil2lx1 <- 3:(lx1 + 2)
  nil2lx2 <- 3:(lx2 + 2)
  nn <- rowSums(expand.grid(n1 = nn1, n2 = nn2))
  dim(nn) <- c(lnn1, lnn2)
  nils <- array(0, dim = c(lx1 + 4, lx2 + 4))
  divdepfac2D <- pmax(nils, 1 - (nn + k) / K)[nil2lx1, nil2lx2]
  divdepfac2Dmin1 <- pmax(nils, 1 - (nn + k - 1) / K)[nil2lx1, nil2lx2]
  divdepfac2Dplus1 <- pmax(nils, 1 - (nn + k + 1) / K)[nil2lx1, nil2lx2]
  res <- list(lx1 = lx1,
              lx2 = lx2,
              nn = nn,
              divdepfac2D = divdepfac2D,
              divdepfac2Dmin1 = divdepfac2Dmin1,
              divdepfac2Dplus1 = divdepfac2Dplus1)
  return(res)
}


DAISIE_loglik_rhs_precomp2_time <- function(t, parslist) {

  # Baseline DAISIE parameters
  lac0 <- parslist$pars[1]
  mu0  <- parslist$pars[2]
  K0   <- parslist$pars[3]
  gam0 <- parslist$pars[4]
  laa0 <- parslist$pars[5]

  # Time-dependent area parameters
  d <- parslist$pars[6]
  x_hyperpar <- parslist$pars[7]

  area_pars <- parslist$pars[8:14]
  island_ontogeny <- parslist$pars[15]
  sea_level <- parslist$pars[16]
  total_time <- parslist$pars[17]
  peak <- parslist$pars[18]

  k <- parslist$k
  ddep <- parslist$ddep

  lx1 <- parslist$nndd$lx1
  lx2 <- parslist$nndd$lx2
   nn <- parslist$nndd$nn

  nil2lx1 <- 3:(lx1 + 2)
  nil2lx2 <- 3:(lx2 + 2)
  nil1 <- rep(1, lx1)
  nil2 <- rep(1, lx2)

  # Island area at current integration time
  area <- island_area_vector(
    timeval = abs(t),
    area_pars = area_pars,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level,
    total_time = total_time,
    peak = peak
  )

  # Helper: DAISIE rate functions may return a vector,
  # but rhs3 needs matrix indexing.
  as_nn_matrix <- function(v) {
    dim(v) <- dim(nn)
    return(v)
  }

  # Cladogenesis rate at n + k species
  lac <- as_nn_matrix(
     get_clado_rate_per_capita(
      lac = lac0,
      d = d,
      A = area,
      K = K0,
      num_spec = nn + k
    )
  )

  # Cladogenesis rate at n + k - 1 species
  lac_min1 <- as_nn_matrix(
     get_clado_rate_per_capita(
      lac = lac0,
      d = d,
      A = area,
      K = K0,
      num_spec = nn + k - 1
    )
  )

  # Cladogenesis rate at n + k + 1 species
  lac_plus1 <- as_nn_matrix(
     get_clado_rate_per_capita(
      lac = lac0,
      d = d,
      A = area,
      K = K0,
      num_spec = nn + k + 1
    )
  )

  # Extinction rate at current area
  mu <- as.numeric(
     get_ext_rate_per_capita(
      mu = mu0,
      x = x_hyperpar,
      A = area,
      extcutoff = 1000000
    )
  )[1]

  # Immigration rate at n + k species
  gam <- as_nn_matrix(
     get_immig_rate_per_capita(
      gam = gam0,
      A = area,
      K = K0,
      num_spec = nn + k
    )
  )

  # Immigration rate at n + k + 1 species
  gam_plus1 <- as_nn_matrix(
     get_immig_rate_per_capita(
      gam = gam0,
      A = area,
      K = K0,
      num_spec = nn + k + 1
    )
  )


  laa <- laa0

  cp <- list(
    lx1 = lx1,
    lx2 = lx2,

    a1 = lac[nil2lx1, nil2lx2][, 1] *
      nn[nil2lx1 - 1],

    a2 = mu * nn[nil2lx1 + 1],

    a3 = (lac_plus1[nil2lx1, nil2lx2][, 1] + mu) *
      nn[nil2lx1] +
      gam_plus1[nil2lx1, nil2lx2][, 1] +
      lac_plus1[nil2lx1, nil2lx2][, 1] +
      laa + mu,

    b1 = laa * matrix(1, lx1, lx2),

    b2 = lac_min1[nil2lx1, nil2lx2],

    b3 = mu * matrix(1, lx1, lx2),

    b4 = mu * rep(1, lx1),

    b5 = laa * rep(1, lx1),

    b6 = lac_min1[nil2lx1, nil2lx2][, 1],

    b7 = lac_min1[nil2lx1, nil2lx2] *
      nn[nil2lx1 - 1, nil2],

    b8 = lac_min1[nil2lx1, nil2lx2] *
      nn[nil1, nil2lx2 - 1],

    b9 = mu * nn[nil2lx1 + 1, nil2],

    b10 = mu * nn[nil1, nil2lx2 + 1],

    b11 = (lac[nil2lx1, nil2lx2] + mu) *
      nn[nil2lx1, nil2lx2] +
      gam[nil2lx1, nil2lx2],

    c1 = gam[nil2lx1, nil2lx2],

    c2 = gam_plus1[nil2lx1, nil2lx2][, 1],

    c3 = lac[nil2lx1, nil2lx2] *
      nn[nil2lx1 - 1, nil2],

    c4 = lac[nil2lx1, nil2lx2] *
      nn[nil1, nil2lx2 - 1],

    c5 = mu * nn[nil2lx1 + 1, nil2],

    c6 = mu * nn[nil1, nil2lx2 + 1],

    c7 = (lac_plus1[nil2lx1, nil2lx2] + mu) *
      nn[nil2lx1, nil2lx2] +
      laa +
      lac_plus1[nil2lx1, nil2lx2] +
      mu
  )

  return(cp)
}

DAISIE_loglik_rhs_time3 <- function(t,x,cp)
{
  rhs <- 3
  lx1 <- cp$lx1
  lx2 <- cp$lx2
  x1 <- x[1:lx1]
  x2 <- x[(lx1 + 1):(lx1 + lx1 * lx2)]
  x3 <- x[(lx1 + lx1 * lx2 + 1):(lx1 + 2 * lx1 * lx2)]
  dim(x2) <- c(lx1 ,lx2)
  dim(x3) <- c(lx1, lx2)
  xx1 <- rep(0,lx1 + 3)
  xx2 <- array(0,dim = c(lx1 + 3, lx2 + 3))
  xx3 <- array(0,dim = c(lx1 + 3, lx2 + 3))
  nil2lx1 <- 3:(lx1 + 2)
  nil2lx2 <- 3:(lx2 + 2)
  xx1[nil2lx1] <- x1
  xx2[nil2lx1, nil2lx2] <- x2
  xx3[nil2lx1, nil2lx2] <- x3
  dx1 <-
    cp$a1 * xx1[nil2lx1 - 1] +
    cp$a2 * xx1[nil2lx1 + 1] -
    cp$a3 * xx1[nil2lx1]
  dx2 <-
    cp$b1 * xx3[nil2lx1, nil2lx2 - 1] +
    cp$b2 * xx3[nil2lx1, nil2lx2 - 2] +
    cp$b3 * xx3[nil2lx1, nil2lx2] +
    cp$b7 * xx2[nil2lx1 - 1, nil2lx2]  +
    cp$b8 * xx2[nil2lx1, nil2lx2 - 1] +
    cp$b9 * xx2[nil2lx1 + 1, nil2lx2] +
    cp$b10 * xx2[nil2lx1, nil2lx2 + 1] -
    cp$b11 * xx2[nil2lx1, nil2lx2]
  dx2[,1] <- dx2[,1] +
    cp$b4 * xx1[nil2lx1] +
    cp$b5 * xx1[nil2lx1 - 1] +
    cp$b6 * xx1[nil2lx1 - 2]
  dim(dx2) <- c(lx1 * lx2, 1)
  dx3 <-
    cp$c1 * xx2[nil2lx1, nil2lx2] +
    cp$c3 * xx3[nil2lx1 - 1, nil2lx2] +
    cp$c4 * xx3[nil2lx1, nil2lx2 - 1]  +
    cp$c5 * xx3[nil2lx1 + 1, nil2lx2] +
    cp$c6 * xx3[nil2lx1, nil2lx2 + 1] -
    cp$c7 * xx3[nil2lx1, nil2lx2]
  dx3[,1] <- dx3[,1] +
    cp$c2 * xx1[nil2lx1]
  dim(dx3) <- c(lx1 * lx2, 1)
  return(list(c(dx1,dx2,dx3)))
}


DAISIE_integrate_time <- function(initprobs,
                                  tvec,
                                  rhs_func,
                                  pars,
                                  rtol,
                                  atol,
                                  method) {

  if (startsWith(method, "odeint::")) {
    stop(
      "odeint methods are not supported for the R time-dependent RHS. ",
      "Use methode = 'lsode' or methode = 'lsoda'."
    )
  }

  if (length(tvec) != 2 || any(is.na(tvec))) {
    stop(
      "tvec must contain exactly two non-NA time points. Current tvec = ",
      paste(tvec, collapse = ", ")
    )
  }

  if (abs(tvec[2] - tvec[1]) < .Machine$double.eps) {
    return(initprobs)
  }

  if (identical(rhs_func, DAISIE_loglik_rhs_time3)) {

    rhs_func_to_use <- function(t, x, parslist) {

      cp <- DAISIE_loglik_rhs_precomp2_time(
        t = t,
        parslist = parslist
      )

      DAISIE_loglik_rhs_time3(
        t = t,
        x = x,
        cp = cp
      )
    }

  } else {

    rhs_func_to_use <- rhs_func
  }

  y <- deSolve::ode(
    y = initprobs,
    times = tvec,
    func = rhs_func_to_use,
    parms = pars,
    atol = atol,
    rtol = rtol,
    method = method
  )

  return(as.numeric(y[nrow(y), -1]))
}



DAISIE_loglik_CS_M1_time <- DAISIE_loglik_time <- function(pars1,
                                                           pars2,
                                                           brts,
                                                           stac,
                                                           missnumspec,
                                                           methode = "odeint::runge_kutta_cash_karp54",
                                                           abstolint = 1E-16,
                                                           reltolint = 1E-10,
                                                           verbose = 0,
                                                           CS_version = list(model = 1, function_to_optimize = 'DAISIE', sampling = 'n')) {
  # stac = status of the clade formed by the immigrant
  #  . stac == 1 : immigrant is present but has not formed an extant clade
  #  . stac == 2 : immigrant is not present but has formed an extant clade
  #  . stac == 3 : immigrant is present and has formed an extant clade
  #  . stac == 4 : immigrant is present but has not formed an extant clade,
  #  and it is known when it immigrated.
  #  . stac == 5 : immigrant is not present and has not formed an extant clade,
  #  but only an endemic species
  #  . stac == 6 : like 2, but with max colonization time
  #  . stac == 7 : like 3, but with max colonization time
  #  . stac == 8 : like 1, but with min colonization time
  #  . stac == 9 : like 5, but with min colonization time
  # warn if laa becomes Inf
  if (any(is.infinite(pars1[c(1,2,4,5)])) ) {
    if (verbose) {
      message('One of the parameters is infinite.')
    }
  }
  if(is.na(pars2[4]))
  {
    pars2[4] = 0
  }
  ddep <- pars2[2]
  K <- pars1[3]
  if (!is.na(pars2[5])) {
    K <- K * pars1[8]
  }
  if(length(pars1) == 6) {
    probability_of_init_presence <- pars1[6]
    pars1 <- pars1[-6]
  } else {
    probability_of_init_presence <- 0
  }
  brts <- -sort(abs(as.numeric(brts)),decreasing = TRUE)
  if(length(brts) == 1 & sum(brts == 0) == 1)
  {
    stop('The branching times contain only a 0. This means the island emerged at the present which is not allowed.');
    loglik <- -Inf
    return(loglik)
  }
  if (sum(brts == 0) == 0) {
    brts[length(brts) + 1] <- 0
  }
  # for stac = 0, brts will contain origin of island and 0; length = 2;
  # no. species should be 0
  # for stac = 1, brts will contain origin of island, maximum colonization time
  # (usually island age) and 0; length = 3; no. species should be 1
  # for stac = 2, brts will contain origin of island, colonization event,
  # branching times, 0; no. species should be no. branching times + 1
  # for stac = 3, brts will contain origin of island, colonization event,
  # branching times, 0; no. species should be no. branching times + 2
  # for stac = 4, brts will contain origin of island, colonization event and 0;
  # length = 3; no. species should be 1
  # for stac = 5, brts will contain origin of island, maximum colonization time
  # (usually island age), and 0; length = 2; number of species should be 1 (+ missing species)
  # for stac = 6, brts will contain origin of island, maximum colonization time
  # (usually island age), branching times and 0;
  # number of species should be no. branching times + 1
  # for stac = 7, brts will contain origin of island, maximum colonization time
  #  usually island age), branching times and 0;
  #  number of species should be no. branching times + 2
  # for stac = 8, brts will contain origin of island, maximum colonization time
  #  usually island age), minimum colonization time and 0; length = 4;
  #  number of species should be 1
  # for stac = 9, brts will contain origin of island, maximum colonization time
  #  usually island age), minimum colonization time and 0; length = 4;
  #  number of species should be 1
  S <- 0 * (stac == 0) + (stac == 1 || stac == 4 || stac == 5 || stac == 8 || stac == 9) +
    (length(brts) - 2) * (stac == 2) + (length(brts) - 1) * (stac == 3) +
    (length(brts) - 2) * (stac == 6) + (length(brts) - 1) * (stac == 7)
  S2 <- S - (stac == 1) - (stac == 3) - (stac == 4) - (stac == 7)
  loglik <- -lgamma(S2 + missnumspec + 1) +
    lgamma(S2 + 1) + lgamma(missnumspec + 1)
  if (min(pars1) < 0) {
    warning("One or more parameters are negative.")
    loglik <- -Inf
    return(loglik)
  }
  if ((ddep == 1 | ddep == 11) & ceiling(K) < (S + missnumspec)) {
    warning('The proposed value of K is incompatible with the number of species
          in the clade. Likelihood for this parameter set
          will be set to -Inf. \n')
    loglik <- -Inf
    return(loglik)
  }
  lac <- pars1[1]
  mu <- pars1[2]
  if(lac > 10^3 & missnumspec == 0 & length(pars1) == 5) {
    if(verbose) message('High lambda detected; approximation used.')
    loglik <- DAISIE_loglik_high_lambda(pars1, -brts, stac)
  } else {
    if (ddep == 1 | ddep == 11) {
      lx <- min(
        1 + max(missnumspec, ceiling(K)),
        DDD::roundn(pars2[1]) + missnumspec
      )
    } else {
      lx <- DDD::roundn(pars2[1]) + missnumspec
    }
    if(loglik > -Inf)
    {
      # in all cases we integrate from the origin of the island to the colonization event
      # (stac 2, 3, 4), the first branching point (stac = 6, 7), to the maximum colonization
      # time (stac = 1, 5, 8, 9) or to the present (stac = 0)
      probs <- rep(0,2 * lx + 1)
      probs[1] <- 1 - probability_of_init_presence #Q^k_n
      probs[lx + 1] <- probability_of_init_presence #Q^{M,k}_n
      k1 <- 0
      probs <- DAISIE_integrate(probs,brts[1:2],DAISIE_loglik_rhs_time,c(pars1,k1,ddep),rtol = reltolint,atol = abstolint,method = methode)
      cp <- checkprobs2(lv = 2 * lx, loglik, probs, verbose); loglik = cp[[1]]; probs = cp[[2]]
      if(stac == 0)
      {
        # for stac = 0, the integration was from the origin of the island until
        # the present so we can immediately evaluate the probability of no clade
        # being present and no immigrant species.
        loglik <- loglik + log(probs[1])
      } else
      {
        if (stac %in% c(1, 5:9) )
        {
          # for stac = 1, we now integrate from the maximum colonization time
          # (usually the island age + tiny time unit) until the present, where
          # we set all probabilities where the immigrant is already present to 0
          # and we evaluate the probability of the immigrant species being
          # present, but there can be missing species.
          # for stac = 5, we do exactly the same, but we evaluate the
          # probability of an endemic species being present alone.
          # for stac = 6 and 7, integration is from the maximum colonization
          # time until the first branching time. This is the same as we did for
          # stac = 1, 5.
          # for stac = 8 and 9, integration is from the maximum colonization
          # time until the minimum colonization time.
          # In all cases we are dealing with a maximum colonization time which
          # means that any colonization that took place before this maximum
          # colonization time (including presence in the non-oceanic scenario)
          # does not count and should be followed by another colonization.
          # To allow this we introduce more sets of equations.
          epss <- 1.01E-5 #We're taking the risk
          if (abs(brts[2] - brts[1]) >= epss)
          {
            if(!is.list(CS_version)) CS_version <- as.list(CS_version)
            function_to_optimize <- CS_version$function_to_optimize
            if(is.null(function_to_optimize)) function_to_optimize <- 'DAISIE'
            if(function_to_optimize == 'DAISIE_approx')
            {
              probs2 <- rep(0, 4 * lx)
              probs2[(2 * lx + 1):(4 * lx)] <- probs[1:(2 * lx)]
              probs2[1:(2 * lx)] <- 0
              #probs2[1:(2 * lx)] <- probs[1:(2 * lx)]
              #probs2[(2 * lx + 1):(4 * lx)] <- 0
              probs <- probs2
              rm(probs2)
              probs <- DAISIE_integrate(probs,brts[2:3],DAISIE_loglik_rhs_time1,c(pars1,k1,ddep),rtol = reltolint,atol = abstolint,method = methode)
              cp <- checkprobs2(lx, loglik, probs, verbose); loglik <- cp[[1]]; probs <- cp[[2]]
              if (stac %in% c(1, 5))
              {
                #loglik <- loglik + log(probs[(stac == 1) * lx + (stac == 5) + 1 + missnumspec])
                loglik <- loglik + log((missnumspec + 1) * probs[(stac == 1) * lx + (stac == 5) + 1 + missnumspec])
              } else if (stac %in% c(6, 7, 8, 9))
              {
                probs2 <- rep(0, 3 * lx)
                probs2[1:(lx - 1)] <- (1:(lx - 1)) * probs[2:lx]
                probs2[(lx + 1):(2 * lx - 1)] <- (1:(lx - 1)) * probs[(lx + 2):(2 * lx)]
                probs2[2 * lx + 1] <- probs[(lx + 1)]
                probs2[(2 * lx + 2):(3 * lx)] <- 0
                probs <- probs2
                rm(probs2)
                if (stac %in% c(8, 9))
                {
                  k1 <- 1
                  probs = DAISIE_integrate(probs,c(brts[3:4]),DAISIE_loglik_rhs_time2,c(pars1,k1,ddep),rtol = reltolint,atol = abstolint,method = methode)
                  cp = checkprobs2(lx, loglik, probs, verbose); loglik = cp[[1]]; probs = cp[[2]]
                  loglik = loglik + log(probs[(stac == 8) * (2 * lx + 1) + (stac == 9) + missnumspec])
                }
              }
            } else
            {
              lx1 <- lx
              lx2 <- lx
              probs2 <- rep(0,lx1 + 2 * lx1 * lx2)
              probs2[1:lx1] <- probs[(lx + 1):(2 * lx)]
              probs2[(lx1 + 1):(lx1 + lx1)] <- probs[1:lx]
              probs <- probs2
              rm(probs2)
              nndd <- nndivdep_CS(lx1 = lx1, lx2 = lx2, K = K, k = 0)
              parslist <- list(pars = pars1, k = 0, ddep = ddep, nndd = nndd)
              probs <- DAISIE_integrate(probs,brts[2:3],DAISIE_loglik_rhs_time3,parslist,rtol = reltolint,atol = abstolint,method = methode)
              probs2 <- rep(0, 3 * lx)
              probs3 <- probs[(lx1 + 1):(lx1 + lx1 * lx2)]
              dim(probs3) <- c(lx1,lx2)
              probs3[1:lx1,1:(lx2 - 1)] <- probs3[1:lx1,2:lx2] * matrix(1:(lx2 - 1),lx1,lx2 - 1,byrow = T)
              probs4 <- probs[(lx1 + lx1 * lx2 + 1):(lx1 + 2 * lx1 * lx2)]
              dim(probs4) <- c(lx1,lx2)
              probs5 <- probs4
              probs4[1:lx1,1:(lx2 - 1)] <- probs4[1:lx1,2:lx2] * matrix(1:(lx2 - 1),lx1,lx2 - 1,byrow = T)
              for(cnt in 2:(lx + 1)) {
                probs2[cnt - 1] <- sum(probs3[row(probs3) + col(probs3) == cnt])
                probs2[lx + cnt - 1] <- sum(probs4[row(probs4) + col(probs4) == cnt])
                probs2[2 * lx + cnt - 1] <- sum(probs5[row(probs5) + col(probs5) == cnt])
              }
              probs <- probs2
              rm(probs2, probs3, probs4, probs5)
              if (stac %in% c(1, 5))
              {
                loglik <- loglik + log(probs[(stac == 1) * 2 * lx + 1 + missnumspec])
              }
            }
          } else
          { #max age equals island age
            probs2 <- rep(0, 4 * lx)
            probs2[1:(2 * lx)] <- probs[1:(2 * lx)]
            probs2[(2 * lx + 1):(4 * lx)] <- 0
            probs <- probs2
            rm(probs2)
            probs <- DAISIE_integrate(probs,brts[2:3],DAISIE_loglik_rhs_time1,c(pars1,k1,ddep),rtol = reltolint,atol = abstolint,method = methode)
            cp <- checkprobs2(lx, loglik, probs, verbose); loglik <- cp[[1]]; probs <- cp[[2]]
            if (stac %in% c(1, 5))
            {
              loglik <- loglik + log((missnumspec + 1) * probs[(stac == 1) * lx + (stac == 5) + 1 + missnumspec])
            } else if (stac %in% c(6, 7, 8, 9))
            {
              probs2 <- rep(0, 3 * lx)
              probs2[1:(lx - 1)] <- (1:(lx - 1)) * probs[2:lx]
              probs2[(lx + 1):(2 * lx - 1)] <- (1:(lx - 1)) * probs[(lx + 2):(2 * lx)]
              probs2[2 * lx + 1] <- probs[(lx + 1)]
              probs2[(2 * lx + 2):(3 * lx)] <- probs[(lx + 2):(2 * lx)]
              probs <- probs2
              rm(probs2)
            }
          }
          if (stac %in% c(8, 9))
          {
            k1 <- 1
            probs <- DAISIE_integrate(probs,c(brts[3:4]),DAISIE_loglik_rhs_time2,c(pars1,k1,ddep),rtol = reltolint,atol = abstolint,method = methode)
            cp <- checkprobs2(lx, loglik, probs, verbose); loglik <- cp[[1]]; probs = cp[[2]]
            loglik <- loglik + log(probs[(stac == 8) * (2 * lx + 1) + (stac == 9) + missnumspec])
          }
        }
        else if (stac %in% c(2, 3, 4) )
        {
          # for stac = 2, 3, 4, integration is then from the colonization
          # event until the first branching time (stac = 2 and 3) or the present
          # (stac = 4). We add a set of equations for Q_M,n, the probability
          # that the process is compatible with the data, and speciation has not
          # happened; during this time immigration is not allowed because it
          # would alter the colonization time.
          t <- brts[2]
          gamvec <- divdepvec(
            lac_or_gam = "gam",
            pars1 = pars1,
            t = t,
            lx = lx,
            k1 = k1,
            ddep = ddep * (ddep == 11 | ddep == 21)
          )
          probs[(2 * lx + 1):(3 * lx)] <- gamvec[1:lx] * probs[1:lx] +
            gamvec[2:(lx + 1)] * probs[(lx + 1):(2 * lx)]
          probs[1:(2 * lx)] <- 0
          k1 <- 1
          probs <- DAISIE_integrate(probs,c(brts[2:3]),DAISIE_loglik_rhs_time2,c(pars1,k1,ddep),rtol = reltolint,atol = abstolint,method = methode)
          cp <- checkprobs2(lx,loglik,probs, verbose); loglik <- cp[[1]]; probs = cp[[2]]
          if (stac == 4)
            # if stac = 4, we're done and we take an element from Q_M,n
          {
            loglik <- loglik + log(probs[2 * lx + 1 + missnumspec])
          }
        }
        if (stac %in% c(2, 3, 6, 7) )
        {
          # at the first branching point all probabilities of states Q_M,n are
          # transferred to probabilities where only endemics are present. Then
          # go through the branching points.
          S1 <- length(brts) - 1
          startk <- 3
          if(S1 >= startk)
          {
            t <- brts[startk]
            lacvec <- divdepvec(
              lac_or_gam = "lac",
              pars1 = pars1,
              t = t,
              lx = lx + stac %in% c(6,7),
              k1 = k1,
              ddep = ddep
            )
            #if(stac %in% c(2,3,6,7)) {
            probs2 <- rep(0, 2 * lx + 1)
            probs2[1:lx] <- lacvec[1:lx] * (probs[1:lx] + probs[(2 * lx + 1):(3 * lx)])
            probs2[(lx + 1):(2 * lx)] <- lacvec[2:(lx + 1)] * probs[(lx + 1):(2 * lx)]
            probs <- probs2
            rm(probs2)
            #} else { # stac in c(6,7)
            #  probs2 <- probs
            #  probs2[1:(lx - 1)] <- lacvec[2:lx] *
            #    ((1:(lx - 1)) * probs[2:lx] + probs[(2 * lx + 1):(3 * lx - 1)])
            #  #probs2[1:(lx - 1)] <- lacvec[2:lx] *
            #  #  ((1:(lx - 1)) * probs[2:lx] + probs[(lx + 1):(2 * lx - 1)])
            #  probs2[(lx + 1):(2 * lx - 1)] <- lacvec[3:(lx + 1)] * (1:(lx - 1)) *
            #    probs[(lx + 2):(2 * lx)]
            #  probs2[lx] <- 0
            #  probs2[2 * lx] <- 0
            #  probs <- probs2
            #  rm(probs2)
            #}
            #probs <- probs[-c((2 * lx + 2):(3 * lx))]
            #probs[2 * lx + 1] <- 0
            # After speciation, colonization is allowed again (re-immigration)
            # all probabilities of states with the immigrant present are set to
            # zero and all probabilities of states with endemics present are
            # transported to the state with the colonist present waiting for
            # speciation to happen. We also multiply by the (possibly diversity-
            # dependent) immigration rate.
            for (k in startk:S1)
            {
              k1 <- k - 1
              probs <- DAISIE_integrate(probs,brts[k:(k+1)],DAISIE_loglik_rhs_time,c(pars1,k1,ddep),rtol = reltolint,atol = abstolint,method = methode)
              cp <- checkprobs2(lx, loglik, probs, verbose); loglik <- cp[[1]]; probs = cp[[2]]
              if(k < S1)
              {
                # speciation event
                t <- brts[k + 1]
                lacvec <- divdepvec(
                  lac_or_gam = "lac",
                  pars1 = pars1,
                  t = t,
                  lx = lx,
                  k1 = k1,
                  ddep = ddep
                )
                probs[1:(2 * lx)] <- c(lacvec[1:lx], lacvec[2:(lx + 1)]) *
                  probs[1:(2 * lx)]
              }
            }
          }
          # we evaluate the probability of the phylogeny with any missing species at the present without (stac = 2 or stac = 6) or with (stac = 3 or stac = 7) the immigrant species
          loglik <- loglik + log(probs[(stac %in% c(3, 7)) * lx + 1 + missnumspec])
        }
      }
    }
  }
  if (length(pars1) == 11) {
    print_parameters_and_loglik(pars = c(stac, pars1[1:10]), # should this be 6:10, or 6:11?
                                loglik = loglik,
                                verbose = pars2[4],
                                type = 'clade_loglik')
  } else {
    print_parameters_and_loglik(pars = c(stac, pars1[1:5]),
                                loglik = loglik,
                                verbose = pars2[4],
                                type = 'clade_loglik')
  }
  if (is.na(loglik)) {
    wrn1 <- paste("NA in loglik encountered for a clade with stac = ",stac, sep = '')
    wrn2 <- paste(" and branching times are:\n", paste(brts, collapse = ", "))
    wrn3 <- paste("\nChanging loglik to -Inf.\n")
    cat(paste(wrn1, wrn2, wrn3, sep = ""))
    loglik <- -Inf
  }
  loglik <- as.numeric(loglik)
  return(loglik)
}
