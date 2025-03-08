#' CS iteration control
#'
#' Sets or retrieves the max. number of iterations used by the odeint solver.
#'
#' @param max_steps \code{max_steps}: sets max. iterations to \code{max_steps}. \cr
#' @return current max. iterations
#'
#' @export DAISIE_CS_max_steps
DAISIE_CS_max_steps <- function(max_steps) {
  return(.Call("daisie_odeint_cs_max_steps", max_steps))
}


#` adams_bashforth and adams_bashforth_moulton integration control
#'
#' Sets or retrieves the factor to calculate the step-size used by the odeint::adams_bashforth[_moulton] solvers.
#'
#' @param factor sets step-size to \code{factor * (t1 - t0)}. \cr
#' @return current factor
#'
#' @export DAISIE_abm_factor
DAISIE_abm_factor <- function(factor) {
  return(.Call("daisie_odeint_abm_factor", factor))
}

DAISIE_loglik_rhs_precomp <- function(pars,lx)
{
  lac = pars[1]
  mu = pars[2]
  K = pars[3]
  gam = pars[4]
  laa = pars[5]
  kk = pars[6]
  ddep = pars[7]

  nn = -2:(lx + 2 * kk + 1)
  lnn = length(nn)
  nn = pmax(rep(0, lnn), nn)

  if(ddep == 0)
  {
    laavec = laa * rep(1,lnn)
    lacvec = lac * rep(1,lnn)
    muvec = mu * rep(1,lnn)
    gamvec = gam * rep(1,lnn)
  } else if(ddep == 1)
  {
    laavec = laa * rep(1,lnn)
    lacvec = pmax(rep(0,lnn),lac * (1 - nn/K))
    muvec = mu * rep(1,lnn)
    gamvec = gam * rep(1,lnn)
  } else if(ddep == 2)
  {
    laavec = laa * rep(1,lnn)
    lacvec = pmax(rep(0,lnn),lac * exp(-nn/K))
    muvec = mu * rep(1,lnn)
    gamvec = gam * rep(1,lnn)
  } else if(ddep == 11)
  {
    laavec = laa * rep(1,lnn)
    #lacvec = pmax(rep(0,lnn),lac * (1 - nn/K))
    lacvec <- get_clado_rate_per_capita(lac = lac,
                                        d = 0,
                                        num_spec = nn,
                                        K = K,
                                        A = 1)
    muvec = mu * rep(1,lnn)
    muvec <- rep(1,lnn) * get_ext_rate_per_capita(mu = mu,
                                                  x = 0)
    #gamvec = pmax(rep(0,lnn),gam * (1 - nn/K))
    gamvec <- get_immig_rate_per_capita(gam = gam,
                                        num_spec = nn,
                                        K = K,
                                        A = 1)

  } else if(ddep == 21)
  {
    laavec = laa * rep(1,lnn)
    lacvec = pmax(rep(0,lnn),lac * exp(-nn/K))
    muvec = mu * rep(1,lnn)
    gamvec = pmax(rep(0,lnn),gam * exp(-nn/K))
  } else if(ddep == 3)
  {
    laavec = laa * rep(1,lnn)
    lacvec = lac * rep(1,lnn)
    muvec = mu * (1 + nn/K)
    gamvec = gam * rep(1,lnn)
  }
  return(c(laavec, lacvec, muvec, gamvec, nn, kk))
}

DAISIE_loglik_rhs <- function(t, x, parsvec) {
  rhs <- 0
  kk <- parsvec[length(parsvec)]
  lx <- (length(x) - 1)/2
  lnn <- lx + 4 + 2 * kk
  laavec <- parsvec[1:lnn]
  lacvec <- parsvec[(lnn + 1):(2 * lnn)]
  muvec <- parsvec[(2 * lnn + 1):(3 * lnn)]
  gamvec <- parsvec[(3 * lnn + 1):(4 * lnn)]
  nn <- parsvec[(4 * lnn + 1):(5 * lnn)]

  xx1 = c(0,0,x[1:lx],0) #Q^0_n
  xx2 = c(0,0,x[(lx + 1):(2 * lx)],0) #Q^{M,0},n
  xx3 = x[2 * lx + 1] #relict

  nil2lx = 3:(lx + 2)

  il1 = nil2lx+kk-1
  il2 = nil2lx+kk+1
  il3 = nil2lx+kk
  il4 = nil2lx+kk-2

  in1 = nil2lx+2*kk-1
  in2 = nil2lx+1
  in3 = nil2lx+kk

  ix1 = nil2lx-1
  ix2 = nil2lx+1
  ix3 = nil2lx
  ix4 = nil2lx-2

  dx1 <- laavec[il1 + 1] * xx2[ix1] +
    lacvec[il4 + 1] * xx2[ix4] +
    muvec[il2 + 1] * xx2[ix3] +
    lacvec[il1] * nn[in1] * xx1[ix1] +
    muvec[il2] * nn[in2] * xx1[ix2] +
    -(muvec[il3] + lacvec[il3]) * nn[in3] * xx1[ix3] +
    -gamvec[il3] * xx1[ix3]
  dx2 <- gamvec[il3] * xx1[ix3] +
    lacvec[il1 + 1] * nn[in1] * xx2[ix1] +
    muvec[il2 + 1] * nn[in2] * xx2[ix2] +
    -(muvec[il3 + 1] + lacvec[il3 + 1]) * nn[in3 + 1] * xx2[ix3] +
    -laavec[il3 + 1] * xx2[ix3]
  dx3 <- 0
  return(list(c(dx1,dx2,dx3)))
}

DAISIE_loglik_rhs1 <- function(t, x, parsvec) {
  rhs <- 1
  kk <- parsvec[length(parsvec)]
  lx <- (length(x))/4
  lnn <- lx + 4 + 2 * kk
  laavec <- parsvec[1:lnn]
  lacvec <- parsvec[(lnn + 1):(2 * lnn)]
  muvec <- parsvec[(2 * lnn + 1):(3 * lnn)]
  gamvec <- parsvec[(3 * lnn + 1):(4 * lnn)]
  nn <- parsvec[(4 * lnn + 1):(5 * lnn)]

  xx1 <- c(0,0,x[1:lx],0) #Q^0_n
  xx2 <- c(0,0,x[(lx + 1):(2 * lx)],0) #Q^{M,0}_n
  xx3 <- c(0,0,x[(2 * lx + 1):(3 * lx)],0) #Q^0_{M,n}
  xx4 <- c(0,0,x[(3 * lx + 1):(4 * lx)],0) #Q^{M,0}_{M,n}

  nil2lx <- 3:(lx + 2)

  il1 <- nil2lx+kk-1
  il2 <- nil2lx+kk+1
  il3 <- nil2lx+kk
  il4 <- nil2lx+kk-2

  in1 <- nil2lx+2*kk-1
  in2 <- nil2lx+1
  in3 <- nil2lx+kk
  in4 <- nil2lx-1

  ix1 <- nil2lx-1
  ix2 <- nil2lx+1
  ix3 <- nil2lx
  ix4 <- nil2lx-2

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
    -(muvec[il3 + 1] + lacvec[il3 + 1]) * nn[in3 + 1] * xx2[ix3] +
    -laavec[il3 + 1] * xx2[ix3]

  dx3 <- lacvec[il1] * nn[in1] * xx3[ix1] +
    laavec[il1 + 1] * xx4[ix1] +
    lacvec[il4 + 1] * xx4[ix4] +
    muvec[il2] * nn[in2] * xx3[ix2] +
    muvec[il3 + 1] * xx4[ix3] +
    -(lacvec[il3] + muvec[il3]) * nn[in3] * xx3[ix3] +
    -gamvec[il3] * xx3[ix3]

  dx4 <- lacvec[il1 + 1] * nn[in1] * xx4[ix1] +
    muvec[il2 + 1] * nn[in2] * xx4[ix2] +
    -(laavec[il3 + 1] * lacvec[il3 + 1] + muvec[il3 + 1]) * nn[in3 + 1] * xx4[ix3] +
    -gamvec[il3 + 1] * xx4[ix3]

  return(list(c(dx1,dx2,dx3,dx4)))
}

DAISIE_loglik_rhs2 <- function(t, x, parsvec) {
  rhs <- 2
  kk <- parsvec[length(parsvec)]
  lx <- (length(x))/3
  lnn <- lx + 4 + 2 * kk
  laavec <- parsvec[1:lnn]
  lacvec <- parsvec[(lnn + 1):(2 * lnn)]
  muvec <- parsvec[(2 * lnn + 1):(3 * lnn)]
  gamvec <- parsvec[(3 * lnn + 1):(4 * lnn)]
  nn <- parsvec[(4 * lnn + 1):(5 * lnn)]

  xx1 = c(0,0,x[1:lx],0) #Q^1_n
  xx2 = c(0,0,x[(lx + 1):(2 * lx)],0) #Q^{M,1}_n
  xx3 = c(0,0,x[(2 * lx + 1):(3 * lx)],0) # Q^1_{M,n}

  nil2lx = 3:(lx + 2)

  il1 = nil2lx+kk-1
  il2 = nil2lx+kk+1
  il3 = nil2lx+kk
  il4 = nil2lx+kk-2

  in1 = nil2lx+2*kk-1
  in2 = nil2lx+1
  in3 = nil2lx+kk
  in4 = nil2lx-1

  ix1 = nil2lx-1
  ix2 = nil2lx+1
  ix3 = nil2lx
  ix4 = nil2lx-2

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
  dx3 <- lacvec[il1] * nn[in4] * xx3[ix1] +
    muvec[il2] * nn[in2] * xx3[ix2] +
    -(lacvec[il3] + muvec[il3]) * nn[in3] * xx3[ix3] +
    -(laavec[il3] + gamvec[il3]) * xx3[ix3]

  return(list(c(dx1,dx2,dx3)))
}

checkprobs <- function(lv, loglik, probs, verbose) {
  probs <- probs * (probs > 0)
  if (is.na(sum(probs[1:lv])) || is.nan(sum(probs))) {
    loglik <- -Inf
  } else if (sum(probs[1:lv]) <= 0) {
    loglik <- -Inf
  } else {
    loglik <- loglik + log(sum(probs[1:lv]))
    probs[1:lv] <- probs[1:lv] / sum(probs[1:lv])
  }
  if (verbose) {
    message("Numerical issues encountered.")
  }
  return(list(loglik, probs))
}

checkprobs2 <- function(lv, loglik, probs, verbose) {
  probs <- probs * (probs > 0)
  if (is.na(sum(probs)) || is.nan(sum(probs))) {
    loglik <- -Inf
  } else if (sum(probs) <= 0) {
    loglik <- -Inf
  } else {
    sp <- sum(sort(probs))
    loglik = loglik + log(sp)
    probs = probs/sp
  }
  if (verbose) {
    message("Numerical issues encountered \n")
  }
  return(list(loglik, probs))
}

divdepvec <- function(lac_or_gam,
                      pars1,
                      t = NULL,
                      lx,
                      k1,
                      ddep) {
  island_ontogeny <- pars1[15]
  if (!is.na(island_ontogeny)) {

    # ddep is NOT yet used in the time dependent case

    # lac0 <- parsvec[1]
    # mu0 <- parsvec[2]
    # K0 <- parsvec[3]
    # gam0 <- parsvec[4]
    # laa0 <- parsvec[5]
    # d <- parsvec[6]
    # x <- parsvec[7]
    # area_pars <- parsvec[8:14]
    # island_ontogeny <- parsvec[15]
    # sea_level <- parsvec[16]
    # total_time <- parsvec[17]
    # peak <- parsvec[18]
    area <- island_area_vector(
      timeval = abs(t),
      area_pars = pars1[8:14],
      island_ontogeny = island_ontogeny,
      sea_level = pars1[16],
      total_time = pars1[17],
      peak = pars1[18]
    )
    if (lac_or_gam == "lac") {
      divdepvector <- get_clado_rate_per_capita(
        lac = pars1[1],
        d = pars1[6],
        num_spec = (0:lx) + k1,
        K = pars1[3],
        A = area
      )
    } else if (lac_or_gam == "gam") {
      divdepvector <- get_immig_rate_per_capita(
        gam = pars1[4],
        num_spec = (0:lx) + k1,
        K = pars1[3],
        A = area
      )
    }
  } else {
    lacgam <- ifelse(lac_or_gam == "lac", pars1[1], pars1[4])
    divdepvector <- divdepvec1(
      lacgam = lacgam,
      K = pars1[3],
      lx = lx,
      k1 = k1,
      ddep = ddep
    )
    #divdepvector <- get_immig_rate_per_capita(gam = lacgam,
    #                                          num_spec = k1 + (0:lx),
    #                                          K = pars1[3],
    #                                          A = 1)
  }
  return(divdepvector)
}

divdepvec1 <- function(lacgam, K, lx, k1, ddep) {
  if (ddep == 1 | ddep == 11) {
    vec <- pmax(rep(0, lx + 1), lacgam * (1 - ( (0:lx) + k1) / K))
  } else {
    if (ddep == 2 | ddep == 21) {
      vec <- pmax(rep(0, lx + 1), lacgam * exp(-((0:lx) + k1) / K))
    } else {
      if (ddep == 0 | ddep == 3) {
        vec <- lacgam * rep(1, lx + 1)
      }
    }
  }
  return(vec)
}

DAISIE_loglik_CS_M1 <- DAISIE_loglik <- function(pars1,
                                                 pars2,
                                                 brts,
                                                 stac,
                                                 missnumspec,
                                                 methode = "lsodes",
                                                 abstolint = 1E-16,
                                                 reltolint = 1E-10,
                                                 verbose) {
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
  if (any(is.infinite(pars1)) ) {
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
    loglik = -Inf
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
    message("One or more parameters are negative.")
    loglik <- -Inf
    return(loglik)
  }
  if ((ddep == 1 | ddep == 11) & ceiling(K) < (S + missnumspec)) {
    if (verbose) {
      message('The proposed value of K is incompatible with the number of species
          in the clade. Likelihood for this parameter set
          will be set to -Inf. \n')
    }
    loglik <- -Inf
    return(loglik)
  }
  lac <- pars1[1]
  if(lac == Inf & missnumspec == 0 & length(pars1) == 5) {
      if(verbose) warning('Infinite lambda detected')
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
      probs = DAISIE_integrate(probs,brts[1:2],DAISIE_loglik_rhs,c(pars1,k1,ddep),rtol = reltolint,atol = abstolint,method = methode)
      cp = checkprobs2(lv = 2 * lx, loglik, probs, verbose); loglik = cp[[1]]; probs = cp[[2]]
      if(stac == 0)
      {
        # for stac = 0, the integration was from the origin of the island until
        # the present so we can immediately evaluate the probability of no clade
        # being present and no immigrant species.
        loglik = loglik + log(probs[1])
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
          # To allow this we introduce a third and fourth set of equations for
          # the probability that colonization might have happened before but
          # recolonization has not taken place yet (Q_M,n and Q^M_{M,n}).
          epss <- 1.01E-5 #We're taking the risk
          if (abs(brts[2] - brts[1]) >= epss) {
            probs[(2 * lx + 1):(4 * lx)] <- probs[1:(2 * lx)]
            probs[1:(2 * lx)] <- 0
          } else { #max age equals island age
            probs[(2 * lx + 1):(4 * lx)] <- 0
          }

          probs <- DAISIE_integrate(probs,brts[2:3],DAISIE_loglik_rhs1,c(pars1,k1,ddep),rtol = reltolint,atol = abstolint,method = methode)
          cp <- checkprobs2(lx, loglik, probs, verbose); loglik <- cp[[1]]; probs <- cp[[2]]
          if (stac %in% c(1, 5))
          {
            loglik <- loglik + log(probs[(stac == 1) * lx + (stac == 5) + 1 + missnumspec])
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
              probs = DAISIE_integrate(probs,c(brts[3:4]),DAISIE_loglik_rhs2,c(pars1,k1,ddep),rtol = reltolint,atol = abstolint,method = methode)
              cp = checkprobs2(lx, loglik, probs, verbose); loglik = cp[[1]]; probs = cp[[2]]
              loglik = loglik + log(probs[(stac == 8) * (2 * lx + 1) + (stac == 9) + missnumspec])
            }
          }
        } else if (stac %in% c(2, 3, 4) )
        {
          # for stac = 2, 3, 4, integration is then from the colonization
          # event until the first branching time (stac = 2 and 3) or the present
          # (stac = 4). We add a set of equations for Q_M,n, the probability
          # that the process is compatible with the data, and speciation has not
          # happened; during this time immigration is not allowed because it
          # would alter the colonization time.
          t <- brts[2]
          gamvec = divdepvec(
            lac_or_gam = "gam",
            pars1 = pars1,
            t = t,
            lx = lx,
            k1 = k1,
            ddep = ddep * (ddep == 11 | ddep == 21)
          )
          probs[(2 * lx + 1):(3 * lx)] = gamvec[1:lx] * probs[1:lx] +
            gamvec[2:(lx + 1)] * probs[(lx + 1):(2 * lx)]
          probs[1:(2 * lx)] = 0
          k1 <- 1
          probs = DAISIE_integrate(probs,c(brts[2:3]),DAISIE_loglik_rhs2,c(pars1,k1,ddep),rtol = reltolint,atol = abstolint,method = methode)
          cp = checkprobs2(lx,loglik,probs, verbose); loglik = cp[[1]]; probs = cp[[2]]
          if (stac == 4)
            # if stac = 4, we're done and we take an element from Q_M,n
          {
            loglik = loglik + log(probs[2 * lx + 1 + missnumspec])
          }
        }
        if (stac %in% c(2, 3, 6, 7) )
        {
          # at the first branching point all probabilities of states Q_M,n are
          # transferred to probabilities where only endemics are present. Then
          # go through the branching points.
          S1 = length(brts) - 1
          startk = 3
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
            #if(stac %in% c(2,3,6,7))
            #{
              probs[1:lx] <- lacvec[1:lx] * (probs[1:lx] + probs[(2 * lx + 1):(3 * lx)])
              probs[(lx + 1):(2 * lx)] <- lacvec[2:(lx + 1)] * probs[(lx + 1):(2 * lx)]
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
            probs <- probs[-c((2 * lx + 2):(3 * lx))]
            probs[2 * lx + 1] <- 0
            # After speciation, colonization is allowed again (re-immigration)
            # all probabilities of states with the immigrant present are set to
            # zero and all probabilities of states with endemics present are
            # transported to the state with the colonist present waiting for
            # speciation to happen. We also multiply by the (possibly diversity-
            # dependent) immigration rate.
            for (k in startk:S1)
            {
              k1 <- k - 1
              probs <- DAISIE_integrate(probs,brts[k:(k+1)],DAISIE_loglik_rhs,c(pars1,k1,ddep),rtol = reltolint,atol = abstolint,method = methode)
              cp <- checkprobs2(lx, loglik, probs, verbose); loglik = cp[[1]]; probs = cp[[2]]
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

 if (length(pars1) == 11) { # CHANGE
    print_parameters_and_loglik(pars = c(stac, pars1[5:10]), # should this be 6:10, or 6:11?
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
    message("NA in loglik encountered. Changing to -Inf.")
    loglik <- -Inf
  }
  loglik <- as.numeric(loglik)
  #testit::assert(is.numeric(loglik))
  return(loglik)
}

DAISIE_loglik_CS_choice <- function(
    pars1,
    pars2,
    datalist = NULL,
    brts,
    stac,
    missnumspec,
    methode = "lsodes",
    CS_version = 1,
    abstolint = 1E-16,
    reltolint = 1E-10,
    verbose = FALSE
)
{
  if (CS_version[[1]] == 1) {
    loglik <- DAISIE_loglik(
      pars1 = pars1,
      pars2 = pars2,
      brts = brts,
      stac = stac,
      missnumspec = missnumspec,
      methode = methode,
      abstolint = abstolint,
      reltolint = reltolint,
      verbose = verbose
    )
  } else if (CS_version[[1]] == 2) {
    loglik <- DAISIE_loglik_integrate(
      pars1 = pars1,
      pars2 = pars2,
      brts = brts,
      stac = stac,
      missnumspec = missnumspec,
      methode = methode,
      CS_version = CS_version,
      abstolint = abstolint,
      reltolint = reltolint,
      verbose = verbose
    )
  } else if (CS_version[[1]] == 0) {
    loglik <- DAISIE_loglik_IW_M1(
      pars1 = pars1,
      pars2 = pars2,
      datalist = datalist,
      brts = brts,
      stac = stac,
      missnumspec = missnumspec,
      methode = methode,
      abstolint = abstolint,
      reltolint = reltolint,
      verbose = verbose
    )
  }
  return(loglik)
}

approximate_logp0 <- function(gamma, mu, t)
{
  logp0 <- -log(mu + gamma) + log(mu + gamma * exp(-(mu + gamma) * t))
  return(logp0)
}

#' @name DAISIE_loglik_CS
#' @aliases DAISIE_loglik_all DAISIE_loglik_CS
#' @title Computes the loglikelihood of the DAISIE model with clade-specific
#' diversity-dependence given data and a set of model parameters
#' @description Computes the loglikelihood of the DAISIE model with clade-specific
#' diversity-dependence given colonization and branching times for lineages on
#' an island, and a set of model parameters. The output is a loglikelihood value
#' @inheritParams default_params_doc
#' @param pars1 Contains the model parameters: \cr \cr
#' \code{pars1[1]} corresponds to lambda^c (cladogenesis rate) \cr
#' \code{pars1[2]} corresponds to mu (extinction rate) \cr
#' \code{pars1[3]} corresponds to K (clade-level carrying capacity) \cr
#' \code{pars1[4]} corresponds to gamma (immigration rate) \cr
#' \code{pars1[5]} corresponds to lambda^a (anagenesis rate) \cr
#' \code{pars1[6]} corresponds to lambda^c (cladogenesis rate) for an optional
#' subset of the species \cr
#' \code{pars1[7]} corresponds to mu (extinction rate) for an optional subset of the species\cr
#' \code{pars1[8]} corresponds to K (clade-level carrying capacity) for an optional subset of the
#' species\cr
#' \code{pars1[9]} corresponds to gamma (immigration rate) for an optional subset of the species\cr
#' \code{pars1[10]} corresponds to lambda^a (anagenesis rate) for an optional subset of the species\cr
#' \code{pars1[11]} corresponds to p_f (fraction of mainland species that belongs to the second
#' subset of species\cr
#' The elements 6:10 and 11 are optional, that is, pars1
#' should either contain 5, 10 or 11 elements. If 10, then the fraction of
#' potential colonists of type 2 is computed from the data. If 11, then
#' pars1[11] is used, overruling any information in the data.
#' @param pars2 Contains the model settings \cr \cr
#' \code{pars2[1]} corresponds to lx = length of ODE variable x \cr
#' \code{pars2[2]} corresponds to ddmodel = diversity-dependent model, model of diversity-dependence, which can be one
#' of\cr \cr
#' ddmodel = 0 : no diversity dependence \cr
#' ddmodel = 1 : linear dependence in speciation rate \cr
#' ddmodel = 11: linear dependence in speciation rate and in immigration rate \cr
#' ddmodel = 2 : exponential dependence in speciation rate\cr
#' ddmodel = 21: exponential dependence in speciation rate and in immigration rate\cr\cr
#' \code{pars2[3]} corresponds to cond = setting of conditioning\cr \cr
#' cond = 0 : conditioning on island age \cr
#' cond = 1 : conditioning on island age and non-extinction of the island biota \cr \cr
#' cond > 1 : conditioning on island age and having at least cond colonizations on the island \cr \cr
#' \code{pars2[4]} sets the level of verbosity. When equal to 0, no output is generated. At higher values
#' (1 or 2) more output will be generated.
#' @param datalist Data object containing information on colonisation and
#' branching times. This object can be generated using the DAISIE_dataprep
#' function, which converts a user-specified data table into a data object, but
#' the object can of course also be entered directly. It is an R list object
#' with the following elements.\cr The first element of the list has two or
#' three components:
#' \cr \cr \code{$island_age} - the island age \cr
#' Then, depending on whether a distinction between types is made, we have:\cr
#' \code{$not_present} - the number of mainland lineages that are not present
#' on the island \cr
#' or:\cr
#' \code{$not_present_type1} - the number of mainland lineages of type 1 that are not present on the island \cr
#' \code{$not_present_type2} - the number of mainland lineages of type 2 that
#' are not present on the island \cr \cr
#' The remaining elements of the list
#' each contains information on a single colonist lineage on the island and has
#' 5 components:\cr \cr
#' \code{$colonist_name} - the name of the species or
#' clade that colonized the island \cr
#' \code{$branching_times} - island age and
#' stem age of the population/species in the case of Non-endemic,
#' Non-endemic_MaxAge and Endemic anagenetic species. For cladogenetic species
#' these should be island age and branching times of the radiation including
#' the stem age of the radiation.\cr
#' \code{$stac} - the status of the colonist \cr \cr
#' * Non_endemic_MaxAge: 1 \cr
#' * Endemic: 2 \cr
#' * Endemic&Non_Endemic: 3 \cr
#' * Non_Endemic: 4 \cr
#' * Endemic_Singleton_MaxAge: 5 \cr
#' * Endemic_Clade_MaxAge: 6 \cr
#' * Endemic&Non_Endemic_Clade_MaxAge: 7 \cr \cr
#' * Non_endemic_MaxAge_MinAge: 8 \cr
#' * Endemic_Singleton_MaxAge_MinAge: 9 \cr
#' \code{$missing_species} - number of island species that were not sampled for
#' particular clade (only applicable for endemic clades) \cr
#' \code{$type1or2} - whether the colonist belongs to type 1 or type 2 \cr
#' @param methode Method of the ODE-solver. See package deSolve for details.
#' Default is "lsodes"
#' @param abstolint Absolute tolerance of the integration
#' @param reltolint Relative tolerance of the integration
#' @return The loglikelihood
#' @author Rampal S. Etienne & Bart Haegeman
#' @seealso \code{\link{DAISIE_ML}}, \code{\link{DAISIE_sim_cr}},
#' \code{\link{DAISIE_sim_time_dep}},
#' \code{\link{DAISIE_sim_cr_shift}}
#' @references Valente, L.M., A.B. Phillimore and R.S. Etienne (2015).
#' Equilibrium and non-equilibrium dynamics simultaneously operate in the
#' Galapagos islands. Ecology Letters 18: 844-852.
#' @keywords internal
#' @examples
#'
#' utils::data(Galapagos_datalist_2types)
#' pars1 = c(0.195442017,0.087959583,Inf,0.002247364,0.873605049,
#'           3755.202241,8.909285094,14.99999923,0.002247364,0.873605049,0.163)
#' pars2 = c(100,11,0,1)
#' DAISIE_loglik_all(pars1,pars2,Galapagos_datalist_2types)
#'
#' @export DAISIE_loglik_CS
#' @export DAISIE_loglik_all
DAISIE_loglik_CS <- DAISIE_loglik_all <- function(
    pars1,
    pars2,
    datalist,
    methode = "lsodes",
    CS_version = 1,
    abstolint = 1E-16,
    reltolint = 1E-10) {
  if (length(pars1) == 14) {
    if (datalist[[1]]$island_age > pars1[11]) {
      stop(
        "The island age in the area parameters is inconsistent with the island
        data."
      )
    }
    peak <- calc_peak(
      total_time = datalist[[1]]$island_age,
      area_pars = create_area_pars(
        max_area = pars1[8],
        current_area = pars1[9],
        proportional_peak_t = pars1[10],
        total_island_age = pars1[11],
        sea_level_amplitude = pars1[12],
        sea_level_frequency = pars1[13],
        island_gradient_angle = pars1[14]
      )
    )
    pars1 <- c(
      pars1,
      island_ontogeny = pars2[5],
      sea_level = pars2[6],
      datalist[[1]]$island_age,
      peak
    )
  }

  pars1 <- as.numeric(pars1)
  cond <- pars2[3]
  if (length(pars1) == 6) {
    endpars1 <- 6
  } else {
    endpars1 <- 5
  }

  if(length(pars1) %in% c(5,6) | !is.na(pars2[5])) {
    if(!is.na(pars2[5]))
    {
      endpars1 <- length(pars1)
    }
    logp0 <- DAISIE_loglik_CS_choice(
      pars1 = pars1,
      pars2 = pars2,
      brts = datalist[[1]]$island_age,
      stac = 0,
      missnumspec = 0,
      methode = methode,
      CS_version = CS_version,
      abstolint = abstolint,
      reltolint = reltolint
    )
    if(logp0 >= 0 & pars1[2]/pars1[1] > 100)
    {
      logp0 <- approximate_logp0(gamma = pars1[4],
                                 mu = pars1[2],
                                 t = datalist[[1]]$island_age)
    }
    if(logp0 >= 0)
    {
      message('Positive values of loglik encountered without possibility for approximation. Setting loglik to -Inf.')
      loglik <- -Inf
      print_parameters_and_loglik(pars = pars,
                                  loglik = loglik,
                                  verbose = pars2[4],
                                  type = 'island_loglik')
      return(loglik)
    }
    if (is.null(datalist[[1]]$not_present)) {
      loglik <- (datalist[[1]]$not_present_type1 +
                   datalist[[1]]$not_present_type2) * logp0
      numimm <- (datalist[[1]]$not_present_type1 +
                   datalist[[1]]$not_present_type2) + length(datalist) - 1
    } else {
      loglik <- datalist[[1]]$not_present * logp0
      numimm <- datalist[[1]]$not_present + length(datalist) - 1
    }
    logcond <- logcondprob(numcolmin = cond,numimm = numimm,logp0 = logp0)
    if (length(datalist) > 1) {
      for (i in 2:length(datalist)) {
        datalist[[i]]$type1or2 <- 1
      }
    }
  } else {
    numimm <- datalist[[1]]$not_present_type1 +
      datalist[[1]]$not_present_type2 + length(datalist) - 1
    numimm_type2 <- length(
      which(unlist(datalist)[which(names(unlist(datalist)) == "type1or2")] == 2)
    )
    numimm_type1 <- length(datalist) - 1 - numimm_type2

    if (is.na(pars1[11]) == FALSE && length(pars1) == 11) {
      if (pars1[11] < numimm_type2 / numimm |
          pars1[11] > (1 - numimm_type1 / numimm)) {
        return(-Inf)
      }
      datalist[[1]]$not_present_type2 <- max(
        0,
        round(pars1[11] * numimm) - numimm_type2
      )
      datalist[[1]]$not_present_type1 <- numimm - (length(datalist) - 1) -
        datalist[[1]]$not_present_type2
    }
    logp0_type1 <- DAISIE_loglik_CS_choice(
      pars1 = pars1[1:5],
      pars2 = pars2,
      brts = datalist[[1]]$island_age,
      stac = 0,
      missnumspec = 0,
      methode = methode,
      CS_version = CS_version,
      abstolint = abstolint,
      reltolint = reltolint
    )
    if(logp0_type1 >= 0 & pars1[2]/pars1[1] > 100)
    {
      logp0_type1 <- approximate_logp0(gamma = pars1[4], mu = pars1[2], t = datalist[[1]]$island_age)
    }
    logp0_type2 <- DAISIE_loglik_CS_choice(
      pars1 = pars1[6:10],
      pars2 = pars2,
      brts = datalist[[1]]$island_age,
      stac = 0,
      missnumspec = 0,
      methode = methode,
      CS_version = CS_version,
      abstolint = abstolint,
      reltolint = reltolint
    )
    if(logp0_type2 >= 0 & pars1[7]/pars1[6] > 100)
    {
      logp0_type2 <- approximate_logp0(gamma = pars1[9], mu = pars1[7], t = datalist[[1]]$island_age)
    }
    loglik <- datalist[[1]]$not_present_type1 * logp0_type1 +
      datalist[[1]]$not_present_type2 * logp0_type2
    #logcond <- (cond == 1) *
    #  log(1 - exp((datalist[[1]]$not_present_type1 + numimm_type1) *
    #                logp0_type1 +
    #                (datalist[[1]]$not_present_type2 + numimm_type2) *
    #                logp0_type2))
    logcond <- logcondprob(numcolmin = cond,
                           numimm = c(datalist[[1]]$not_present_type1 + numimm_type1,datalist[[1]]$not_present_type2 + numimm_type2),
                           logp0 = c(logp0_type1,logp0_type2) )
  }
  loglik <- loglik - logcond

  if(length(datalist) > 1)
  {
    for(i in 2:length(datalist))
    {
      if(datalist[[i]]$type1or2 == 1)
      {
        pars <- pars1[1:endpars1]
      } else {
        pars <- pars1[6:10]
      }
      loglik <- loglik + DAISIE_loglik_CS_choice(
        pars1 = pars,
        pars2 = pars2,
        datalist = datalist[[i]],
        brts = datalist[[i]]$branching_times,
        stac = datalist[[i]]$stac,
        missnumspec = datalist[[i]]$missing_species,
        methode = methode,
        CS_version = CS_version,
        abstolint = abstolint,
        reltolint = reltolint)
    }
  }

  print_parameters_and_loglik(pars = pars,
                              loglik = loglik,
                              verbose = pars2[4],
                              parnames = c("lambda^c", "mu", "K", "gamma", "lambda^a", "prob_init_pres"),
                              type = 'island_loglik')
  return(loglik)
}

print_parameters_and_loglik <- function(pars,
                                        loglik,
                                        verbose,
                                        parnames = c("lambda^c", "mu", "K", "gamma", "lambda^a"),
                                        type = 'island_loglik',
                                        distance_dep = NULL)
{
  if (isTRUE(verbose >= 2)) {
    if(type == 'clade_loglik') {
      s1a <- sprintf("Status of colonist: %d", pars[1])
      s1b <- sprintf("Parameters:")
      s1 <- paste(s1a, s1b, sep = ', ')
      s2 <- paste(sprintf("%f", pars[-1]), collapse = ', ')
      s12 <- paste(s1, s2, collapse = ' ')
      s3 <- paste(sprintf("Loglikelihood: %f", loglik), collapse = '')
      message(paste(s12, s3, sep = ', '))
    } else {
      if(type == 'global_ML') {
        s1 <- s1output(pars, distance_dep)
        s3 <- sprintf("Maximum Loglikelihood: %f", loglik)
        message(paste(s1, s3, sep = '\n'))
      } else {
        if(is.null(ncol(pars))) {
          lpars <- length(pars)
        } else {
          lpars <- ncol(pars)
        }
        if(lpars != length(parnames))
        {
          warning('The vectors of parameters and parameter names have different lengths.')
          parnames <- NULL
        } else {
          parnames <- paste(parnames, collapse = ', ')
        }
        if(type == 'island_ML') {
          s1 <- sprintf("Maximum likelihood parameters: ")
          s2 <- paste(sprintf("%f", pars), collapse = ', ')
          s3 <- sprintf("Maximum Loglikelihood: %f", loglik)
          message(paste(s1, parnames, s2, s3, sep = '\n'))
        } else {
          if(type == 'multiple_island_ML') {
            s1 <- sprintf("Maximum likelihood parameters: ")
            s2 <- parnames
            for(i in 1:nrow(pars)) {
               s2 <- paste(s2,paste(sprintf("%f", pars[i,]), collapse = ', '), sep = '\n')
            }
            s3 <- sprintf("Maximum Loglikelihood: %f", loglik)
            message(paste(s1, s2, s3, sep = '\n'))
          } else {
            if(type == 'island_loglik') {
              s1 <- sprintf("Parameters: ")
              s2 <- paste(sprintf("%f", pars), collapse = ', ')
              s3 <- sprintf("Loglikelihood: %f", loglik)
              message(paste(s1, parnames, s2, s3, sep = '\n'))
            } else {
              stop('Type of printing output unknown')
            }
          }
        }
      }
    }
  }
}

s1output <- function(MLpars1,distance_dep)
{
  s1 <- switch(distance_dep,
               power = sprintf('Maximum likelihood parameter estimates:\n
               lambda_c = %f * A^ %f\n
               mu = %f * A^ -%f\n
               K = %f * A^ %f\n
               M * gamma = %f * d^ -%f\n
               lambda_a = %f * d^ %f\n',MLpars1[1],MLpars1[2],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7],MLpars1[8],MLpars1[9],MLpars1[10]),
               sigmoidal_col = sprintf('Maximum likelihood parameter estimates:\n
               lambda_c = %f * A^ %f\n
               mu = %f * A^ -%f\n
               K = %f * A^ %f\n
               M * gamma = %f * (1 - (d/%f)^%f / (1 + (d/%f)^%f )\n
               lambda_a = %f * d^ %f\n',MLpars1[1],MLpars1[2],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7],MLpars1[11],MLpars1[8],MLpars1[11],MLpars1[8],MLpars1[9],MLpars1[10]),
               sigmoidal_ana = sprintf('Maximum likelihood parameter estimates:\n
               lambda_c = %f * A^ %f\n
               mu = %f * A^ -%f\n
               K = %f * A^ %f\n
               M * gamma = %f * d^ -%f\n
               lambda_a = %f * (d/%f)^%f / (1 + (d/%f)^%f )\n',MLpars1[1],MLpars1[2],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7],MLpars1[8],MLpars1[9],MLpars1[11],MLpars1[10],MLpars1[11],MLpars1[10]),
               sigmoidal_clado = sprintf('Maximum likelihood parameter estimates:\n
               lambda_c = %f * (d/%f)^%f / (1 + (d/%f)^%f )\n
               mu = %f * A^ -%f\n
               K = %f * A^ %f\n
               M * gamma = %f * d^ -%f\n
               lambda_a = %f * d^ %f\n',MLpars1[1],MLpars1[11],MLpars1[2],MLpars1[11],MLpars1[2],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7],MLpars1[8],MLpars1[9],MLpars1[10]),
               sigmoidal_col_ana = sprintf('Maximum likelihood parameter estimates:\n
               lambda_c = %f * A^ %f\n
               mu = %f * A^ -%f\n
               K = %f * A^ %f\n
               M * gamma = %f * (1 - (d/%f)^%f / (1 + (d/%f)^%f )\n
               lambda_a = %f * (d/%f)^%f / (1 + (d/%f)^%f )\n',MLpars1[1],MLpars1[2],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7],MLpars1[11],MLpars1[8],MLpars1[11],MLpars1[8],MLpars1[9],MLpars1[12],MLpars1[10],MLpars1[12],MLpars1[10]),
               area_additive_clado = sprintf('Maximum likelihood parameter estimates:\n
               lambda_c = %f * A^ %f * d^ %f \n
               mu = %f * A^ -%f\n
               K = %f * A^ %f\n
               M * gamma = %f * d^ -%f\n
               lambda_a = %f * d^ %f\n',MLpars1[1],MLpars1[2],MLpars1[11],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7],MLpars1[8],MLpars1[9],MLpars1[10]),
               area_interactive_clado = sprintf('Maximum likelihood parameter estimates:\n
               lambda_c = %f * A^ (%f + %f * log(d)) \n
               mu = %f * A^ -%f\n
               K = %f * A^ %f\n
               M * gamma = %f * d^ -%f\n
               lambda_a = %f * d^ %f\n',MLpars1[1],MLpars1[2],MLpars1[11],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7],MLpars1[8],MLpars1[9],MLpars1[10]),
               area_interactive_clado0 = sprintf('Maximum likelihood parameter estimates:\n
               lambda_c = %f * A^ (%f + %f * log(d)) \n
               mu = %f * A^ -%f\n
               K = %f * A^ %f\n
               M * gamma = %f * d^ -%f\n
               lambda_a = %f * d^ %f\n',MLpars1[1],MLpars1[2],MLpars1[11],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7],MLpars1[8],MLpars1[9],MLpars1[10]),
               area_interactive_clado1 = sprintf('Maximum likelihood parameter estimates:\n
               lambda_c = %f * A^ (%f + d/%f)) \n
               mu = %f * A^ -%f\n
               K = %f * A^ %f\n
               M * gamma = %f * d^ -%f\n
               lambda_a = %f * d^ %f\n',MLpars1[1],MLpars1[2],MLpars1[11],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7],MLpars1[8],MLpars1[9],MLpars1[10]),
               area_interactive_clado2 = sprintf('Maximum likelihood parameter estimates:\n
               lambda_c = %f * A^ (%f + d/(d + %f)) \n
               mu = %f * A^ -%f\n
               K = %f * A^ %f\n
               M * gamma = %f * d^ -%f\n
               lambda_a = %f * d^ %f\n',MLpars1[1],MLpars1[2],MLpars1[11],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7],MLpars1[8],MLpars1[9],MLpars1[10]),
               area_interactive_clado3 = sprintf('Maximum likelihood parameter estimates:\n
               lambda_c = %f * (A + d/%f)^ %f\n \n
               mu = %f * A^ -%f\n
               K = %f * A^ %f\n
               M * gamma = %f * d^ -%f\n
               lambda_a = %f * d^ %f\n',MLpars1[1],MLpars1[11],MLpars1[2],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7],MLpars1[8],MLpars1[9],MLpars1[10]),
               area_interactive_clado4 = sprintf('Maximum likelihood parameter estimates:\n
               lambda_c = %f * A^ (%f * d/(d + %f)) \n
               mu = %f * A^ -%f\n
               K = %f * A^ %f\n
               M * gamma = %f * d^ -%f\n
               lambda_a = %f * d^ %f\n',MLpars1[1],MLpars1[2],MLpars1[11],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7],MLpars1[8],MLpars1[9],MLpars1[10])
  )
  return(s1)
}

DAISIE_integrate <- function(initprobs,
                             tvec,
                             rhs_func,
                             pars,
                             rtol,
                             atol,
                             method) {
  if (length(pars) <= 7) {
    return(DAISIE_integrate_const(
      initprobs,
      tvec,
      rhs_func,
      pars,
      rtol,
      atol,
      method)
    )
  } else {
    return(DAISIE_integrate_time(
      initprobs,
      tvec,
      rhs_func,
      pars,
      rtol,
      atol,
      method)
    )
  }
}

DAISIE_integrate_const <- function(initprobs,tvec,rhs_func,pars,rtol,atol,method)
{
  # During code coverage, 'function_as_text' may become:
  #
  # if (TRUE) {
  #   covr::count("DAISIE_loglik_CS.R:58:3:58:25:3:25:4657:4657")
  #   lx <- (length(x) - 1)/2
  # }
  #
  # It is the 'lx <- [something]' part that we are interested in
  #
  # Use a regular expression to extract if the part that we are interested
  # in is present
  function_as_text <- as.character(body(rhs_func)[2])
  do_fun_1 <- grepl(pattern = "rhs <- 0", x = function_as_text)
  do_fun_2 <- grepl(pattern = "rhs <- 1", x = function_as_text)
  do_fun_3 <- grepl(pattern = "rhs <- 2", x = function_as_text)

  if (do_fun_1)
  {
    lx <- (length(initprobs) - 1)/2
    parsvec <- c(DAISIE_loglik_rhs_precomp(pars,lx))
    y <- DAISIE_ode_cs(
      initprobs,
      tvec,
      parsvec,
      atol,
      rtol,
      method,
      runmod = "daisie_runmod"
    )
    #y <- deSolve::ode(
    #    y = initprobs,
    #    times = tvec,
    #    func = DAISIE_loglik_rhs1,
    #    parms = parsvec,
    #    rtol = rtol,
    #    atol = atol,
    #    method = method
    #  )[2, -1]
  } else if (do_fun_2)
  {
    lx <- (length(initprobs))/4
    parsvec <- c(DAISIE_loglik_rhs_precomp(pars,lx))
    y <- DAISIE_ode_cs(initprobs,
                       tvec,
                       parsvec,
                       atol,
                       rtol,
                       method,
                       runmod = "daisie_runmod1")

  } else if (do_fun_3)
  {
    lx <- (length(initprobs))/3
    parsvec <- c(DAISIE_loglik_rhs_precomp(pars,lx))
    y <- DAISIE_ode_cs(initprobs,
                       tvec,
                       parsvec,
                       atol,
                       rtol,
                       method,
                       runmod = "daisie_runmod2")
    #y <- deSolve::ode(
    #  y = initprobs,
    #  times = tvec,
    #  func = DAISIE_loglik_rhs2,
    #  parms = parsvec,
    #  rtol = rtol,
    #  atol = atol,
    #  method = method
    #)[2, -1]
  } else
  {
    stop(
      "The integrand function is written incorrectly. ",
      "Value of 'function_as_text':", function_as_text
    )
  }
  return(y)
}

# replacement for DAISIE_ode_FORTRAN
# returns striped deSolve result
DAISIE_ode_cs <- function(
    initprobs,
    tvec,
    parsvec,
    atol,
    rtol,
    methode,
    runmod = "daisie_runmod") {
  N <- length(initprobs)
  kk <- parsvec[length(parsvec)]
  if (runmod == "daisie_runmod") {
    lx <- (N - 1) / 2
    rhs_func <- DAISIE_loglik_rhs
  } else if (runmod == "daisie_runmod1")
  {
    lx <- N / 4
    rhs_func <- DAISIE_loglik_rhs1
  } else if (runmod == "daisie_runmod2") {
    lx <- N / 3
    rhs_func <- DAISIE_loglik_rhs2
  }
  if (startsWith(methode, "odeint")) {
    probs <- .Call("daisie_odeint_cs", runmod, initprobs, tvec, lx, kk, parsvec[-length(parsvec)], methode, atol, rtol)
  } else if (startsWith(methode, "deSolve_R::")) {
    methode <- substring(methode,12)
    y <- deSolve::ode(y = initprobs,
                      times = tvec,
                      func = rhs_func,
                      parms = parsvec,
                      atol = atol,
                      rtol = rtol,
                      method = methode)[,1:(N + 1)]
    probs <- y[-1,-1]
  } else {
    y <- deSolve::ode(y = initprobs, parms = c(lx + 0.,kk + 0.), rpar = parsvec[-length(parsvec)],
                      times = tvec, func = runmod, initfunc = "daisie_initmod",
                      ynames = c("SV"), dimens = N + 2, nout = 1, outnames = c("Sum"),
                      dllname = "DAISIE",atol = atol, rtol = rtol, method = methode)[,1:(N + 1)]
    probs <- y[-1,-1]  # strip 1st row and 1st column
  }
  return(probs)
}

# left in to cover the case this function is called from outside DAISIE_loglik_CS.R
DAISIE_ode_FORTRAN <- function(
    initprobs,
    tvec,
    parsvec,
    atol,
    rtol,
    methode,
    runmod = "daisie_runmod") {
  N <- length(initprobs)
  kk <- parsvec[length(parsvec)]
  if (runmod == "daisie_runmod") {
    lx <- (N - 1) / 2
  } else if (runmod == "daisie_runmod1") {
    lx <- N / 4
  } else if (runmod == "daisie_runmod2") {
    lx <- N / 3
  }
  probs <- deSolve::ode(y = initprobs, parms = c(lx + 0.,kk + 0.), rpar = parsvec[-length(parsvec)],
                        times = tvec, func = runmod, initfunc = "daisie_initmod",
                        ynames = c("SV"), dimens = N + 2, nout = 1, outnames = c("Sum"),
                        dllname = "DAISIE",atol = atol, rtol = rtol, method = methode)[,1:(N + 1)]
  return(probs)
}

logcondprob <- function(numcolmin, numimm, logp0, fac = 2) {
  if(numcolmin > sum(numimm)) {
    stop('The minimum number of colonizations cannot be smaller than the number of immigrants.')
  }
  maxi <- min(sum(numimm),fac * numcolmin)
  logcond <- 0
  if(numcolmin >= 1) {
    if(numcolmin == 1 && length(logp0) == 2) {
      message('With two types, conditioning on at least one colonization
              implies at least two colonizations. Therefore, the minimum
              number of colonizations is changed to 2.\n')
      numcolmin <- 2
    }
    lognotp0 <- rep(NA,length(logp0))
    for(i in 1:length(logp0))
    {
      if(exp(logp0[i]) == 1 & logp0[i] < 0)
      {
        lognotp0[i] <- log(-logp0[i])
      } else
      {
        lognotp0[i] <- log1p(-exp(logp0[i]))
      }
    }
    logpc <- matrix(0,nrow = maxi + 1,ncol = length(logp0))
    for(i in 0:maxi) {
      logpc[i + 1,] <- lgamma(numimm + 1) - lgamma(i + 1) - lgamma(numimm - i + 1) +
        (numimm - i) * logp0 + i * lognotp0
    }
    pc <- exp(logpc)
    if(length(logp0) == 2) {
      condprob <- pc[1,1] + pc[1,2] - pc[1,1] * pc[1,2]
      if(numcolmin > 2) {
        for(i in 2:(numcolmin - 1)) {
          condprob <- condprob + sum(pc[2:i,1] * pc[i:2,2])
        }
      }
      if(condprob >= 1) {
        logcond <- log(sum(pc[2:numcolmin,1] * pc[numcolmin:2,2]))
        message('A simple approximation of logcond must be made. Results may be unreliable.')
      } else {
        logcond <- log1p(-condprob)
      }
    } else {
      if(sum(pc) >= 1) {
        logcond <- log(sum(pc[(numcolmin + 1):(maxi + 1)]))
        message('An approximation of logcond must be made. Results may be unreliable.')
      } else {
        logcond <- log1p(-sum(pc[-((numcolmin + 1):(maxi + 1))]))
      }
    }
  }
  return(logcond)
}

#' @name DAISIE_logp0
#' @title Computes the log probability of no species present under the DAISIE
#' model with clade-specific diversity-dependence
#' @description Computes the log probability of no species present under the DAISIE
#' model with clade-specific diversity-dependence. The output is a log value.
#' @inheritParams default_params_doc
#' @param pars2 Contains the model settings \cr \cr
#' \code{pars2[1]} corresponds to lx = length of ODE variable x \cr
#' \code{pars2[2]} corresponds to ddmodel = diversity-dependent model, model of diversity-dependence, which can be one
#' of\cr \cr
#' ddmodel = 0 : no diversity dependence \cr
#' ddmodel = 1 : linear dependence in speciation rate \cr
#' ddmodel = 11: linear dependence in speciation rate and in immigration rate \cr
#' ddmodel = 2 : exponential dependence in speciation rate\cr
#' ddmodel = 21: exponential dependence in speciation rate and in immigration rate\cr\cr
#' @param island_age the island age \cr
#' @return The logarithm of the probability
#' @author Rampal S. Etienne & Bart Haegeman
#' @keywords internal
#' @export DAISIE_logp0
DAISIE_logp0 <- function(pars1,
                         pars2,
                         island_age,
                         methode = "odeint::runge_kutta_fehlberg78",
                         abstolint = 1E-16,
                         reltolint = 1E-10) {
  logp0 <- DAISIE_loglik_CS_choice(
    pars1 = pars1,
    pars2 = pars2,
    brts = island_age,
    stac = 0,
    missnumspec = 0,
    methode = methode,
    CS_version = 1,
    abstolint = abstolint,
    reltolint = reltolint)
  if(logp0 >= 0 & pars1[2]/pars1[1] > 100)
  {
    logp0 <- approximate_logp0(gamma = pars1[4],
                               mu = pars1[2],
                               t = island_age)
  }
  return(logp0)
}
