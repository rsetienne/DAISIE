DAISIE_loglik_CS_shift <- function(pars1,
                                   pars2,
                                   brts,
                                   stac,
                                   missnumspec,
                                   methode = "lsodes",
                                   abstolint = 1E-16,
                                   reltolint = 1E-10,
                                   verbose,
                                   tshift) {
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
    pars2[4] <- 0
  }
  ddep <- pars2[2]
  K <- pars1[3]
  #  if (!is.na(pars2[5])) {
  #    K <- K * pars1[8]
  #  }
  if(length(pars1) == 6) {
    probability_of_init_presence <- pars1[6]
    pars1 <- pars1[-6]
  } else {
    probability_of_init_presence <- 0
  }
  brts <- -sort(abs(as.numeric(brts)),decreasing = TRUE)
  tshift <- -abs(tshift)
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
  #  lac <- pars1[1]
  #  if(lac == Inf & missnumspec == 0 & length(pars1) == 5) {
  #    if(verbose) warning('Infinite lambda detected')
  #    loglik <- DAISIE_loglik_high_lambda(pars1, -brts, stac)
  #  } else {
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
    probs <- DAISIE_integrate_shift(probs,times = brts[1:2],tshift = tshift,loglik_rhs = DAISIE_loglik_rhs,pars1 = pars1,k1 = k1,ddep = ddep,rtol = reltolint,atol = abstolint,method = methode)
    cp <- checkprobs2(lv = 2 * lx, loglik, probs, verbose); loglik <- cp[[1]]; probs = cp[[2]]
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
        # To allow this we introduce a third set of equations for the
        # probability that colonization might have happened before but
        # recolonization has not taken place yet (Q_M,n).
        epss <- 1.01E-5 #We're taking the risk
        if (abs(brts[2] - brts[1]) >= epss) {
          probs[(2 * lx + 1):(4 * lx)] <- probs[1:(2 * lx)]
          probs[1:(2 * lx)] <- 0
        } else {
          probs[(2 * lx + 1):(4 * lx)] <- 0
        }

        probs <- DAISIE_integrate_shift(probs,times = brts[2:3],tshift = tshift,loglik_rhs = DAISIE_loglik_rhs1,pars1 = pars1,k1 = k1,ddep = ddep,rtol = reltolint,atol = abstolint,method = methode)
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
            probs <- DAISIE_integrate_shift(probs,times = brts[3:4],tshift = tshift,loglik_rhs = DAISIE_loglik_rhs2,pars1 = pars1,k1 = k1,ddep = ddep,rtol = reltolint,atol = abstolint,method = methode)
            cp <- checkprobs2(lx, loglik, probs, verbose); loglik <- cp[[1]]; probs = cp[[2]]
            loglik <- loglik + log(probs[(stac == 8) * (2 * lx + 1) + (stac == 9) + missnumspec])
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
        gamvec <- divdepvec(
          lac_or_gam = "gam",
          pars1 = ifelse(rep(t,5) <= rep(tshift,5), pars1[1:5], pars1[6:10]),
          t = t,
          lx = lx,
          k1 = k1,
          ddep = ddep * (ddep == 11 | ddep == 21)
        )
        probs[(2 * lx + 1):(3 * lx)] = gamvec[1:lx] * probs[1:lx] +
          gamvec[2:(lx + 1)] * probs[(lx + 1):(2 * lx)]
        probs[1:(2 * lx)] = 0
        k1 <- 1
        probs <- DAISIE_integrate_shift(probs,times = brts[2:3],tshift = tshift,loglik_rhs = DAISIE_loglik_rhs2,pars1 = pars1,k1 = k1,ddep = ddep,rtol = reltolint,atol = abstolint,method = methode)
        cp <- checkprobs2(lx,loglik,probs, verbose); loglik = cp[[1]]; probs = cp[[2]]
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
            pars1 = ifelse(rep(t,5) <= rep(tshift,5), pars1[1:5], pars1[6:10]),
            t = t,
            lx = lx + stac %in% c(6,7),
            k1 = k1,
            ddep = ddep
          )
          if(stac %in% c(2,3))
          {
            probs[1:lx] <- lacvec[1:lx] * (probs[1:lx] + probs[(2 * lx + 1):(3 * lx)])
            probs[(lx + 1):(2 * lx)] <- lacvec[2:(lx + 1)] * probs[(lx + 1):(2 * lx)]
          } else { # stac in c(6,7)
            probs[1:(lx - 1)] <- lacvec[2:lx] *
              ((1:(lx - 1)) * probs[2:lx] + probs[(2 * lx + 1):(3 * lx - 1)])
            probs[(lx + 1):(2 * lx - 1)] <- lacvec[3:(lx + 1)] * (1:(lx - 1)) *
              probs[(lx + 2):(2 * lx)]
            probs[lx] <- 0
            probs[2 * lx] <- 0
          }
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
            probs <- DAISIE_integrate_shift(probs,times = brts[k:(k+1)],tshift = tshift,loglik_rhs = DAISIE_loglik_rhs,pars1 = pars1,k1 = k1,ddep = ddep,rtol = reltolint,atol = abstolint,method = methode)
            cp <- checkprobs2(lx, loglik, probs, verbose); loglik <- cp[[1]]; probs = cp[[2]]
            if(k < S1)
            {
              # speciation event
              t <- brts[k + 1]
              lacvec <- divdepvec(
                lac_or_gam = "lac",
                pars1 = ifelse(rep(t,5) <= rep(tshift,5), pars1[1:5], pars1[6:10]),
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
  #}

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
    stop('stopping')
  }
  loglik <- as.numeric(loglik)
  #testit::assert(is.numeric(loglik))
  return(loglik)
}

DAISIE_integrate_shift <- function(probs,times,tshift,loglik_rhs,pars1,k1,ddep,rtol,atol,method) {
  if(tshift > times[1] && tshift >= times[2]) {
    probs <- DAISIE_integrate(probs,c(times[1:2]),loglik_rhs,c(pars1[1:5],k1,ddep),rtol = rtol,atol = atol,method = method)
  } else
    if(tshift > times[1] && tshift < times[2]) {
      probs <- DAISIE_integrate(probs,c(times[1],tshift),loglik_rhs,c(pars1[1:5],k1,ddep),rtol = rtol,atol = atol,method = method)
      probs <- DAISIE_integrate(probs,c(tshift,times[2]),loglik_rhs,c(pars1[6:10],k1,ddep),rtol = rtol,atol = atol,method = method)
    } else
      if(tshift <= times[1] && tshift < times[2]) {
        probs <- DAISIE_integrate(probs,c(times[1:2]),loglik_rhs,c(pars1[6:10],k1,ddep),rtol = rtol,atol = atol,method = method)
      }
  return(probs)
}
