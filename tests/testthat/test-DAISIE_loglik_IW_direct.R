# Tests of rhs = "direct" of DAISIE_loglik_IW(): the right-hand side and the event step of the island-wide likelihood that are
# computed on the fly (src/DAISIE_IW_direct.cpp), against the version with tables (rhs = "tables", the default).
#
# (a) the right-hand side, on random vectors with many exact zeros, for many parameter sets and small state spaces
# (b) the event step for every kind of event, against a copy of the code of DAISIE_loglik_IW() before rhs = "direct" existed
# (c) the loglik on several data sets
# (d) the maximum likelihood estimates on a small example
# (e) the default (rhs = "tables") has not changed
# and the numbering of the states and the error for state spaces with more than 2^31 - 1 states.  The last part of the file is
# about rhs = "compact" (see there).

# ---- helpers ---------------------------------------------------------------------------------------------------------

# the arguments of the tables version for one ODE segment (as in the loop of DAISIE_loglik_IW())
iw_test_tables_pars <- function(pars1, k, lxm, lxe, s) {
  sysdim <- 2 ^ s
  if (sysdim > 1) {
    l0ki <- create_l0ki(dec2binmat(s), lxm, lxe, sysdim)
  } else {
    l0ki <- list(l0 = 0, ki = NULL)
  }
  nndd <- nndivdep(lxm = lxm, lxe = lxe, sysdim = sysdim, Kprime = pars1[3], k = k, M = pars1[6], l0 = l0ki$l0)
  parslist <- list(pars = pars1, k = k, ddep = 11, dime = list(lxm = lxm, lxe = lxe, sysdim = sysdim), l0ki = l0ki, nndd = nndd)
  DAISIE_IW_pars(parslist)
}

# the maximum of the absolute differences relative to the largest absolute value
iw_test_reldiff <- function(a, b) {
  max(abs(a - b)) / max(abs(a), .Machine$double.xmin)
}

# runs expr with the given number of threads of the solver and puts the old number back
iw_test_with_threads <- function(num_threads, expr) {
  old <- DAISIE_IW_num_threads()
  on.exit(DAISIE_IW_num_threads(old), add = TRUE)
  DAISIE_IW_num_threads(num_threads)
  force(expr)
}

# a random state vector, a third of which is exactly 0
iw_test_random_state <- function(n) {
  x <- stats::runif(n)
  x[sample.int(n, floor(n / 3))] <- 0
  x
}

# data sets: the island age, the number of mainland species that are not on the island, and per colonist the branching times
# (stac 4: non-endemic, 2: endemic).  D1 has all kinds of events, D2 has five clades that are open at the same time (32
# configurations), D3 has radiations of which the events are interleaved with colonisations.
iw_test_data <- function(island_age, not_present, clades) {
  datalist <- list(list(island_age = island_age, not_present = not_present))
  for (i in seq_along(clades)) {
    datalist[[i + 1]] <- list(colonist_name = paste0("clade", i),
                              branching_times = c(island_age, clades[[i]]$bt),
                              stac = clades[[i]]$stac,
                              missing_species = 0)
  }
  datalist
}
iw_test_D1 <- function() {
  iw_test_data(5, 40, list(list(stac = 4, bt = 4.0), list(stac = 2, bt = 3.0), list(stac = 2, bt = c(2.5, 1.5, 0.7)),
                           list(stac = 4, bt = 0.5)))
}
iw_test_D2 <- function() {
  iw_test_data(5, 40, lapply(c(4.5, 3.5, 2.5, 1.5, 0.5), function(t) list(stac = 4, bt = t)))
}
iw_test_D3 <- function() {
  iw_test_data(6, 60, list(list(stac = 2, bt = c(5.5, 4.5, 2.0)), list(stac = 4, bt = 5.0), list(stac = 2, bt = c(4.0, 3.5, 1.0)),
                           list(stac = 2, bt = 3.0), list(stac = 4, bt = 0.7)))
}

iw_test_fehlberg <- "odeint::runge_kutta_fehlberg78"

# ---- the numbering of the states -------------------------------------------------------------------------------------

test_that("the states are numbered as in the array of the tables version, in R and in C++", {
  # 1-based position of (n1 = i, n2 = j, configuration a) in the array with dim c(lxm, lxe, 2^s)
  lxm <- 5
  lxe <- 4
  for (s in 0:3) {
    sysdim <- 2 ^ s
    arr <- array(seq_len(lxm * lxe * sysdim), dim = c(lxm, lxe, sysdim))
    g <- expand.grid(i = 0:(lxm - 1), j = 0:(lxe - 1), a = 0:(sysdim - 1))
    idx <- .Call("daisie_iw_state_index", g$i, g$j, g$a, lxm, lxe)
    expect_identical(idx, as.numeric(arr[cbind(g$i + 1, g$j + 1, g$a + 1)]))
    expect_identical(idx, 1 + g$i + lxm * (g$j + lxe * g$a))
  }
  # in 64 bit: states beyond 2^31 are not wrapped around
  lxm <- 2 ^ 16
  lxe <- 2 ^ 16 + 3
  i <- c(0, lxm - 1, 17, lxm - 1)
  j <- c(0, lxe - 1, 5, lxe - 1)
  a <- c(0, 31, 20, 1)
  idx <- .Call("daisie_iw_state_index", i, j, a, lxm, lxe)
  expect_identical(idx, 1 + i + lxm * (j + lxe * a))
  expect_true(idx[2] > 2 ^ 31)
})

test_that("a configuration is a bit pattern of the open clades: l0 and ki of the tables version follow it", {
  for (s in 1:4) {
    sysdim <- 2 ^ s
    l0ki <- create_l0ki(dec2binmat(s), 2, 3, sysdim)
    # l0 = number of open clades that are not occupied = s - number of bits that are set
    popcount <- vapply(0:(sysdim - 1), function(a) sum(bitwAnd(a, 2 ^ (0:(s - 1))) > 0), numeric(1))
    expect_equal(l0ki$l0[1, 1, ], s - popcount)
    # ki[a, c] = 1 if and only if c is a with one bit that is set cleared
    ki <- matrix(0, sysdim, sysdim)
    for (a in 1:(sysdim - 1)) {
      for (b in 0:(s - 1)) {
        if (bitwAnd(a, 2 ^ b) > 0) {
          ki[1 + a, 1 + bitwXor(a, 2 ^ b)] <- 1
        }
      }
    }
    expect_identical(l0ki$ki, ki)
  }
})

# ---- (a) the right-hand side -----------------------------------------------------------------------------------------

test_that("(a) the right-hand side with tables (C++) and direct agree on random vectors with many zeros", {
  set.seed(20260101)
  # lac, mu, K', gam, laa, M ; K' small (the diversity-dependence is 0 for large state), M small (pmax(0, M - n) is active),
  # K' = Inf, and zero rates
  parsets <- list(
    c(0.4, 0.2, 12, 0.05, 0.3, 20),
    c(0.5, 0.1, 5, 0.1, 0.2, 6),
    c(1.3, 0.7, 30.5, 0.9, 0.05, 8),
    c(0.4, 0.2, Inf, 0.05, 0.3, 25),
    c(0, 0.2, 9, 0.05, 0.3, 15),
    c(0.4, 0, 9, 0.05, 0.3, 15),
    c(0.4, 0.2, 9, 0, 0.3, 15),
    c(0.4, 0.2, 9, 0.05, 0, 15),
    c(0, 0, 9, 0, 0, 15),
    c(2, 3, 1, 4, 5, 3)
  )
  settings <- list(c(lxm = 7, lxe = 9), c(lxm = 2, lxe = 2), c(lxm = 3, lxe = 8), c(lxm = 8, lxe = 3), c(lxm = 2, lxe = 5))
  maxdiff <- 0
  for (pars1 in parsets) {
    for (setting in settings) {
      lxm <- setting[["lxm"]]
      lxe <- setting[["lxe"]]
      for (s in 0:4) {
        for (k in c(0, 1, 4, 6)) {
          if (s > k) next                          # s open clades need at least s colonisations
          n <- lxm * lxe * 2 ^ s
          x <- iw_test_random_state(n)
          tables <- .Call("daisie_iw_rhs_eval", x, iw_test_tables_pars(pars1, k, lxm, lxe, s))
          direct <- .Call("daisie_iw_rhs_eval_direct", x, DAISIE_IW_direct_pars(pars1, pars1[3], k, lxm, lxe, s))
          expect_length(direct, n)
          d <- iw_test_reldiff(tables, direct)
          maxdiff <- max(maxdiff, d)
          if (!(d < 1E-12)) {
            fail(paste0("rhs differs: relative difference ", d, " for pars1 = ", paste(pars1, collapse = ","),
                        ", lxm = ", lxm, ", lxe = ", lxe, ", s = ", s, ", k = ", k))
          }
        }
      }
    }
  }
  expect_lt(maxdiff, 1E-12)
})

test_that("(a) the right-hand side direct agrees with the R function used for the non-C++ solvers", {
  skip_if_not_installed("tensor")
  set.seed(20260102)
  pars1 <- c(0.4, 0.2, 9, 0.05, 0.3, 12)
  for (s in 0:3) {
    for (k in c(s, s + 2)) {
      lxm <- 6
      lxe <- 5
      x <- iw_test_random_state(lxm * lxe * 2 ^ s)
      cp <- iw_test_tables_pars(pars1, k, lxm, lxe, s)
      r <- DAISIE_loglik_rhs_IW(0, x, cp)[[1]]
      direct <- .Call("daisie_iw_rhs_eval_direct", x, DAISIE_IW_direct_pars(pars1, pars1[3], k, lxm, lxe, s))
      expect_lt(iw_test_reldiff(as.numeric(r), direct), 1E-12)
    }
  }
})

test_that("(a) the right-hand side direct does not depend on the number of threads", {
  set.seed(20260103)
  pars1 <- c(0.4, 0.2, 12, 0.05, 0.3, 20)
  lxm <- 9
  lxe <- 11
  s <- 4
  x <- iw_test_random_state(lxm * lxe * 2 ^ s)
  p <- DAISIE_IW_direct_pars(pars1, pars1[3], 5, lxm, lxe, s)
  one <- iw_test_with_threads(1, .Call("daisie_iw_rhs_eval_direct", x, p))
  two <- iw_test_with_threads(2, .Call("daisie_iw_rhs_eval_direct", x, p))
  expect_identical(one, two)
})

test_that("(a) NaN in the state: direct gives NaN where there is a flow from it, the tables version also elsewhere", {
  # the tables version multiplies with the dense matrix ki, so 0 * NaN = NaN reaches all configurations of a pair (n1, n2);
  # direct only adds what flows.  The state itself stays NaN in both (the term with c7), so checkprobs2() sees it in both.
  pars1 <- c(0.4, 0.2, 9, 0.05, 0.3, 12)
  lxm <- 4
  lxe <- 5
  s <- 2
  set.seed(20260104)
  x <- iw_test_random_state(lxm * lxe * 2 ^ s)
  x[7] <- NaN
  tables <- .Call("daisie_iw_rhs_eval", x, iw_test_tables_pars(pars1, 3, lxm, lxe, s))
  direct <- .Call("daisie_iw_rhs_eval_direct", x, DAISIE_IW_direct_pars(pars1, pars1[3], 3, lxm, lxe, s))
  expect_true(is.nan(direct[7]))
  expect_true(is.nan(tables[7]))
  expect_true(all(is.nan(tables)[is.nan(direct)]))
  # elsewhere they agree
  ok <- !is.nan(tables) & !is.nan(direct)
  expect_lt(iw_test_reldiff(tables[ok], direct[ok]), 1E-12)
})

# ---- (b) the event step ----------------------------------------------------------------------------------------------

# the event step as it is in the loop of DAISIE_loglik_IW() (a verbatim copy of what that function does between two ODE
# segments, before rhs = "direct" existed).  type 1: colonisation, 2: speciation, 3: first speciation in a clade, of which the
# bit tocollapse (1-based position in expandvec) is removed.
iw_test_event_reference <- function(probs, type, pars1, k, lxm, lxe, s, tocollapse = 1) {
  lac <- pars1[1]
  gam <- pars1[4]
  M <- pars1[6]
  Kprime <- pars1[3]
  sysdim <- 2 ^ s
  totdim <- lxm * lxe * sysdim
  if (sysdim > 1) {
    l0ki <- create_l0ki(dec2binmat(s), lxm, lxe, sysdim)
  } else {
    l0ki <- list(l0 = 0, ki = NULL)
  }
  nndd <- nndivdep(lxm = lxm, lxe = lxe, sysdim = sysdim, Kprime = Kprime, k = k, M = M, l0 = l0ki$l0)
  dim(probs) <- c(lxm, lxe, sysdim)
  divdepfac <- nndd$divdepfac
  divdepfacplus1 <- nndd$divdepfacplus1
  mfac <- nndd$mfac
  oneminmfac <- nndd$oneminmfac
  if (type == 1) {
    probs2 <- array(0, dim = dim(probs))
    probs2[1:(lxm - 1), , ] <- probs[2:lxm, , ]
    probs <- gam * divdepfac * oneminmfac * probs[, , 1:sysdim] +
      gam * divdepfacplus1 * mfac * probs2[, , 1:sysdim]
    probs <- c(probs, rep(0, totdim))
  } else {
    probs <- lac * divdepfac * probs[, , 1:sysdim]
    if (type == 3) {
      sr <- selectrows(sysdim, 2 ^ (tocollapse - 1))
      probs <- probs[, , sr[, 1]] + probs[, , sr[, 2]]
      sysdim <- sysdim / 2
      dim(probs) <- c(lxm, lxe, sysdim)
    }
  }
  as.numeric(probs)
}

test_that("(b) the event step direct agrees with the event step of the tables version, for every kind of event", {
  set.seed(20260105)
  parsets <- list(
    c(0.4, 0.2, 12, 0.05, 0.3, 20),
    c(0.5, 0.1, 5, 0.1, 0.2, 9),                   # M small: pmax is active; the mainland pool is nearly exhausted
    c(1.3, 0.7, 30.5, 0.9, 0.05, 8),
    c(0.4, 0.2, Inf, 0.05, 0.3, 25),
    c(0, 0.2, 9, 0.05, 0.3, 15),
    c(0.4, 0.2, 9, 0, 0.3, 15)
  )
  settings <- list(c(lxm = 7, lxe = 9), c(lxm = 2, lxe = 2), c(lxm = 6, lxe = 3), c(lxm = 3, lxe = 6))
  maxdiff <- 0
  nchecks <- 0
  for (pars1 in parsets) {
    M <- pars1[6]
    for (setting in settings) {
      lxm <- setting[["lxm"]]
      lxe <- setting[["lxe"]]
      for (s in 0:4) {
        for (k in unique(c(s, s + 1, s + 3))) {
          if (k >= M - 1) next                    # M - l0 must be positive (as it is in the likelihood)
          n <- lxm * lxe * 2 ^ s
          y <- iw_test_random_state(n)
          p <- DAISIE_IW_direct_pars(pars1, pars1[3], k, lxm, lxe, s)
          # colonisation
          ref <- iw_test_event_reference(y, 1, pars1, k, lxm, lxe, s)
          dir <- .Call("daisie_iw_event", y, 1L, p, 0L)
          expect_length(dir, 2 * n)
          maxdiff <- max(maxdiff, iw_test_reldiff(ref, dir))
          # speciation
          ref <- iw_test_event_reference(y, 2, pars1, k, lxm, lxe, s)
          dir <- .Call("daisie_iw_event", y, 2L, p, 0L)
          expect_length(dir, n)
          maxdiff <- max(maxdiff, iw_test_reldiff(ref, dir))
          nchecks <- nchecks + 2
          # first speciation in a clade: for each of the open clades
          if (s >= 1) {
            for (tocollapse in 1:s) {
              ref <- iw_test_event_reference(y, 3, pars1, k, lxm, lxe, s, tocollapse)
              dir <- .Call("daisie_iw_event", y, 3L, p, as.integer(tocollapse - 1))
              expect_length(dir, n / 2)
              maxdiff <- max(maxdiff, iw_test_reldiff(ref, dir))
              nchecks <- nchecks + 1
            }
          }
        }
      }
    }
  }
  expect_gt(nchecks, 500)
  expect_lt(maxdiff, 1E-13)
})

test_that("(b) the event step direct: wrong input is an error and not a crash", {
  pars1 <- c(0.4, 0.2, 9, 0.05, 0.3, 12)
  p <- DAISIE_IW_direct_pars(pars1, pars1[3], 2, 4, 5, 2)
  expect_error(.Call("daisie_iw_event", rep(0.1, 7), 1L, p, 0L))               # wrong length of the state vector
  expect_error(.Call("daisie_iw_rhs_eval_direct", rep(0.1, 7), p))
})

test_that("(b) the normalisation after an event is checkprobs2()", {
  set.seed(20260106)
  for (verbose in c(FALSE, TRUE)) {
    # usual vectors, vectors with negative values, with exact zeros
    for (trial in 1:20) {
      y <- stats::rnorm(50, mean = 0.1, sd = 0.2) * 10 ^ stats::runif(1, -8, 2)
      y[sample.int(50, 10)] <- 0
      old <- checkprobs2(NA, 1.5, y, verbose)
      new <- DAISIE_IW_check_probs(y, 1.5, verbose)
      expect_equal(new[[1]], old[[1]], tolerance = 1E-14)
      expect_equal(new[[2]], old[[2]], tolerance = 1E-14)
    }
    # all zero or negative: loglik is -Inf
    for (y in list(rep(0, 10), c(-1, -2, 0, -3))) {
      expect_identical(DAISIE_IW_check_probs(y, 1.5, FALSE)[[1]], checkprobs2(NA, 1.5, y, FALSE)[[1]])
      expect_identical(DAISIE_IW_check_probs(y, 1.5, FALSE)[[1]], -Inf)
    }
    # NaN: loglik is -Inf
    y <- c(0.1, NaN, 0.3)
    expect_identical(DAISIE_IW_check_probs(y, 1.5, FALSE)[[1]], checkprobs2(NA, 1.5, y, FALSE)[[1]])
    expect_identical(DAISIE_IW_check_probs(y, 1.5, FALSE)[[1]], -Inf)
  }
  expect_message(DAISIE_IW_check_probs(c(0, 0), 0, TRUE), "Numerical issues")
  expect_silent(DAISIE_IW_check_probs(c(0, 0), 0, FALSE))
})

# ---- (c) the loglik --------------------------------------------------------------------------------------------------

test_that("(c) loglik direct = loglik tables on hand-made data sets with all kinds of events", {
  pars1 <- c(0.3, 0.2, 12, 0.02, 0.3)
  for (cond in 0:1) {
    pars2 <- c(12, 11, cond, 0)
    for (datalist in list(iw_test_D1(), iw_test_D2(), iw_test_D3())) {
      tables <- DAISIE_loglik_IW(pars1, pars2, datalist, methode = iw_test_fehlberg)
      direct <- DAISIE_loglik_IW(pars1, pars2, datalist, methode = iw_test_fehlberg, rhs = "direct")
      expect_true(is.finite(tables))
      expect_equal(direct, tables, tolerance = 1E-8)
    }
  }
})

test_that("(c) loglik direct = loglik tables for K' = Inf, another M, and other parameters", {
  datalist <- iw_test_D3()
  for (pars1 in list(c(0.3, 0.2, Inf, 0.02, 0.3), c(0.8, 0.5, 20, 0.1, 0.05), c(0.2, 0.1, 15, 0.005, 1))) {
    tables <- DAISIE_loglik_IW(pars1, c(10, 11, 0, 0), datalist, methode = iw_test_fehlberg)
    direct <- DAISIE_loglik_IW(pars1, c(10, 11, 0, 0), datalist, methode = iw_test_fehlberg, rhs = "direct")
    expect_true(is.finite(tables))
    expect_equal(direct, tables, tolerance = 1E-8)
  }
})

test_that("(c) loglik direct = loglik tables with the other solvers of odeint", {
  skip_on_cran()
  pars1 <- c(0.3, 0.2, 12, 0.02, 0.3)
  datalist <- iw_test_D1()
  for (methode in c("odeint::runge_kutta_cash_karp54", "odeint::runge_kutta_dopri5", "odeint::bulirsch_stoer",
                    "odeint::rosenbrock4")) {
    tables <- DAISIE_loglik_IW(pars1, c(8, 11, 0, 0), datalist, methode = methode, abstolint = 1E-12, reltolint = 1E-10)
    direct <- DAISIE_loglik_IW(pars1, c(8, 11, 0, 0), datalist, methode = methode, abstolint = 1E-12, reltolint = 1E-10,
                               rhs = "direct")
    expect_equal(direct, tables, tolerance = 1E-8, label = methode)
  }
})

test_that("(c) loglik direct is the same with one and with two threads", {
  skip_on_cran()
  pars1 <- c(0.3, 0.2, 12, 0.02, 0.3)
  one <- iw_test_with_threads(1, DAISIE_loglik_IW(pars1, c(12, 11, 0, 0), iw_test_D2(), methode = iw_test_fehlberg, rhs = "direct"))
  two <- iw_test_with_threads(2, DAISIE_loglik_IW(pars1, c(12, 11, 0, 0), iw_test_D2(), methode = iw_test_fehlberg, rhs = "direct"))
  expect_identical(one, two)
})

test_that("(c) loglik direct = loglik tables for the frogs data", {
  skip_on_cran()
  utils::data(frogs_datalist, package = "DAISIE")
  pars1 <- c(0.2, 0.1, 1000.1, 0.001, 0.3)
  pars2 <- c(15, 11, 0, 0)
  tables <- DAISIE_loglik_IW(pars1, pars2, frogs_datalist, methode = iw_test_fehlberg)
  direct <- DAISIE_loglik_IW(pars1, pars2, frogs_datalist, methode = iw_test_fehlberg, rhs = "direct")
  expect_equal(direct, tables, tolerance = 1E-8)
  expect_equal(direct, -211.0636352905, tolerance = 1E-8)
})

test_that("(c) loglik direct = loglik tables for a data set with a recolonisation", {
  skip_on_cran()
  # data set 897 from the CS simulations, as in test-integration_DAISIE.R; a recolonisation is a colonisation of a clade that
  # has already been colonised, so for the likelihood it is an event of kind 1 (the branch for 'a recolonisation without
  # speciation' of DAISIE_loglik_IW() cannot be reached, and rhs = "direct" does not have it)
  utils::data(frogs_sim_datalist, package = "DAISIE")
  pars1 <- c(0.4012, 0.17, 131.9, 0.00035, 0)
  pars2 <- c(20, 11, 1, 0)
  datalist2 <- frogs_sim_datalist
  datalist2[[2]]$all_colonisations[[2]]$event_times <- c(datalist2[[2]]$all_colonisations[[2]]$event_times, 5)
  for (datalist in list(frogs_sim_datalist, datalist2)) {
    tables <- DAISIE_loglik_IW(pars1, pars2, datalist, methode = iw_test_fehlberg, abstolint = 1E-12, reltolint = 1E-10)
    direct <- DAISIE_loglik_IW(pars1, pars2, datalist, methode = iw_test_fehlberg, abstolint = 1E-12, reltolint = 1E-10,
                               rhs = "direct")
    expect_true(is.finite(tables))
    expect_equal(direct, tables, tolerance = 1E-8)
  }
})

test_that("(c) loglik direct = loglik tables for simulated islands", {
  skip_on_cran()
  set.seed(42)
  islands <- DAISIE_sim(time = 3, M = 100, pars = c(0.5, 0.3, 8, 0.01, 0.4), replicates = 6, divdepmodel = "IW", cond = 1,
                        verbose = FALSE, plot_sims = FALSE, use_rcpp = TRUE)
  pars1 <- c(0.4, 0.25, 7, 0.012, 0.3)
  pars2 <- c(12, 11, 1, 0)
  for (i in seq_along(islands)) {
    isl <- islands[[i]]
    tables <- DAISIE_loglik_IW(pars1, pars2, isl, methode = iw_test_fehlberg)
    direct <- DAISIE_loglik_IW(pars1, pars2, isl, methode = iw_test_fehlberg, rhs = "direct")
    expect_equal(direct, tables, tolerance = 1E-8, label = paste("simulated island", i))
  }
})

test_that("(c) non-finite parameters give the same result (NaN) with both", {
  datalist <- iw_test_D1()
  for (pars1 in list(c(0.3, 0.2, 12, Inf, 0.3), c(Inf, 0.2, 12, 0.02, 0.3), c(0.3, 0.2, 12, 0.02, Inf))) {
    tables <- DAISIE_loglik_IW(pars1, c(10, 11, 0, 0), datalist, methode = iw_test_fehlberg)
    direct <- DAISIE_loglik_IW(pars1, c(10, 11, 0, 0), datalist, methode = iw_test_fehlberg, rhs = "direct")
    expect_identical(direct, tables)
  }
})

test_that("(c) rhs = direct with a likelihood that is 0: -Inf as in the tables version", {
  # K' smaller than the number of species: the likelihood is -Inf before any ODE is solved
  datalist <- iw_test_D2()
  pars1 <- c(0.3, 0.2, 2, 0.02, 0.3)
  expect_message(tables <- DAISIE_loglik_IW(pars1, c(12, 11, 0, 0), datalist, methode = iw_test_fehlberg))
  expect_message(direct <- DAISIE_loglik_IW(pars1, c(12, 11, 0, 0), datalist, methode = iw_test_fehlberg, rhs = "direct"))
  expect_identical(direct, tables)
  # negative parameter
  expect_warning(direct <- DAISIE_loglik_IW(c(0.3, -0.2, 12, 0.02, 0.3), c(12, 11, 0, 0), datalist, methode = iw_test_fehlberg,
                                            rhs = "direct"))
  expect_identical(direct, -Inf)
})

# ---- the choice of rhs and the error for impossible combinations -----------------------------------------------------

test_that("rhs = direct is refused for the solvers that are not C++", {
  pars1 <- c(0.3, 0.2, 12, 0.02, 0.3)
  for (methode in c("lsodes", "ode45", "lsoda")) {
    expect_error(DAISIE_loglik_IW(pars1, c(8, 11, 0, 0), iw_test_D1(), methode = methode, rhs = "direct"), "odeint")
  }
  expect_error(DAISIE_loglik_IW(pars1, c(8, 11, 0, 0), iw_test_D1(), methode = iw_test_fehlberg, rhs = "something"),
               "should be one of")
})

test_that("rhs = direct stops with a clear message if the state space has more than 2^31 - 1 states", {
  # the number of states is lxm x lxe x 2^s
  expect_silent(DAISIE_IW_check_dim(2 ^ 31 - 1, 1, 1))
  expect_error(DAISIE_IW_check_dim(2 ^ 31, 1, 1), "2\\^31 - 1")
  expect_error(DAISIE_IW_check_dim(1000, 1000, 2 ^ 12), "more than")           # 4.1e9 states
  expect_error(DAISIE_IW_check_dim(50000, 50000, 2), "pars2\\[1\\]")
  expect_silent(DAISIE_IW_check_dim(1000, 1000, 2 ^ 11))                         # 2.048e9 states
  expect_identical(DAISIE_IW_max_states(), 2147483647L)
  # the largest number of open clades is found from the events
  expect_identical(DAISIE_IW_max_open_clades(c(0, 1, 1, 2, 1, 2, 2), 7), 2)
  # through the likelihood: a large lx and five open clades; the check is before anything is allocated
  datalist <- iw_test_D2()
  datalist[[1]]$not_present <- 1E5
  expect_error(
    DAISIE_loglik_IW(c(0.3, 0.2, 1E6, 0.02, 0.3), c(50000, 11, 0, 0), datalist, methode = iw_test_fehlberg, rhs = "direct"),
    "2\\^31 - 1")
})

# ---- (d) maximum likelihood ------------------------------------------------------------------------------------------

test_that("(d) the ML estimates are the same with rhs = tables and rhs = direct", {
  skip_on_cran()
  datalist <- iw_test_data(5, 20, list(list(stac = 4, bt = 4.0), list(stac = 2, bt = c(3.0, 1.5)), list(stac = 4, bt = 0.5)))
  args <- list(datalist = datalist, initparsopt = c(0.5, 0.3), idparsopt = c(1, 2), parsfix = c(12, 0.05, 0.3),
               idparsfix = c(3, 4, 5), res = 8, ddmodel = 11, cond = 0, methode = iw_test_fehlberg, tolint = c(1E-12, 1E-10),
               tol = c(1E-3, 1E-4, 1E-6), maxiter = 200)
  tables <- suppressMessages(do.call(DAISIE_ML_IW, args))
  direct <- suppressMessages(do.call(DAISIE_ML_IW, c(args, list(rhs = "direct"))))
  expect_true(is.finite(tables$loglik))
  expect_equal(direct$loglik, tables$loglik, tolerance = 1E-8)
  expect_equal(unlist(direct[1:5]), unlist(tables[1:5]), tolerance = 1E-6)
  expect_equal(direct$conv, tables$conv)
  # rhs = "tables" is what the function does if the argument is not given
  default <- suppressMessages(do.call(DAISIE_ML_IW, args))
  expect_identical(default, tables)
})

# ---- (e) the default has not changed ---------------------------------------------------------------------------------

test_that("(e) the default is rhs = tables, and gives what the version before rhs existed gave", {
  expect_identical(eval(formals(DAISIE_loglik_IW)$rhs), "tables")
  expect_identical(eval(formals(DAISIE_ML_IW)$rhs), "tables")
  pars1 <- c(0.3, 0.2, 12, 0.02, 0.3)
  pars2 <- c(12, 11, 0, 0)
  d1 <- iw_test_D1()
  # the argument that is not given is the same as rhs = "tables"
  expect_identical(DAISIE_loglik_IW(pars1, pars2, d1, methode = iw_test_fehlberg),
                   DAISIE_loglik_IW(pars1, pars2, d1, methode = iw_test_fehlberg, rhs = "tables"))
  # values computed with the code before rhs was introduced (git 6e3c1ab of rsetienne/DAISIE, branch develop); R 4.3.3, Linux,
  # GCC, x86-64.  The tolerance allows other platforms to differ in the last digits.
  expect_equal(DAISIE_loglik_IW(pars1, pars2, d1), -28.011763791047176, tolerance = 1E-12)
  expect_equal(DAISIE_loglik_IW(pars1, c(12, 11, 1, 0), d1), -27.961783764079872, tolerance = 1E-12)
  expect_equal(DAISIE_loglik_IW(pars1, pars2, iw_test_D2(), methode = iw_test_fehlberg), -31.830508973137299, tolerance = 1E-12)
  expect_equal(DAISIE_loglik_IW(pars1, pars2, iw_test_D3(), methode = iw_test_fehlberg), -41.717724796772835, tolerance = 1E-12)
  expect_equal(DAISIE_loglik_IW(pars1, pars2, iw_test_D3(), methode = "ode45"), -41.717724737315294, tolerance = 1E-12)
  expect_equal(DAISIE_loglik_IW(c(0.3, 0.2, Inf, 0.02, 0.3), pars2, d1, methode = "odeint::runge_kutta_dopri5"),
               -28.198588763278369, tolerance = 1E-12)
  skip_on_cran()
  utils::data(frogs_datalist, package = "DAISIE")
  expect_equal(DAISIE_loglik_IW(c(0.2, 0.1, 1000.1, 0.001, 0.3), c(15, 11, 0, 0), frogs_datalist, methode = iw_test_fehlberg),
               -211.06363529050506, tolerance = 1E-12)
})

# ---- rhs = "compact" ---------------------------------------------------------------------------------------------------
#
# rhs = "compact" is rhs = "direct" that does not store the states with i + j > ceiling(K') - k (k = number of events so far),
# which have probability exactly 0 for a finite K'.  Per configuration of the open clades the stored states are numbered per j
# with i = 0..len[j]-1 (DAISIE_IW_layout() in R, daisie_iw_layout() in C++).  The tests here compare it with the box layout of
# rhs = "direct" (a) for the numbering, (b) for the right-hand side and the event step on vectors that are 0 outside the
# triangle, (c) for the loglik and the ML estimates, and test the errors.

# the compact positions (1-based) of states (i, j, a) in R from DAISIE_IW_layout(); NA for the states that are not stored
iw_test_compact_index <- function(i, j, a, lxm, lxe, tri) {
  lay <- DAISIE_IW_layout(lxm, lxe, tri)
  stored <- i + j <= tri
  idx <- rep(NA_real_, length(i))
  idx[stored] <- 1 + a[stored] * lay$P + lay$off[j[stored] + 1] + i[stored]
  idx
}

# the box states (rows of the data frame) of lxm x lxe x nconf, and for tri the logical vector of the stored ones
iw_test_box <- function(lxm, lxe, nconf) {
  expand.grid(i = 0:(lxm - 1), j = 0:(lxe - 1), a = 0:(nconf - 1))
}

test_that("compact: the layout of the stored states is the same in R and in C++, and is the triangle i + j <= tri", {
  lxm <- 5
  lxe <- 4
  g <- iw_test_box(lxm, lxe, 4)
  for (tri in c(Inf, 100, lxm + lxe - 2, lxm + lxe - 3, 5, 3, 1, 0)) {
    lay <- .Call("daisie_iw_layout", lxm, lxe, tri)
    ref <- DAISIE_IW_layout(lxm, lxe, tri)
    expect_equal(lay$len, ref$len)
    expect_equal(lay$off, ref$off)
    expect_equal(lay$P, ref$P)
    # P is the number of states of the box with i + j <= tri
    expect_equal(lay$P, sum(g$a == 0 & g$i + g$j <= tri))
    idx <- .Call("daisie_iw_state_index_tri", g$i, g$j, g$a, lxm, lxe, tri)
    expect_identical(is.na(idx), g$i + g$j > tri)
    expect_identical(idx, iw_test_compact_index(g$i, g$j, g$a, lxm, lxe, tri))
    # every stored state has its own place, and the places are 1..P * 4
    expect_identical(sort(idx[!is.na(idx)]), as.numeric(seq_len(lay$P * 4)))
    # no state is stored if tri is large: the numbering of rhs = "direct" and of the tables version
    if (tri >= lxm + lxe - 2) {
      expect_equal(lay$P, lxm * lxe)
      expect_identical(idx, .Call("daisie_iw_state_index", g$i, g$j, g$a, lxm, lxe))
    }
  }
  # in 64 bit: a position beyond 2^31 is not wrapped around
  lxm <- 2 ^ 16
  lxe <- 2 ^ 16
  tri <- 2 ^ 16 + 10
  lay <- .Call("daisie_iw_layout", lxm, lxe, tri)
  expect_gt(lay$P, 2 ^ 31)
  i <- c(0, 5, 100, 11)
  j <- c(0, 7, 20, 2 ^ 16 - 1)
  a <- c(0, 1, 3, 2)
  keep <- i + j <= tri
  idx <- .Call("daisie_iw_state_index_tri", i[keep], j[keep], a[keep], lxm, lxe, tri)
  expect_identical(idx, iw_test_compact_index(i[keep], j[keep], a[keep], lxm, lxe, tri))
  expect_gt(idx[3], 2 ^ 32)
})

test_that("compact: the bound tri follows from K' and the number of events", {
  expect_identical(DAISIE_IW_tri(12, 0, TRUE), 12)
  expect_identical(DAISIE_IW_tri(12.2, 3, TRUE), 10)
  expect_identical(DAISIE_IW_tri(Inf, 3, TRUE), Inf)
  expect_identical(DAISIE_IW_tri(12, 3, FALSE), Inf)
})

test_that("compact: the right-hand side is 0 outside the triangle, and inside it is the right-hand side of the box", {
  set.seed(20260201)
  parsets <- list(
    c(0.4, 0.2, 12, 0.05, 0.3, 20),
    c(0.5, 0.1, 5, 0.1, 0.2, 6),
    c(1.3, 0.7, 30.5, 0.9, 0.05, 8),
    c(0.4, 0.2, Inf, 0.05, 0.3, 25),
    c(0, 0.2, 9, 0.05, 0.3, 15),
    c(0.4, 0, 9, 0.05, 0.3, 15),
    c(0.4, 0.2, 9, 0, 0.3, 15),
    c(0.4, 0.2, 9, 0.05, 0, 15),
    c(0.4, 0.2, 9.5, 0.05, 0.3, 15)               # K' that is not a whole number: ceiling(K') - k
  )
  settings <- list(c(lxm = 7, lxe = 9), c(lxm = 2, lxe = 2), c(lxm = 3, lxe = 8), c(lxm = 8, lxe = 3), c(lxm = 2, lxe = 5))
  maxdiff <- 0
  nchecks <- 0
  ncompact <- 0
  for (pars1 in parsets) {
    for (setting in settings) {
      lxm <- setting[["lxm"]]
      lxe <- setting[["lxe"]]
      for (s in 0:4) {
        for (k in c(0, 1, 4, 6)) {
          if (s > k) next
          tri <- DAISIE_IW_tri(pars1[3], k, TRUE)
          if (tri < 0) next
          g <- iw_test_box(lxm, lxe, 2 ^ s)
          stored <- g$i + g$j <= tri
          x <- iw_test_random_state(nrow(g))
          x[!stored] <- 0
          box <- .Call("daisie_iw_rhs_eval_direct", x, DAISIE_IW_direct_pars(pars1, pars1[3], k, lxm, lxe, s))
          # the invariant behind the compact storage: the box right-hand side does not feed the states outside the triangle
          expect_true(all(box[!stored] == 0))
          pc <- DAISIE_IW_direct_pars(pars1, pars1[3], k, lxm, lxe, s, tri)
          idx <- .Call("daisie_iw_state_index_tri", g$i[stored], g$j[stored], g$a[stored], lxm, lxe, tri)
          xc <- numeric(DAISIE_IW_layout(lxm, lxe, tri)$P * 2 ^ s)
          xc[idx] <- x[stored]
          compact <- .Call("daisie_iw_rhs_eval_direct", xc, pc)
          expect_length(compact, length(xc))
          d <- iw_test_reldiff(compact[idx], box[stored])
          maxdiff <- max(maxdiff, d)
          nchecks <- nchecks + 1
          ncompact <- ncompact + (sum(stored) < nrow(g))
          if (!(d < 1E-14)) {
            fail(paste0("compact rhs differs: relative difference ", d, " for pars1 = ", paste(pars1, collapse = ","),
                        ", lxm = ", lxm, ", lxe = ", lxe, ", s = ", s, ", k = ", k))
          }
        }
      }
    }
  }
  expect_gt(nchecks, 400)
  expect_gt(ncompact, 100)                       # in many of them states are really not stored
  expect_lt(maxdiff, 1E-14)
})

test_that("compact: the right-hand side does not depend on the number of threads", {
  set.seed(20260202)
  pars1 <- c(0.4, 0.2, 12, 0.05, 0.3, 20)
  lxm <- 9
  lxe <- 11
  s <- 4
  k <- 5
  tri <- DAISIE_IW_tri(pars1[3], k, TRUE)
  x <- iw_test_random_state(DAISIE_IW_layout(lxm, lxe, tri)$P * 2 ^ s)
  p <- DAISIE_IW_direct_pars(pars1, pars1[3], k, lxm, lxe, s, tri)
  one <- iw_test_with_threads(1, .Call("daisie_iw_rhs_eval_direct", x, p))
  two <- iw_test_with_threads(2, .Call("daisie_iw_rhs_eval_direct", x, p))
  expect_identical(one, two)
})

test_that("compact: the event step agrees with the event step of the tables version, for every kind of event", {
  set.seed(20260203)
  parsets <- list(
    c(0.4, 0.2, 12, 0.05, 0.3, 20),
    c(0.5, 0.1, 5, 0.1, 0.2, 9),
    c(1.3, 0.7, 30.5, 0.9, 0.05, 8),
    c(0.4, 0.2, Inf, 0.05, 0.3, 25),
    c(0, 0.2, 9, 0.05, 0.3, 15),
    c(0.4, 0.2, 9, 0, 0.3, 15),
    c(0.4, 0.2, 9.5, 0.05, 0.3, 15)
  )
  settings <- list(c(lxm = 7, lxe = 9), c(lxm = 2, lxe = 2), c(lxm = 6, lxe = 3), c(lxm = 3, lxe = 6))
  maxdiff <- 0
  nchecks <- 0
  ncompact <- 0
  # embeds the compact vector of the layout tri in the box
  for (pars1 in parsets) {
    M <- pars1[6]
    for (setting in settings) {
      lxm <- setting[["lxm"]]
      lxe <- setting[["lxe"]]
      for (s in 0:4) {
        for (k in unique(c(s, s + 1, s + 3))) {
          if (k >= M - 1) next
          tri <- DAISIE_IW_tri(pars1[3], k, TRUE)
          trid <- if (is.finite(tri)) tri - 1 else Inf        # the bound after the event
          if (is.finite(tri) && trid < 0) next
          g <- iw_test_box(lxm, lxe, 2 ^ s)
          stored <- g$i + g$j <= tri
          y <- iw_test_random_state(nrow(g))
          y[!stored] <- 0
          idx <- .Call("daisie_iw_state_index_tri", g$i[stored], g$j[stored], g$a[stored], lxm, lxe, tri)
          yc <- numeric(DAISIE_IW_layout(lxm, lxe, tri)$P * 2 ^ s)
          yc[idx] <- y[stored]
          p <- DAISIE_IW_direct_pars(pars1, pars1[3], k, lxm, lxe, s, tri)
          Pd <- DAISIE_IW_layout(lxm, lxe, trid)$P
          check <- function(type, tocollapse = 1, sd) {
            ref <- iw_test_event_reference(y, type, pars1, k, lxm, lxe, s, tocollapse)
            gd <- iw_test_box(lxm, lxe, sd)
            storedd <- gd$i + gd$j <= trid
            expect_true(all(ref[!storedd] == 0))            # the result is 0 outside the triangle of the next segment
            cmp <- .Call("daisie_iw_event", yc, type, p, as.integer(tocollapse - 1))
            expect_length(cmp, sd * Pd)
            idd <- .Call("daisie_iw_state_index_tri", gd$i[storedd], gd$j[storedd], gd$a[storedd], lxm, lxe, trid)
            d <- iw_test_reldiff(ref[storedd], cmp[idd])
            maxdiff <<- max(maxdiff, d)
            nchecks <<- nchecks + 1
            ncompact <<- ncompact + (sum(storedd) < nrow(gd))
          }
          check(1L, 1, 2 * 2 ^ s)
          check(2L, 1, 2 ^ s)
          if (s >= 1) {
            for (tocollapse in 1:s) {
              check(3L, tocollapse, 2 ^ s / 2)
            }
          }
        }
      }
    }
  }
  expect_gt(nchecks, 500)
  expect_gt(ncompact, 200)
  expect_lt(maxdiff, 1E-13)
})

test_that("compact: the event step with a wrong length of the vector is an error and not a crash", {
  pars1 <- c(0.4, 0.2, 9, 0.05, 0.3, 12)
  p <- DAISIE_IW_direct_pars(pars1, pars1[3], 2, 4, 5, 2, 7)
  expect_error(.Call("daisie_iw_event", rep(0.1, 7), 1L, p, 0L))
  expect_error(.Call("daisie_iw_rhs_eval_direct", rep(0.1, 7), p))
})

test_that("compact: loglik compact = loglik direct = loglik tables on hand-made data sets, K' finite and infinite", {
  # K' = 9 is the number of events of D3, so that the triangle shrinks to the state (0, 0) in the last segment
  for (K in c(9, 12, 12.5, 20, Inf)) {
    pars1 <- c(0.3, 0.2, K, 0.02, 0.3)
    for (cond in 0:1) {
      pars2 <- c(12, 11, cond, 0)
      for (datalist in list(iw_test_D1(), iw_test_D2(), iw_test_D3())) {
        tables <- DAISIE_loglik_IW(pars1, pars2, datalist, methode = iw_test_fehlberg)
        direct <- DAISIE_loglik_IW(pars1, pars2, datalist, methode = iw_test_fehlberg, rhs = "direct")
        compact <- DAISIE_loglik_IW(pars1, pars2, datalist, methode = iw_test_fehlberg, rhs = "compact")
        expect_true(is.finite(tables))
        expect_equal(compact, direct, tolerance = 1E-12)
        expect_equal(compact, tables, tolerance = 1E-8)
      }
    }
  }
})

test_that("compact: loglik compact = loglik tables for other M, other parameters, and a K' that is too small (-Inf)", {
  datalist <- iw_test_D3()
  for (pars1 in list(c(0.8, 0.5, 20, 0.1, 0.05), c(0.2, 0.1, 15, 0.005, 1), c(0.3, 0.2, 10.5, 0.02, 0.3))) {
    for (lx in c(6, 10, 25)) {
      tables <- DAISIE_loglik_IW(pars1, c(lx, 11, 0, 0), datalist, methode = iw_test_fehlberg)
      compact <- DAISIE_loglik_IW(pars1, c(lx, 11, 0, 0), datalist, methode = iw_test_fehlberg, rhs = "compact")
      expect_true(is.finite(tables))
      expect_equal(compact, tables, tolerance = 1E-8)
    }
  }
  # K' smaller than the number of events: -Inf with the same message as with the tables
  pars1 <- c(0.3, 0.2, 8, 0.02, 0.3)
  expect_message(tables <- DAISIE_loglik_IW(pars1, c(12, 11, 0, 0), datalist, methode = iw_test_fehlberg))
  expect_message(compact <- DAISIE_loglik_IW(pars1, c(12, 11, 0, 0), datalist, methode = iw_test_fehlberg, rhs = "compact"))
  expect_identical(compact, tables)
  expect_identical(compact, -Inf)
})

test_that("compact: a larger lx than K' + 1 makes no difference (DAISIE_loglik_IW() cuts the box at K' + 1)", {
  datalist <- iw_test_D2()
  pars1 <- c(0.3, 0.2, 5, 0.02, 0.3)
  small <- DAISIE_loglik_IW(pars1, c(12, 11, 0, 0), datalist, methode = iw_test_fehlberg, rhs = "compact")
  for (lx in c(30, 200, 50000)) {
    expect_equal(DAISIE_loglik_IW(pars1, c(lx, 11, 0, 0), datalist, methode = iw_test_fehlberg, rhs = "compact"), small,
                 tolerance = 1E-12)
  }
  expect_equal(small, DAISIE_loglik_IW(pars1, c(12, 11, 0, 0), datalist, methode = iw_test_fehlberg), tolerance = 1E-8)
  # a box that is smaller than K' - k does not store less (all states are in the triangle): the same as "direct"
  pars1 <- c(0.3, 0.2, 40, 0.02, 0.3)
  expect_equal(DAISIE_loglik_IW(pars1, c(8, 11, 0, 0), datalist, methode = iw_test_fehlberg, rhs = "compact"),
               DAISIE_loglik_IW(pars1, c(8, 11, 0, 0), datalist, methode = iw_test_fehlberg, rhs = "direct"), tolerance = 1E-12)
})

test_that("compact: loglik compact = loglik tables with the other solvers of odeint", {
  skip_on_cran()
  pars1 <- c(0.3, 0.2, 12, 0.02, 0.3)
  datalist <- iw_test_D1()
  for (methode in c("odeint::runge_kutta_cash_karp54", "odeint::runge_kutta_dopri5", "odeint::bulirsch_stoer")) {
    tables <- DAISIE_loglik_IW(pars1, c(8, 11, 0, 0), datalist, methode = methode, abstolint = 1E-12, reltolint = 1E-10)
    compact <- DAISIE_loglik_IW(pars1, c(8, 11, 0, 0), datalist, methode = methode, abstolint = 1E-12, reltolint = 1E-10,
                                rhs = "compact")
    expect_equal(compact, tables, tolerance = 1E-8, label = methode)
  }
})

test_that("compact: loglik compact is the same with one and with two threads", {
  skip_on_cran()
  pars1 <- c(0.3, 0.2, 12, 0.02, 0.3)
  one <- iw_test_with_threads(1, DAISIE_loglik_IW(pars1, c(12, 11, 0, 0), iw_test_D2(), methode = iw_test_fehlberg, rhs = "compact"))
  two <- iw_test_with_threads(2, DAISIE_loglik_IW(pars1, c(12, 11, 0, 0), iw_test_D2(), methode = iw_test_fehlberg, rhs = "compact"))
  expect_identical(one, two)
})

test_that("compact: loglik compact = loglik tables for the frogs data (a large K' and a smaller one)", {
  skip_on_cran()
  utils::data(frogs_datalist, package = "DAISIE")
  for (K in c(1000.1, 100)) {
    pars1 <- c(0.2, 0.1, K, 0.001, 0.3)
    pars2 <- c(15, 11, 0, 0)
    tables <- DAISIE_loglik_IW(pars1, pars2, frogs_datalist, methode = iw_test_fehlberg)
    compact <- DAISIE_loglik_IW(pars1, pars2, frogs_datalist, methode = iw_test_fehlberg, rhs = "compact")
    expect_true(is.finite(tables))
    expect_equal(compact, tables, tolerance = 1E-8)
  }
})

test_that("compact: loglik compact = loglik tables for a data set with a recolonisation", {
  skip_on_cran()
  utils::data(frogs_sim_datalist, package = "DAISIE")
  pars1 <- c(0.4012, 0.17, 131.9, 0.00035, 0)
  pars2 <- c(20, 11, 1, 0)
  datalist2 <- frogs_sim_datalist
  datalist2[[2]]$all_colonisations[[2]]$event_times <- c(datalist2[[2]]$all_colonisations[[2]]$event_times, 5)
  for (datalist in list(frogs_sim_datalist, datalist2)) {
    tables <- DAISIE_loglik_IW(pars1, pars2, datalist, methode = iw_test_fehlberg, abstolint = 1E-12, reltolint = 1E-10)
    compact <- DAISIE_loglik_IW(pars1, pars2, datalist, methode = iw_test_fehlberg, abstolint = 1E-12, reltolint = 1E-10,
                                rhs = "compact")
    expect_true(is.finite(tables))
    expect_equal(compact, tables, tolerance = 1E-8)
  }
})

test_that("compact: loglik compact = loglik tables for simulated islands with a small K'", {
  skip_on_cran()
  set.seed(42)
  islands <- DAISIE_sim(time = 3, M = 100, pars = c(0.5, 0.3, 8, 0.01, 0.4), replicates = 6, divdepmodel = "IW", cond = 1,
                        verbose = FALSE, plot_sims = FALSE, use_rcpp = TRUE)
  pars1 <- c(0.4, 0.25, 7, 0.012, 0.3)
  pars2 <- c(12, 11, 1, 0)
  for (i in seq_along(islands)) {
    isl <- islands[[i]]
    tables <- DAISIE_loglik_IW(pars1, pars2, isl, methode = iw_test_fehlberg)
    compact <- DAISIE_loglik_IW(pars1, pars2, isl, methode = iw_test_fehlberg, rhs = "compact")
    expect_equal(compact, tables, tolerance = 1E-8, label = paste("simulated island", i))
  }
})

test_that("compact: non-finite parameters give the same result (NaN) as the tables version", {
  datalist <- iw_test_D1()
  for (pars1 in list(c(0.3, 0.2, 12, Inf, 0.3), c(Inf, 0.2, 12, 0.02, 0.3), c(0.3, 0.2, 12, 0.02, Inf))) {
    tables <- DAISIE_loglik_IW(pars1, c(10, 11, 0, 0), datalist, methode = iw_test_fehlberg)
    compact <- DAISIE_loglik_IW(pars1, c(10, 11, 0, 0), datalist, methode = iw_test_fehlberg, rhs = "compact")
    expect_identical(compact, tables)
  }
  expect_warning(compact <- DAISIE_loglik_IW(c(0.3, -0.2, 12, 0.02, 0.3), c(12, 11, 0, 0), iw_test_D2(),
                                             methode = iw_test_fehlberg, rhs = "compact"))
  expect_identical(compact, -Inf)
})

test_that("compact: impossible combinations are refused with a clear message", {
  pars1 <- c(0.3, 0.2, 12, 0.02, 0.3)
  for (methode in c("lsodes", "ode45", "lsoda")) {
    expect_error(DAISIE_loglik_IW(pars1, c(8, 11, 0, 0), iw_test_D1(), methode = methode, rhs = "compact"), "odeint")
  }
  # the sparse Jacobian of rosenbrock4 is in the box layout: rhs = "direct" works with it, "compact" does not
  expect_error(DAISIE_loglik_IW(pars1, c(8, 11, 0, 0), iw_test_D1(), methode = "odeint::rosenbrock4", rhs = "compact"),
               "rosenbrock4")
  expect_true(is.finite(DAISIE_loglik_IW(pars1, c(8, 11, 0, 0), iw_test_D1(), methode = "odeint::rosenbrock4", rhs = "direct")))
  expect_error(.Call("daisie_odeint_iw_direct", rep(0.1, 10), c(1, 0), DAISIE_IW_direct_pars(c(0.3, 0.2, 12, 0.02, 0.3, 20), 12, 2, 4, 5, 0, 3),
                     "odeint::rosenbrock4", 1E-12, 1E-10), "rosenbrock4")
})

test_that("compact: the limit of 2^31 - 1 states is for the stored states", {
  # P states per configuration and 2^s configurations
  expect_silent(DAISIE_IW_check_dim(50000, 50000, 2 ^ 5, P = 100))
  expect_silent(DAISIE_IW_check_dim(50000, 50000, 2, P = 2 ^ 30 - 1))
  expect_error(DAISIE_IW_check_dim(50000, 50000, 2, P = 2 ^ 30), "more than")
  expect_error(DAISIE_IW_check_dim(50000, 50000, 2 ^ 5, P = 2 ^ 27), "states per configuration")
  # through the likelihood: a large K' (the box is not cut), the check is before anything is allocated, as for "direct"
  datalist <- iw_test_D2()
  datalist[[1]]$not_present <- 1E5
  expect_error(DAISIE_loglik_IW(c(0.3, 0.2, 1E6, 0.02, 0.3), c(50000, 11, 0, 0), datalist, methode = iw_test_fehlberg,
                                rhs = "compact"), "2\\^31 - 1")
  # (a data set that is too large for "direct" and small enough for "compact" would need about 13 GB to run, so that case is
  # only tested through DAISIE_IW_check_dim() above)
})

test_that("compact: the ML estimates are the same with rhs = tables and rhs = compact", {
  skip_on_cran()
  datalist <- iw_test_data(5, 20, list(list(stac = 4, bt = 4.0), list(stac = 2, bt = c(3.0, 1.5)), list(stac = 4, bt = 0.5)))
  args <- list(datalist = datalist, initparsopt = c(0.5, 0.3), idparsopt = c(1, 2), parsfix = c(12, 0.05, 0.3),
               idparsfix = c(3, 4, 5), res = 8, ddmodel = 11, cond = 0, methode = iw_test_fehlberg, tolint = c(1E-12, 1E-10),
               tol = c(1E-3, 1E-4, 1E-6), maxiter = 200)
  tables <- suppressMessages(do.call(DAISIE_ML_IW, args))
  compact <- suppressMessages(do.call(DAISIE_ML_IW, c(args, list(rhs = "compact"))))
  expect_true(is.finite(tables$loglik))
  expect_equal(compact$loglik, tables$loglik, tolerance = 1E-8)
  expect_equal(unlist(compact[1:5]), unlist(tables[1:5]), tolerance = 1E-6)
  expect_equal(compact$conv, tables$conv)
})
