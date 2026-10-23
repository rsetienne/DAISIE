test_that("DAISIE_A_ExpEIN equals DAISIE_ExpEIN for one island", {
  p <- c(0.7, 0.3, Inf, 0.05, 0.4)
  pars <- DAISIE_A_create_pars(n_islands = 1, lac = p[1], mu = p[2],
                               gam = p[4], laa = p[5])
  a <- DAISIE_ExpEIN(tvec = c(1, 2, 5), pars = p, M = 100)
  b <- DAISIE_A_ExpEIN(tvec = c(1, 2, 5), pars = pars, M = 100)
  expect_equal(b$ExpE, a$ExpE, tolerance = 1e-7)
  expect_equal(b$ExpI, a$ExpI, tolerance = 1e-7)
  expect_equal(b$ExpN, a$ExpN, tolerance = 1e-7)
})

test_that("lumped system equals the full system for exchangeable islands", {
  for (vm in c("all", "single")) {
    for (cs in c("split", "bud")) {
      pars <- DAISIE_A_create_pars(
        n_islands = 4, lac = 0.4, mu = 0.3, gam = 0.03, laa = 0.2,
        lav = 0.25, disp = 0.15, laj = 0.05, vic_mode = vm,
        clado_single = cs, max_range_size = 3
      )
      full <- DAISIE_A_ExpEIN(c(1, 3), pars, M = 50)
      lumped <- DAISIE_A_ExpEIN(c(1, 3), pars, M = 50, lumped = TRUE)
      sizes <- vapply(DAISIE_A_ranges(pars)$sets, length, numeric(1))
      for (k in 1:3) {
        expect_equal(rowSums(full$ExpE_range[, sizes == k, drop = FALSE]),
                     unname(lumped$ExpE_range[, k]), tolerance = 1e-6)
        expect_equal(rowSums(full$ExpI_range[, sizes == k, drop = FALSE]),
                     unname(lumped$ExpI_range[, k]), tolerance = 1e-6)
      }
      expect_equal(full$ExpN_island, lumped$ExpN_island, tolerance = 1e-6)
      expect_equal(full$ExpSIE_island, lumped$ExpSIE_island,
                   tolerance = 1e-6)
    }
  }
  pars <- DAISIE_A_create_pars(n_islands = 3, lac = c(1, 2, 3), mu = 1,
                               gam = 1, laa = 1)
  expect_error(DAISIE_A_ExpEIN(1, pars, M = 1, lumped = TRUE),
               "exchangeable")
})

test_that("simulated means agree with expected numbers per range", {
  skip_on_cran()
  pars <- DAISIE_A_create_pars(
    n_islands = 3, lac = c(0.3, 0.5, 0.7), mu = c(0.2, 0.3, 0.4),
    gam = c(0.02, 0.03, 0.01), laa = 0.3, lav = 0.4, disp = 0.4, laj = 0.1,
    connectivity = "line"
  )
  reps <- 100
  M <- 50
  set.seed(42)
  sims <- DAISIE_A_sim(time = 2, M = M, pars = pars, replicates = reps,
                       verbose = FALSE)
  rs <- DAISIE_A_ranges(pars)$strings
  counts <- t(sapply(sims, function(rep) {
    sp <- do.call(rbind, c(list(data.frame(range = character(0),
                                           type = character(0))),
                           lapply(rep[-1], function(x) {
                             x$species[, c("range", "type")]
                           })))
    c(table(factor(sp$range[sp$type == "I"], levels = rs)),
      table(factor(sp$range[sp$type != "I"], levels = rs)))
  }))
  expected <- DAISIE_A_ExpEIN(tvec = 2, pars = pars, M = M)
  expected <- c(expected$ExpI_range[1, ], expected$ExpE_range[1, ])
  se <- apply(counts, 2, stats::sd) / sqrt(reps)
  z <- (colMeans(counts) - expected) / se
  z <- z[se > 0]
  # chi-square test of all ranges and types jointly
  expect_gt(stats::pchisq(sum(z^2), df = length(z), lower.tail = FALSE),
            0.001)
})
