test_that("DAISIE_A_create_pars expands scalars and builds matrices", {
  pars <- DAISIE_A_create_pars(
    n_islands = 3, lac = 0.5, mu = c(0.1, 0.2, 0.3), K = 10, gam = 0.02,
    laa = 0.4, laa_isl = c(0.05, 0, 0.1), lav = 0.1, disp = 0.3,
    connectivity = "line"
  )
  expect_true(DAISIE_A_are_pars(pars))
  expect_equal(pars$lac, rep(0.5, 3))
  expect_equal(pars$mu, c(0.1, 0.2, 0.3))
  expect_equal(pars$K, rep(10, 3))
  expect_equal(pars$D, 0.3 * DAISIE_A_connectivity(3, "line"))
  expect_equal(pars$laa_isl, c(0.05, 0, 0.1))
  expect_equal(DAISIE_A_create_pars(n_islands = 2, lac = 1, mu = 1, gam = 1,
                                    laa = 1, laa_isl = 0.3)$laa_isl,
               c(0.3, 0.3))
  expect_equal(pars$adjacency[1, ], c(FALSE, TRUE, FALSE))
  expect_equal(pars$max_range_size, 3L)
  expect_equal(pars$vic_per, "split")
  expect_equal(pars$vic_kappa, 0)
  expect_null(pars$submergence_ages)
  expect_equal(DAISIE_A_create_pars(n_islands = 2, lac = 1, mu = 1, gam = 1,
                                    laa = 1, vic_per = "lineage")$vic_per,
               "lineage")
  # all-NA submergence ages mean no submergence
  expect_null(DAISIE_A_create_pars(n_islands = 2, lac = 1, mu = 1, gam = 1,
                                   laa = 1,
                                   submergence_ages = c(NA, NA))$submergence_ages)
  expect_output(print(pars), "archipelago of 3 islands")
})

test_that("DAISIE_A_create_pars accepts a general connectivity matrix", {
  C <- matrix(c(0, 1, 0.2,
                0.5, 0, 0,
                0, 0, 0), 3, 3, byrow = TRUE)
  pars <- DAISIE_A_create_pars(
    n_islands = 3, lac = 1, mu = 1, gam = 1, laa = 1, disp = 2,
    connectivity = C
  )
  expect_equal(pars$D, 2 * C)
  # adjacency is symmetric: 1-2 and 1-3 adjacent, 2-3 not
  expect_equal(pars$adjacency,
               matrix(c(FALSE, TRUE, TRUE,
                        TRUE, FALSE, FALSE,
                        TRUE, FALSE, FALSE), 3, 3, byrow = TRUE))
})

test_that("DAISIE_A_create_pars catches wrong input", {
  expect_error(DAISIE_A_create_pars(n_islands = 3, lac = c(1, 2), mu = 1,
                                    gam = 1, laa = 1), "length")
  expect_error(DAISIE_A_create_pars(n_islands = 3, lac = 1, mu = -1,
                                    gam = 1, laa = 1), "non-negative")
  expect_error(DAISIE_A_create_pars(n_islands = 2, lac = 1, mu = 1,
                                    gam = 1, laa = 1,
                                    connectivity = matrix(1, 3, 3)),
               "connectivity")
  expect_error(DAISIE_A_create_pars(n_islands = 2, lac = 1, mu = 1,
                                    gam = 1, laa = 1, vic_mode = "foo"))
  expect_error(DAISIE_A_create_pars(n_islands = 2, lac = 1, mu = 1,
                                    gam = 1, laa = 1, vic_per = "foo"))
  expect_error(DAISIE_A_create_pars(n_islands = 2, lac = 1, mu = 1,
                                    gam = 1, laa = 1, vic_kappa = -1),
               "vic_kappa")
  expect_error(DAISIE_A_create_pars(n_islands = 2, lac = 1, mu = 1,
                                    gam = 1, laa = 1,
                                    submergence_ages = 1), "submergence")
  expect_error(DAISIE_A_create_pars(n_islands = 2, lac = 1, mu = 1,
                                    gam = 1, laa = 1, island_ages = c(3, 2),
                                    submergence_ages = c(NA, 2.5)),
               "emerge before")
  # with connected ranges, an island in the middle of a line cannot submerge
  expect_error(DAISIE_A_create_pars(n_islands = 3, lac = 1, mu = 1, gam = 1,
                                    laa = 1, connectivity = "line",
                                    connected_ranges = TRUE,
                                    submergence_ages = c(NA, 1, NA)),
               "neighbours")
  expect_silent(DAISIE_A_create_pars(n_islands = 3, lac = 1, mu = 1, gam = 1,
                                     laa = 1, connectivity = "line",
                                     connected_ranges = TRUE,
                                     submergence_ages = c(1, NA, NA)))
  expect_error(DAISIE_A_create_pars(n_islands = 2, lac = 1, mu = 1,
                                    gam = 1, laa = 1,
                                    island_ages = c(1, 2, 3)), "island_ages")
  expect_error(DAISIE_A_sim(time = 1, M = 1, pars = c(1, 1, 1, 1, 1),
                            replicates = 1), "DAISIE_A_create_pars")
})

test_that("connectivity helpers", {
  expect_equal(DAISIE_A_connectivity(3, "full"),
               matrix(1, 3, 3) - diag(3))
  ring <- DAISIE_A_connectivity(4, "ring")
  expect_equal(rowSums(ring), rep(2, 4))
  expect_equal(sum(DAISIE_A_connectivity(4, "none")), 0)
  d <- matrix(c(0, 1, 2,
                1, 0, 1,
                2, 1, 0), 3, 3)
  C <- DAISIE_A_connectivity_distance(d, scale = 1)
  expect_equal(C[1, 3], exp(-2))
  expect_equal(diag(C), rep(0, 3))
  C2 <- DAISIE_A_connectivity_distance(d, scale = 1, threshold = 1.5)
  expect_equal(C2[1, 3], 0)
})
