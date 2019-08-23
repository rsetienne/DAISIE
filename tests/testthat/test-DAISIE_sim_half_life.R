context("DAISIE_sim_half_life")

test_that("DAISIE_sim_half_life gives output for nonoceanic IW", {
  replicates <- 1
  #change the mu as mu > lambda_c is invalid
  set.seed(1)
  half_life <- DAISIE_sim_half_life(time = 1, replicates = replicates, lac = 1, 
                                    mu = 3, K = 40, gam = 0.1, laa = 1, ssr = 2,
                                    divdepmodel = "IW", ddmodel_sim = 11,
                                    island_type = "nonoceanic", x_s = 0.1,
                                    x_nonend = 0.9, verbose = FALSE)
  expect_true(is.list(half_life))
  expect_length(half_life[[1]], replicates)
  expect_true(is.numeric(half_life[[1]][[1]]))
})


test_that("DAISIE_sim_half_life gives output for nonoceanic CS", {
  replicates <- 1
  #change the mu as mu > lambda_c is invalid
  set.seed(1)
  half_life <- DAISIE_sim_half_life(time = 1, replicates = replicates, lac = 1,
                                    mu = 5, K = 2, gam = 0.1, laa = 1, ssr = 2,
                                    divdepmodel = "CS", ddmodel_sim = 11,
                                    island_type = "nonoceanic", x_s = 0.9,
                                    x_nonend = 0.9, verbose = FALSE)
  expect_true(is.list(half_life))
  expect_length(half_life[[1]], replicates)
  expect_true(is.numeric(half_life[[1]][[1]]))
})

test_that("DAISIE_sim_half_life gives output for oceanic IW", {
  replicates <- 1
  set.seed(1)
  half_life <- DAISIE_sim_half_life(time = 1, replicates = replicates, lac = 1,
                                    mu = 0.9, K = 40, gam = 0.1, laa = 1, 
                                    ssr = 2, divdepmodel = "IW", 
                                    ddmodel_sim = 11, island_type = "oceanic",
                                    x_s = NULL, x_nonend = NULL,
                                    verbose = FALSE)
  expect_true(is.list(half_life))
  expect_length(half_life, replicates)
  expect_true(is.numeric(half_life[[1]][[1]]))
})

test_that("DAISIE_sim_half_life gives output for oceanic CS", {
skip("test is too long")  
  set.seed(1)
  half_life <- DAISIE_sim_half_life(time = 1, replicates = 1, lac = 1, mu = 0.9,
                                    K = 40, gam = 0.1, laa = 1, ssr = 2,
                                    divdepmodel = "CS", ddmodel_sim = 11,
                                    island_type = "oceanic", x_s = NULL,
                                    x_nonend = NULL, verbose = FALSE)
  expect_true(is.list(half_life))
  expect_length(half_life, replicates)
  expect_true(is.numeric(half_life[[1]][[1]]))
})
