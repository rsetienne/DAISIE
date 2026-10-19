# Loops over a group of species or clades must run zero times, not twice, when
# the group is empty (1:0 and (n + 1):n count down in R).

area_pars_const <- function() {
  create_area_pars(
    max_area = 1,
    current_area = 1,
    proportional_peak_t = 0,
    total_island_age = 0,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0
  )
}

# Every mainland species colonises (huge immigration rate, no other events),
# so the number of colonist clades on the island is the number of simulated
# mainland species.
n_clades <- function(island_replicate) {
  length(island_replicate) - 1
}

test_that("DAISIE_sim_cr_cs simulates exactly M species when type 2 pool is empty", {
  skip_on_cran()
  M <- 10
  set.seed(1)
  sim <- DAISIE_sim_cr_cs(
    total_time = 1,
    M = M,
    pars = c(1e-9, 1e-9, 20, 1e4, 1e-9, 1e-9, 1e-9, 20, 1e4, 1e-9),
    replicates = 1,
    nonoceanic_pars = c(0, 0),
    prop_type2_pool = 0.01, # M * 0.01 rounds to 0 type 2 species
    replicates_apply_type2 = FALSE,
    sample_freq = 5,
    hyper_pars = create_hyper_pars(d = 0, x = 0),
    area_pars = area_pars_const(),
    cond = 0,
    verbose = FALSE
  )
  testthat::expect_equal(n_clades(sim[[1]]), M)
})

test_that("DAISIE_sim_min_type2 runs when type 1 pool is empty", {
  skip_on_cran()
  M <- 5
  set.seed(1)
  sim <- DAISIE_sim_min_type2(
    time = 1,
    M = M,
    pars = c(1e-9, 1e-9, 20, 1e4, 1e-9, 1e-9, 1e-9, 20, 1e4, 1e-9),
    replicates = 1,
    prop_type2_pool = 1,
    area_pars = area_pars_const(),
    hyper_pars = create_hyper_pars(d = 0, x = 0),
    verbose = FALSE
  )
  testthat::expect_length(sim, 1)
  testthat::expect_length(sim[[1]], M)
})

test_that("DAISIE_sim_min_type2 refuses an empty type 2 pool instead of looping forever", {
  testthat::expect_error(
    DAISIE_sim_min_type2(
      time = 1,
      M = 10,
      pars = c(0.5, 0.1, 20, 0.01, 0.1, 0.9, 0.1, 20, 0.01, 0.1),
      replicates = 1,
      prop_type2_pool = 0.01,
      area_pars = area_pars_const(),
      hyper_pars = create_hyper_pars(d = 0, x = 0),
      verbose = FALSE
    ),
    regexp = "type 2 pool is empty"
  )
})

test_that("DAISIE_sim_trait_dep simulates M species when M2 = 0", {
  skip_on_cran()
  M <- 5
  set.seed(1)
  trait_pars <- create_trait_pars(
    trans_rate = 0,
    immig_rate2 = 1e4,
    ext_rate2 = 1e-9,
    ana_rate2 = 1e-9,
    clado_rate2 = 1e-9,
    trans_rate2 = 0,
    M2 = 0
  )
  sim <- DAISIE_sim_trait_dep(
    time = 1,
    M = M,
    pars = c(1e-9, 1e-9, 20, 1e4, 1e-9),
    replicates = 1,
    island_ontogeny = "const",
    sea_level = "const",
    trait_pars = trait_pars,
    extcutoff = 1000,
    plot_sims = FALSE,
    verbose = FALSE
  )
  testthat::expect_equal(n_clades(sim[[1]]), M)
})

test_that("DAISIE_DE_loglik_CS handles an island without colonists", {
  pars1 <- c(0.5, 0.3, 20, 0.05, 0.8)
  not_present <- 5
  datalist <- list(list(island_age = 1, not_present = not_present))
  loglik <- DAISIE_DE_loglik_CS(
    pars1 = pars1,
    pars2 = c(20, 11, 0, 0, NA, NA),
    datalist = datalist,
    methode = "ode45",
    equal_extinction = TRUE
  )
  pars1_equal_ext <- pars1
  pars1_equal_ext[3] <- pars1_equal_ext[2] # equal_extinction = TRUE
  logp0 <- DAISIE_DE_logp0(
    island_age = 1,
    pars1 = pars1_equal_ext,
    methode = "ode45",
    reltolint = 1e-10,
    abstolint = 1e-16
  )
  testthat::expect_equal(loglik, not_present * logp0, tolerance = 1e-6)
})
