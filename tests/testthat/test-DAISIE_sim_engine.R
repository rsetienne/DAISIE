# The C++ simulation engine draws random numbers in the same order as the R
# code. Hence, with the same seed, both implementations must give identical
# simulations. These tests check that for each simulation core.

expect_same_sim <- function(f, seeds = 1:5) {
  for (s in seeds) {
    set.seed(s)
    r <- suppressMessages(f(FALSE))
    set.seed(s)
    cpp <- suppressMessages(f(TRUE))
    expect_equal(cpp, r)
  }
}

test_that("C++ and R simulation cores agree: constant rates, non-oceanic", {
  skip_on_cran()
  expect_same_sim(function(use_rcpp) {
    DAISIE:::DAISIE_sim_core_cr(
      time = 5,
      mainland_n = 20,
      pars = c(0.5, 0.3, 15, 0.05, 0.4),
      nonoceanic_pars = c(0.3, 0.5),
      hyper_pars = create_hyper_pars(d = 0.1, x = 0.2),
      area_pars = create_area_pars(1, 1, 0, 0, 0, 0, 0),
      use_rcpp = use_rcpp
    )
  }, seeds = 1:15)
})

test_that("C++ and R simulation cores agree: rate shift", {
  skip_on_cran()
  pars <- c(0.5, 0.3, 15, 0.05, 0.4, 0.7, 0.2, 25, 0.1, 0.2)
  for (shift_times in list(c(4, 2), c(6, 2))) { # the second shifts at time 0
    expect_same_sim(function(use_rcpp) {
      DAISIE:::DAISIE_sim_core_cr_shift(
        time = 6,
        mainland_n = 3,
        pars = pars,
        hyper_pars = create_hyper_pars(d = 0, x = 0),
        area_pars = create_area_pars(1, 1, 0, 0, 0, 0, 0),
        shift_times = shift_times,
        use_rcpp = use_rcpp
      )
    }, seeds = 1:10)
  }
})

test_that("C++ and R simulation cores agree: split-rates (SR) model", {
  skip_on_cran()
  pars <- c(0.5, 0.3, 15, 0.05, 0.4, 0.7, 0.2, 25, 0.1, 0.2, 3)
  expect_same_sim(function(use_rcpp) {
    DAISIE:::DAISIE_SR_sim_core(
      time = 6, mainland_n = 1, pars = pars, use_rcpp = use_rcpp
    )
  }, seeds = 1:10)
})

test_that("C++ and R simulation cores agree: time-dependent", {
  skip_on_cran()
  base <- c(max_area = 5000, current_area = 1000, proportional_peak_t = 0.5,
            total_island_age = 5)
  total_time <- 4
  for (cfg in list(c(1, 0), c(0, 1), c(1, 1))) {
    area_pars <- create_area_pars(
      max_area = base[["max_area"]],
      current_area = base[["current_area"]],
      proportional_peak_t = base[["proportional_peak_t"]],
      total_island_age = base[["total_island_age"]],
      sea_level_amplitude = if (cfg[2] == 1) 20 else 0,
      sea_level_frequency = if (cfg[2] == 1) 10 else 0,
      island_gradient_angle = if (cfg[2] == 1) 20 else 0
    )
    peak <- if (cfg[1] == 1) DAISIE:::calc_peak(total_time, area_pars) else NULL
    Amax <- DAISIE:::get_global_max_area(total_time, area_pars, peak, cfg[1], cfg[2])
    Amin <- DAISIE:::get_global_min_area(total_time, area_pars, peak, cfg[1], cfg[2])
    expect_same_sim(function(use_rcpp) {
      DAISIE:::DAISIE_sim_core_time_dep(
        time = total_time,
        mainland_n = 1,
        pars = c(0.5, 0.3, 10, 0.05, 0.4),
        nonoceanic_pars = c(0, 0),
        island_ontogeny = cfg[1],
        sea_level = cfg[2],
        hyper_pars = create_hyper_pars(d = 0.1, x = 0.2),
        area_pars = area_pars,
        peak = peak,
        Amax = Amax,
        Amin = Amin,
        use_rcpp = use_rcpp
      )
    }, seeds = 1:5)
  }
})

test_that("C++ and R simulation cores agree: trait-dependent", {
  skip_on_cran()
  trait_pars <- create_trait_pars(
    trans_rate = 0.3, immig_rate2 = 0.02, ext_rate2 = 0.5, ana_rate2 = 0.2,
    clado_rate2 = 0.4, trans_rate2 = 0.2, M2 = 3
  )
  trait_pars_2K <- DAISIE:::create_trait_pars_2K(
    trans_rate = 0.3, immig_rate2 = 0.02, ext_rate2 = 0.5, ana_rate2 = 0.2,
    clado_rate2 = 0.4, trans_rate2 = 0.2, M2 = 3, K2 = 8
  )
  for (tp in list(trait_pars, trait_pars_2K)) {
    expect_same_sim(function(use_rcpp) {
      DAISIE:::DAISIE_sim_core_trait_dep(
        time = 5,
        mainland_n = 2,
        pars = c(0.5, 0.3, 15, 0.05, 0.4),
        hyper_pars = create_hyper_pars(d = 0, x = 0),
        area_pars = create_area_pars(1, 1, 0, 0, 0, 0, 0),
        trait_pars = tp,
        use_rcpp = use_rcpp
      )
    }, seeds = 1:10)
  }
})

test_that("use_rcpp is passed on by the top level simulation functions", {
  skip_on_cran()
  area_pars <- create_area_pars(
    max_area = 5000, current_area = 1000, proportional_peak_t = 0.5,
    total_island_age = 5, sea_level_amplitude = 0, sea_level_frequency = 0,
    island_gradient_angle = 0
  )
  for (divdepmodel in c("CS", "IW")) {
    expect_same_sim(function(use_rcpp) {
      DAISIE_sim_time_dep(
        time = 4, M = 4, pars = c(0.5, 0.3, 10, 0.05, 0.4), replicates = 2,
        area_pars = area_pars, hyper_pars = create_hyper_pars(0.1, 0.2),
        island_ontogeny = "beta", divdepmodel = divdepmodel,
        plot_sims = FALSE, verbose = FALSE, use_rcpp = use_rcpp
      )
    }, seeds = 1:2)
    expect_same_sim(function(use_rcpp) {
      DAISIE_sim_cr_shift(
        time = 6, M = 4,
        pars = c(0.5, 0.3, 15, 0.05, 0.4, 0.7, 0.2, 25, 0.1, 0.2),
        replicates = 2, shift_times = c(4, 2), divdepmodel = divdepmodel,
        plot_sims = FALSE, verbose = FALSE, use_rcpp = use_rcpp
      )
    }, seeds = 1:2)
  }
  trait_pars <- create_trait_pars(
    trans_rate = 0.3, immig_rate2 = 0.02, ext_rate2 = 0.5, ana_rate2 = 0.2,
    clado_rate2 = 0.4, trans_rate2 = 0.2, M2 = 3
  )
  expect_same_sim(function(use_rcpp) {
    DAISIE_sim_trait_dep(
      time = 5, M = 3, pars = c(0.5, 0.3, 15, 0.05, 0.4), replicates = 2,
      trait_pars = trait_pars, plot_sims = FALSE, verbose = FALSE,
      use_rcpp = use_rcpp
    )
  }, seeds = 1:2)
})

test_that("C++ time-dependent engine stops if a rate exceeds its bound", {
  skip_on_cran()
  # Amax deliberately too small: the clado/immigration bound is violated
  area_pars <- create_area_pars(
    max_area = 5000, current_area = 1000, proportional_peak_t = 0.5,
    total_island_age = 5, sea_level_amplitude = 0, sea_level_frequency = 0,
    island_gradient_angle = 0
  )
  peak <- DAISIE:::calc_peak(4, area_pars)
  expect_error(
    {
      set.seed(1)
      DAISIE:::DAISIE_sim_core_time_dep(
        time = 4, mainland_n = 5, pars = c(0.5, 0.3, 10, 0.5, 0.4),
        nonoceanic_pars = c(0, 0), island_ontogeny = 1, sea_level = 0,
        hyper_pars = create_hyper_pars(d = 0.1, x = 0.2),
        area_pars = area_pars, peak = peak, Amax = 1, Amin = 0.5,
        area_windows = list(end = 4, Amin = 0.5, Amax = 1),
        use_rcpp = TRUE
      )
    },
    "exceeds its upper bound"
  )
})
