test_that("A relaxed-cladogenesis should run silent with correct output", {
  set.seed(1)
  replicates <- 1
  testthat::expect_silent(
    sim <- DAISIE_sim_relaxed_rate(
      time = 5,
      M = 100,
      pars = c(5, 1, 10, 0.01, 1, 5),
      replicates = replicates,
      relaxed_par = "cladogenesis",
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
  testthat::expect_true(is.list(sim))
  testthat::expect_true(length(sim) == replicates)
  testthat::expect_equal(sim[[1]][[1]]$island_age, 5)
  testthat::expect_equal(sim[[1]][[1]]$not_present, 97)
  testthat::expect_true(is.matrix(sim[[1]][[1]]$stt_all))
  testthat::expect_equal(nrow(sim[[1]][[1]]$stt_all), 26)
  testthat::expect_equal(ncol(sim[[1]][[1]]$stt_all), 5)
  testthat::expect_length(sim[[1]], 4)
  testthat::expect_equal(sim[[1]][[2]]$branching_times,
               c(5, 2.05346943810582, 1.70906507509317, 1.66251474618355,
                 1.58428007690634, 1.51035828260695, 0.93816016264294,
                 0.88268838887265, 0.75636093416664, 0.0135436428997))
  testthat::expect_equal(sim[[1]][[2]]$stac, 2)
  testthat::expect_equal(sim[[1]][[2]]$missing_species, 0)
})

test_that("A relaxed-cladogenesis should cond run silent with correct output", {
  set.seed(Sys.time())
  replicates <- 1
  pars <- c(5, 1, 10, 0.01, 1, 5)
  cond <- 5
  testthat::expect_silent(
    out <- DAISIE_sim_relaxed_rate(
      time = 5,
      M = 100,
      pars = pars,
      replicates = replicates,
      relaxed_par = "cladogenesis",
      plot_sims = FALSE,
      verbose = FALSE,
      cond = cond
    )
  )
testthat::expect_true(out[[1]][[1]]$stt_all[nrow(out[[1]][[1]]$stt_all), 5] >= cond)
})

test_that("A relaxed-cladogenesis should cond run silent with correct output", {
  set.seed(1)
  replicates <- 1
  pars <- c(5, 1, 10, 0.01, 1, 5)
  cond <- 5
  testthat::expect_silent(
    out_cond <- DAISIE_sim_relaxed_rate(
      time = 5,
      M = 100,
      pars = pars,
      replicates = replicates,
      relaxed_par = "cladogenesis",
      plot_sims = FALSE,
      verbose = FALSE,
      cond = cond
    )
  )

  testthat::expect_true(
    out_cond[[1]][[1]]$stt_all[nrow(out_cond[[1]][[1]]$stt_all), 5] >= cond
  )

  set.seed(1)
  testthat::expect_silent(
    out_no_cond <- DAISIE_sim_relaxed_rate(
      time = 5,
      M = 100,
      pars = pars,
      replicates = replicates,
      relaxed_par = "cladogenesis",
      plot_sims = FALSE,
      verbose = FALSE,
      cond = 0
    )
  )
  testthat::expect_true(
    out_no_cond[[1]][[1]]$stt_all[nrow(out_no_cond[[1]][[1]]$stt_all), 5] < 5
  )

})

test_that("A relaxed-extinction should run silent with correct output", {
  set.seed(1)
  replicates <- 1
  testthat::expect_silent(
    sim <- DAISIE_sim_relaxed_rate(
      time = 5,
      M = 100,
      pars = c(1, 0.5, 10, 0.01, 1, 1),
      replicates = replicates,
      relaxed_par = "extinction",
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
  testthat::expect_true(is.list(sim))
  testthat::expect_true(length(sim) == replicates)
  testthat::expect_equal(sim[[1]][[1]]$island_age, 5)
  testthat::expect_equal(sim[[1]][[1]]$not_present, 98)
  testthat::expect_true(is.matrix(sim[[1]][[1]]$stt_all))
  testthat::expect_equal(nrow(sim[[1]][[1]]$stt_all), 26)
  testthat::expect_equal(ncol(sim[[1]][[1]]$stt_all), 5)
  testthat::expect_length(sim[[1]], 3)
  testthat::expect_equal(sim[[1]][[2]]$branching_times,
               c(5, 1.27314736135304, 0.8347663550417, 0.50364880811414,
                 0.00577353632018962))
  testthat::expect_equal(sim[[1]][[2]]$stac, 2)
  testthat::expect_equal(sim[[1]][[2]]$missing_species, 0)
})

test_that("A relaxed-extinction should cond run silent with correct output", {
  set.seed(Sys.time())
  replicates <- 1
  pars <- c(1, 0.5, 10, 0.01, 1, 1)
  cond <- 5
  testthat::expect_silent(
    out <- DAISIE_sim_relaxed_rate(
      time = 5,
      M = 100,
      pars = pars,
      replicates = replicates,
      relaxed_par = "extinction",
      plot_sims = FALSE,
      verbose = FALSE,
      cond = cond
    )
  )
  testthat::expect_true(out[[1]][[1]]$stt_all[nrow(out[[1]][[1]]$stt_all), 5] >= cond)
})

test_that("A relaxed-extinction should cond run silent with correct output", {
  set.seed(1)
  replicates <- 1
  pars <- c(1, 0.5, 10, 0.01, 1, 1)
  cond <- 5
  testthat::expect_silent(
    out_cond <- DAISIE_sim_relaxed_rate(
      time = 5,
      M = 100,
      pars = pars,
      replicates = replicates,
      relaxed_par = "extinction",
      plot_sims = FALSE,
      verbose = FALSE,
      cond = cond
    )
  )

  testthat::expect_true(
    out_cond[[1]][[1]]$stt_all[nrow(out_cond[[1]][[1]]$stt_all), 5] >= cond
  )

  set.seed(1)
  testthat::expect_silent(
    out_no_cond <- DAISIE_sim_relaxed_rate(
      time = 5,
      M = 100,
      pars = pars,
      replicates = replicates,
      relaxed_par = "extinction",
      plot_sims = FALSE,
      verbose = FALSE,
      cond = 0
    )
  )
  testthat::expect_true(
    out_no_cond[[1]][[1]]$stt_all[nrow(out_no_cond[[1]][[1]]$stt_all), 5] < 5
  )

})


test_that("A relaxed-K should run silent with correct output", {
  set.seed(1)
  replicates <- 1
  testthat::expect_silent(
    sim <- DAISIE_sim_relaxed_rate(
      time = 5,
      M = 100,
      pars = c(1, 1, 5, 0.01, 1, 5),
      replicates = replicates,
      relaxed_par = "carrying_capacity",
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
  testthat::expect_true(is.list(sim))
  testthat::expect_true(length(sim) == replicates)
  testthat::expect_equal(sim[[1]][[1]]$island_age, 5)
  testthat::expect_equal(sim[[1]][[1]]$not_present, 99)
  testthat::expect_true(is.matrix(sim[[1]][[1]]$stt_all))
  testthat::expect_equal(nrow(sim[[1]][[1]]$stt_all), 26)
  testthat::expect_equal(ncol(sim[[1]][[1]]$stt_all), 5)
  testthat::expect_length(sim[[1]], 2)
  testthat::expect_equal(sim[[1]][[2]]$branching_times,
               c(5.00000000000000, 2.05346943810582))
  testthat::expect_equal(sim[[1]][[2]]$stac, 2)
  testthat::expect_equal(sim[[1]][[2]]$missing_species, 0)
})

test_that("A relaxed-K should cond run silent with correct output", {
  set.seed(Sys.time())
  replicates <- 1
  pars <- c(1, 1, 5, 0.01, 1, 5)
  cond <- 5
  testthat::expect_silent(
    out <- DAISIE_sim_relaxed_rate(
      time = 5,
      M = 100,
      pars = pars,
      replicates = replicates,
      relaxed_par = "carrying_capacity",
      plot_sims = FALSE,
      verbose = FALSE,
      cond = cond
    )
  )
  testthat::expect_true(out[[1]][[1]]$stt_all[nrow(out[[1]][[1]]$stt_all), 5] >= cond)
})

test_that("A relaxed-K should cond run silent with correct output", {
  set.seed(1)
  replicates <- 1
  pars <- c(1, 1, 5, 0.01, 1, 5)
  cond <- 5
  testthat::expect_silent(
    out_cond <- DAISIE_sim_relaxed_rate(
      time = 5,
      M = 100,
      pars = pars,
      replicates = replicates,
      relaxed_par = "carrying_capacity",
      plot_sims = FALSE,
      verbose = FALSE,
      cond = cond
    )
  )

  testthat::expect_true(
    out_cond[[1]][[1]]$stt_all[nrow(out_cond[[1]][[1]]$stt_all), 5] >= cond
  )

  set.seed(1)
  testthat::expect_silent(
    out_no_cond <- DAISIE_sim_relaxed_rate(
      time = 5,
      M = 100,
      pars = pars,
      replicates = replicates,
      relaxed_par = "carrying_capacity",
      plot_sims = FALSE,
      verbose = FALSE,
      cond = 0
    )
  )
  testthat::expect_true(
    out_no_cond[[1]][[1]]$stt_all[nrow(out_no_cond[[1]][[1]]$stt_all), 5] < 5
  )

})

test_that("A relaxed-immigration should run silent with correct output", {
  set.seed(3141)
  replicates <- 1
  testthat::expect_silent(
    sim <- DAISIE_sim_relaxed_rate(
      time = 5,
      M = 100,
      pars = c(1, 1, 10, 1, 1, 5),
      replicates = replicates,
      relaxed_par = "immigration",
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
  testthat::expect_true(is.list(sim))
  testthat::expect_true(length(sim) == replicates)
  testthat::expect_equal(sim[[1]][[1]]$island_age, 5)
  testthat::expect_equal(sim[[1]][[1]]$not_present, 90)
  testthat::expect_true(is.matrix(sim[[1]][[1]]$stt_all))
  testthat::expect_equal(nrow(sim[[1]][[1]]$stt_all), 26)
  testthat::expect_equal(ncol(sim[[1]][[1]]$stt_all), 5)
  testthat::expect_length(sim[[1]], 11)
  testthat::expect_equal(sim[[1]][[2]]$branching_times,
               c(5, 2.10054118799481))
  testthat::expect_equal(sim[[1]][[2]]$stac, 2)
  testthat::expect_equal(sim[[1]][[2]]$missing_species, 0)
})

test_that("A relaxed-immigration should cond run silent with correct output", {
  set.seed(Sys.time())
  replicates <- 1
  pars <- c(1, 1, 10, 1, 1, 5)
  cond <- 5
  testthat::expect_silent(
    out <- DAISIE_sim_relaxed_rate(
      time = 5,
      M = 100,
      pars = pars,
      replicates = replicates,
      relaxed_par = "immigration",
      plot_sims = FALSE,
      verbose = FALSE,
      cond = cond
    )
  )
  testthat::expect_true(out[[1]][[1]]$stt_all[nrow(out[[1]][[1]]$stt_all), 5] >= cond)
})

test_that("A relaxed-immigration should cond run silent with correct output", {
  set.seed(3)
  replicates <- 1
  pars <- c(1, 1, 10, 0.01, 1, 0.01)
  cond <- 5
  testthat::expect_silent(
    out_cond <- DAISIE_sim_relaxed_rate(
      time = 5,
      M = 100,
      pars = pars,
      replicates = replicates,
      relaxed_par = "immigration",
      plot_sims = FALSE,
      verbose = FALSE,
      cond = cond
    )
  )

  testthat::expect_true(
    out_cond[[1]][[1]]$stt_all[nrow(out_cond[[1]][[1]]$stt_all), 5] >= cond
  )

  set.seed(3)
  testthat::expect_silent(
    out_no_cond <- DAISIE_sim_relaxed_rate(
      time = 5,
      M = 100,
      pars = pars,
      replicates = replicates,
      relaxed_par = "immigration",
      plot_sims = FALSE,
      verbose = FALSE,
      cond = 0
    )
  )
  testthat::expect_true(
    out_no_cond[[1]][[1]]$stt_all[nrow(out_no_cond[[1]][[1]]$stt_all), 5] < 5
  )

})

test_that("A relaxed-anagenesis should run silent with correct output", {
  set.seed(1)
  replicates <- 1
  testthat::expect_silent(
    sim <- DAISIE_sim_relaxed_rate(
      time = 5,
      M = 100,
      pars = c(1, 1, 10, 0.01, 5, 5),
      replicates = replicates,
      relaxed_par = "anagenesis",
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
  testthat::expect_true(is.list(sim))
  testthat::expect_true(length(sim) == replicates)
  testthat::expect_equal(sim[[1]][[1]]$island_age, 5)
  testthat::expect_equal(sim[[1]][[1]]$not_present, 98)
  testthat::expect_true(is.matrix(sim[[1]][[1]]$stt_all))
  testthat::expect_equal(nrow(sim[[1]][[1]]$stt_all), 26)
  testthat::expect_equal(ncol(sim[[1]][[1]]$stt_all), 5)
  testthat::expect_length(sim[[1]], 3)
  testthat::expect_equal(sim[[1]][[2]]$branching_times,
               c(5.00000000000000, 2.05346943810582))
  testthat::expect_equal(sim[[1]][[2]]$stac, 2)
  testthat::expect_equal(sim[[1]][[2]]$missing_species, 0)
})

test_that("A relaxed-anagenesis should cond run silent with correct output", {
  set.seed(Sys.time())
  replicates <- 1
  pars <- c(1, 1, 10, 0.01, 5, 5)
  cond <- 5
  testthat::expect_silent(
    out <- DAISIE_sim_relaxed_rate(
      time = 5,
      M = 100,
      pars = pars,
      replicates = replicates,
      relaxed_par = "anagenesis",
      plot_sims = FALSE,
      verbose = FALSE,
      cond = cond
    )
  )
  testthat::expect_true(out[[1]][[1]]$stt_all[nrow(out[[1]][[1]]$stt_all), 5] >= cond)
})

test_that("A relaxed-anagenesis should cond run silent with correct output", {
  set.seed(3)
  replicates <- 1
  pars <- c(1, 1, 10, 0.01, 5, 5)
  cond <- 5
  testthat::expect_silent(
    out_cond <- DAISIE_sim_relaxed_rate(
      time = 5,
      M = 100,
      pars = pars,
      replicates = replicates,
      relaxed_par = "anagenesis",
      plot_sims = FALSE,
      verbose = FALSE,
      cond = cond
    )
  )

  testthat::expect_true(
    out_cond[[1]][[1]]$stt_all[nrow(out_cond[[1]][[1]]$stt_all), 5] >= cond
  )

  set.seed(1)
  testthat::expect_silent(
    out_no_cond <- DAISIE_sim_relaxed_rate(
      time = 5,
      M = 100,
      pars = pars,
      replicates = replicates,
      relaxed_par = "anagenesis",
      plot_sims = FALSE,
      verbose = FALSE,
      cond = 0
    )
  )
  testthat::expect_true(
    out_no_cond[[1]][[1]]$stt_all[nrow(out_no_cond[[1]][[1]]$stt_all), 5] < 5
  )

})

test_that("Output is silent and correct for a nonoceanic simulation", {
  set.seed(1)
  replicates <- 1
  testthat::expect_silent(
    sim <- DAISIE_sim_relaxed_rate(
      time = 5,
      M = 100,
      pars = c(5, 1, 10, 0.01, 1, 5),
      replicates = replicates,
      relaxed_par = "cladogenesis",
      nonoceanic_pars = c(0.1, 0.9),
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
  testthat::expect_true(is.list(sim))
  testthat::expect_true(length(sim) == replicates)
  #number of immigrants (nonendemics) is greater than zero
  testthat::expect_gt(sim[[1]][[1]]$stt_all[1, 2], 0)
  #number of anagenetic species (endemic) is greater than zero
  testthat::expect_gt(sim[[1]][[1]]$stt_all[1, 3], 0)
})

