test_that("DAISIE_loglik_CS_M1 produces correct output",{

  dataset <- list(
    list(island_age = 1, not_present = 90),
    list(branching_times = c(1.0, 0.99999999), stac = 1, missing_species = 0),
    list(branching_times = c(1.0, 0.99999), stac = 1, missing_species = 0),
    list(branching_times = c(1.0, 0.999989), stac = 1, missing_species = 0),
    list(branching_times = c(1.0, 0.5), stac = 1, missing_species = 0),
    list(branching_times = c(1.0, 0.99999), stac = 5, missing_species = 0),
    list(branching_times = c(1.0, 0.99989), stac = 5, missing_species = 0),
    list(branching_times = c(1.0, 0.5), stac = 5, missing_species = 0),
    list(branching_times = c(1.0, 0.99999, 0.000001), stac = 8,
         missing_species = 0),
    list(branching_times = c(1.0, 0.99989, 0.000001), stac = 8,
         missing_species = 0),
    list(branching_times = c(1.0, 0.5, 0.000001), stac = 8,
         missing_species = 0),
    list(branching_times = c(1.0, 0.99999, 0.000001), stac = 9,
         missing_species = 0),
    list(branching_times = c(1.0, 0.99989, 0.000001), stac = 9,
         missing_species = 0),
    list(branching_times = c(1.0, 0.5, 0.000001), stac = 9,
         missing_species = 0),
    list(branching_times = c(1.0, 0.5, 0.1), stac = 6,
         missing_species = 0),
    list(branching_times = c(1.0, 0.5, 0.1), stac = 7,
         missing_species = 0)
  )
  out_1 <- c()
  out_2 <- c()
  for (i in 2:length(dataset)) {

    out_1[i - 1] <- DAISIE_loglik_CS_M1(
      brts = dataset[[i]]$branching_times,
      stac = dataset[[i]]$stac,
      missnumspec = dataset[[i]]$missing_species,
      pars1 = c(0.1, 0.1, 10.0, 0.1, 0.1),
      pars2 = c(100, 11, 0, 1),
      verbose = FALSE,
      CS_version = list(model = 1, function_to_optimize = 'DAISIE')
    )
    out_2[i - 1] <- DAISIE_loglik_CS_M1(
      brts = dataset[[i]]$branching_times,
      stac = dataset[[i]]$stac,
      missnumspec = dataset[[i]]$missing_species,
      pars1 = c(0.5, 0.1, 10.0, 0.1, 0.1),
      pars2 = c(100, 11, 0, 1),
      verbose = FALSE,
      CS_version = list(model = 1, function_to_optimize = 'DAISIE')
    )
  }
  expected_out_1 <- c(-2.494337,
                      -2.494337,
                      -2.494513,
                      -3.100835,
                      -5.462253,
                      -5.462688,
                      -6.774352,
                      -2.494339,
                      -2.494596,
                      -3.100837,
                      -5.462253,
                      -5.462688,
                      -6.774352,
                      -5.630455,
                      -10.441115)
  expected_out_2 <- c(-2.663271,
                      -2.663271,
                      -2.663444,
                      -3.194598,
                      -5.514034,
                      -5.514461,
                      -6.795040,
                      -2.663273,
                      -2.663513,
                      -3.194601,
                      -5.514034,
                      -5.514461,
                      -6.795040,
                      -4.158192,
                      -8.976662)
  testthat::expect_equal(out_1, expected_out_1, 1e-5)
  testthat::expect_equal(out_2, expected_out_2, 1e-5)
  # Max_ages at island age should be very close to max ages at very close to
  # island age
  testthat::expect_lt(out_1[1] - out_1[2], 1e-3)
  testthat::expect_lt(out_1[2] - out_1[3], 1e-3)
  testthat::expect_lt(out_1[5] - out_1[6], 1e-3)
  testthat::expect_lt(out_1[8] - out_1[9], 1e-3)
  testthat::expect_lt(out_1[11] - out_1[12], 1e-3)
  testthat::expect_lt(out_2[1] - out_2[2], 1e-3)
  testthat::expect_lt(out_2[2] - out_2[3], 1e-3)
  testthat::expect_lt(out_2[5] - out_2[6], 1e-3)
  testthat::expect_lt(out_2[8] - out_2[9], 1e-3)
  testthat::expect_lt(out_2[11] - out_2[12], 1e-3)


  #DE should give the same result as DAISIE if max age very close to island age
  loglik_DE <- DAISIE_DE_loglik_CS(
    pars1 = c(0.1, 0.09, 0.09, 0.1, 0.1),
    pars2 = c(200, 11, 0, 1),
    datalist = dataset[c(1,2)],
    methode = 'ode45',
    reltolint = 1E-16,
    abstolint = 1E-16)
  loglik_DAISIE <- DAISIE_loglik_CS(
    pars1 = c(0.1, 0.09, Inf, 0.1, 0.1),
    pars2 = c(200, 11, 0, 1),
    datalist = dataset[c(1,2)],
    methode = 'ode45',
    reltolint = 1E-16,
    abstolint = 1E-16)
  testthat::expect_equal(loglik_DE,loglik_DAISIE,1E-6)
  #DE should also give the same result when max age is quite different
  loglik_DE <- DAISIE_DE_loglik_CS(
    pars1 = c(0.1, 0.09, 0.09, 0.1, 0.1),
    pars2 = c(200, 11, 0, 1),
    datalist = dataset[c(1,5)],
    methode = 'ode45',
    reltolint = 1E-16,
    abstolint = 1E-16)
  loglik_DAISIE <- DAISIE_loglik_CS(
    pars1 = c(0.1, 0.09, Inf, 0.1, 0.1),
    pars2 = c(200, 11, 0, 1),
    datalist = dataset[c(1,5)],
    methode = 'ode45',
    reltolint = 1E-16,
    abstolint = 1E-16)
  testthat::expect_equal(loglik_DE,loglik_DAISIE,1E-6)
  loglik_DAISIE_approx <- DAISIE_loglik_CS(
    pars1 = c(0.1, 0.09, Inf, 0.1, 0.1),
    pars2 = c(200, 11, 0, 1),
    datalist = dataset[c(1,5)],
    methode = 'ode45',
    reltolint = 1E-16,
    abstolint = 1E-16,
    CS_version = list(1,function_to_optimize = 'DAISIE_approx'))
  testthat::expect_equal(loglik_DAISIE_approx,loglik_DAISIE,1E-6)
})

