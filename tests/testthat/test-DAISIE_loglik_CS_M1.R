test_that("DAISIE_loglik_CS_M1 produces correct output",{

  dataset <- list(
    list(island_age = 1, not_present = 90),
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
         missing_species = 0)
  )
  out_1 <- c()
  out_2 <- c()
  for (i in 2:length(dataset)) {

    invisible(capture.output(out_1[i] <- DAISIE_loglik_CS_M1(
      brts = dataset[[i]]$branching_times,
      stac = dataset[[i]]$stac,
      missnumspec = dataset[[i]]$missing_species,
      pars1 = c(0.1, 0.1, 10.0, 0.1, 0.1),
      pars2 = c(1.0e+02, 1.1e+01, 0.0e+00, 1.0e+00, NA, 0.0e+00, 1.0e-04,
                1.0e-05,
                1.0e-07, 3.0e+03, 9.5e-01, 9.8e-01),
      verbose = FALSE
    )))
    invisible(capture.output(out_2[i] <- DAISIE_loglik_CS_M1(
      brts = dataset[[i]]$branching_times,
      stac = dataset[[i]]$stac,
      missnumspec = dataset[[i]]$missing_species,
      pars1 = c(0.5, 0.1, 10.0, 0.1, 0.1),
      pars2 = c(1.0e+02, 1.1e+01, 0.0e+00, 1.0e+00, NA, 0.0e+00, 1.0e-04,
                1.0e-05,
                1.0e-07, 3.0e+03, 9.5e-01, 9.8e-01),
      verbose = FALSE
    )))
  }
  out_1 <- out_1[-1]
  out_2 <- out_2[-1]
  expected_out_1 <- c(-2.49433738044967,
    -2.49434644748344,
    -3.09980594295824,
    -5.46225301144597,
    -5.46245593750427,
    -6.77188944406840,
    -2.49433857201992,
    -2.49442924792230,
    -3.09980813006604,
    -5.46225301144714,
    -5.46245593750544,
    -6.77188944407274
  )
  expected_out_2 <- c(
    -2.66327113437861,
    -2.66327862430581,
    -3.19372583697100,
    -5.51403419701533,
    -5.51422850883500,
    -6.79274094715768,
    -2.66327252801121,
    -2.66334743269888,
    -3.19372820911879,
    -5.51403419701655,
    -5.51422850883622,
    -6.79274094716205
  )
  testthat::expect_equal(out_1, expected_out_1)
  testthat::expect_equal(out_2, expected_out_2)

  # Max_ages at island age should be very close to max ages at very close to
  # island age
  testthat::expect_lt(out_1[1] - out_1[2], 1e-3)
  testthat::expect_lt(out_1[4] - out_1[5], 1e-3)
  testthat::expect_lt(out_1[7] - out_1[8], 1e-3)
  testthat::expect_lt(out_1[10] - out_1[11], 1e-3)
  testthat::expect_lt(out_2[1] - out_2[2], 1e-3)
  testthat::expect_lt(out_2[4] - out_2[5], 1e-3)
  testthat::expect_lt(out_2[7] - out_2[8], 1e-3)
  testthat::expect_lt(out_2[10] - out_2[11], 1e-3)
})
