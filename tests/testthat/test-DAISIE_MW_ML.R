test_that("DAISIE_MW_ML produces correct output", {

  skip_if(Sys.getenv("CI") == "", message = "Run only on CI")

  utils::data(archipelagos41)

  invisible(capture.output(
    M19_tested <- DAISIE::DAISIE_MW_ML(
      datalist = archipelagos41,
      initparsopt = c(
        0.040073803,
        1.945656546,
        0.150429656,
        67.25643672,
        0.293635061,
        0.059096872,
        0.382688527,
        0.026510781),
      idparsopt = c(1, 3, 4, 7, 8, 9, 10, 11),
      parsfix = c(0, Inf, 0) ,
      idparsfix = c(2, 5, 6),
      res = 100,
      ddmodel = 0,
      methode = 'lsodes',
      cpus = 4,
      parallel = 'no',
      optimmethod = 'subplex',
      tol = c(1E-1, 1E-3, 1E-5),
      distance_type = 'continent',
      distance_dep = 'area_interactive_clado'
    )
  ))


  M19_Nature_expected <- c(
    0.040073803,
    0.0,
    1.945656546,
    0.150429656,
    Inf,
    0.0,
    67.2564367200001,
    0.293635061,
    0.059096872,
    0.382688527,
    0.026510781,
    -3651.48307905794,
    8,
    0
  )

  testthat::expect_equal(
    M19_Nature_expected,
    as.numeric(M19_tested)
  )

})

test_that("DAISIE_MW_ML produces correct output when in parallel", {
  skip('Parallel code cannot be tested automatically. Run it manually.')

  utils::data(archipelagos41)

  M19_Nature_parameters <- c(
    4.007380e-02,
    0.000000e+00,
    1.945657e+00,
    1.504297e-01,
    Inf,
    0.000000e+00,
    67.2564367200001,
    2.936351e-01,
    5.909687e-02,
    3.826885e-01,
    2.651078e-02,
    -3651.6624531664,
    8.000000e+00,
    0.000000e+00
  )

  invisible(capture.output(
    M19_computation <- DAISIE::DAISIE_MW_ML(
      datalist = archipelagos41,
      initparsopt = c(
        0.040073803,
        1.945656546,
        0.150429656,
        67.25643672,
        0.293635061,
        0.059096872,
        0.382688527,
        0.026510781),
      idparsopt = c(1, 3, 4, 7, 8, 9, 10, 11),
      parsfix = c(0, Inf, 0) ,
      idparsfix = c(2, 5, 6),
      res = 100,
      ddmodel = 0,
      methode = 'lsodes',
      cpus = 4,
      parallel = 'local',
      optimmethod = 'subplex',
      tol = c(1E-1, 1E-3, 1E-5),
      distance_type = 'continent',
      distance_dep = 'area_interactive_clado'
    )
  ))

  testthat::expect_equal(
    M19_Nature_parameters,
    as.numeric(M19_computation),
    tol = 0.000001
  )

})

