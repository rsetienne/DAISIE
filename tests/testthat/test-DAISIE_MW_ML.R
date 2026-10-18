test_that("DAISIE_MW_ML produces correct output; fast run", {
  skip_if(Sys.getenv("CI") == "", message = "Run only on CI")
  skip_on_cran()
  skip_on_covr()
  # This test uses very lax tolerances and parameters to run in
  # a short enough time. The following _skipped_ test has the
  # true expected parameters. Run it only locally.
  utils::data(archipelagos41)
  test_data <- archipelagos41[1:6]
  start <- Sys.time()
  invisible(capture.output(M19_tested <- DAISIE_MW_ML(
    datalist = test_data,
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
    res = 10,
    ddmodel = 0,
    methode = 'odeint::runge_kutta_cash_karp54',
    cpus = 1,
    parallel = 'no',
    optimmethod = 'simplex',
    tol = c(1E-1, 1E-2, 1E-3),
    distance_type = 'continent',
    distance_dep = 'area_interactive_clado',
    maxiter = 1000
  )))
  end <- Sys.time()
  M19_Nature_expected <- c(
    0.040073803,
    0,
    1.945656546,
    0.150429656,
    Inf,
    0,
    67.25643672000005,
    0.293635061,
    0.05909687199999999,
    0.382688527,
    0.026510781,
    -3651.733505146444,
    8,
    0)
  M19_expected_parameters <- c(
    0.010826711434425,
    0.0,
    6.7307303049659710,
    0.2122290735744987,
    Inf,
    0.000000e+00,
    447.2304110855639010,
    0.5018081980501639,
    0.1177972324294647,
    0.3010264314360413,
    0.0498243725674074,
    -567.4429002593708447,
    8.000000e+00,
    0.000000e+00
  )
  testthat::expect_equal(
    M19_expected_parameters,
    as.numeric(M19_tested)
  )
})

test_that("DAISIE_MW_ML produces correct output", {
  skip("Takes too long to run; run it only locally")
  utils::data(archipelagos41)
  invisible(capture.output(M19_tested <- DAISIE_MW_ML(
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
    methode = 'odeint::runge_kutta_cash_karp54',
    cpus = 4,
    parallel = 'no',
    optimmethod = 'simplex',
    tol = c(1E-1, 1E-3, 1E-5),
    distance_type = 'continent',
    distance_dep = 'area_interactive_clado'
  )))

  M19_Nature_expected <- c(
    0.040073803,
    0,
    1.945656546,
    0.150429656,
    Inf,
    0,
    67.25643672000005,
    0.293635061,
    0.05909687199999999,
    0.382688527,
    0.026510781,
    -3651.733505146444,
    8,
    0)
  M19_expected_parameters <- c(
    0.039701416238015007,
    0.000000e+00,
    1.9620274492764262,
    0.1507885411238431,
    Inf,
    0.000000e+00,
    70.330803377621763,
    0.29869325763892701,
    0.056674928493858981,
    0.38783023606642458,
    0.026694014758619697,
    -3651.7196509022697,
    8.000000e+00,
    0.000000e+00
  )
  testthat::expect_equal(
    M19_expected_parameters,
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
  M19_expected_parameters <- c(
    0.039701416238015007,
    0.000000e+00,
    1.9620274492764262,
    0.1507885411238431,
    Inf,
    0.000000e+00,
    70.330803377621763,
    0.29869325763892701,
    0.056674928493858981,
    0.38783023606642458,
    0.026694014758619697,
    -3651.7196509022697,
    8.000000e+00,
    0.000000e+00
  )
  invisible(capture.output(
    M19_computation <- DAISIE_MW_ML(
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
      methode = 'odeint::runge_kutta_cash_karp',
      cpus = 4,
      parallel = 'local',
      optimmethod = 'simplex',
      tol = c(1E-1, 1E-3, 1E-5),
      distance_type = 'continent',
      distance_dep = 'area_interactive_clado'
    )
  ))
  testthat::expect_equal(
    M19_expected_parameters,
    as.numeric(M19_computation),
    tol = 1E-5
  )
})

test_that("colonisation decreases with distance in sigmoidal_col and sigmoidal_col_ana", {
  # gamma_0 (k of the sigmoid) = 2, alpha = 3, d0 for colonisation = 100
  pars_user <- c(0.5, 0.1, 0.2, 0.1, 10, 0.1, 2, 3, 0.4, 0.1, 100, 100)
  colonisation_at <- function(distance_dep, distance) {
    n_pars <- if (distance_dep == "sigmoidal_col_ana") 12 else 11
    pars1 <- sign_conventions_MW(pars_user[1:n_pars], distance_dep)
    convert_parameters_MW(
      pars1 = pars1,
      area = 1,
      distance = distance,
      M = 1,
      distance_dep = distance_dep
    )[4]
  }
  distances <- c(10, 100, 1000)
  col <- sapply(distances, function(d) colonisation_at("sigmoidal_col", d))
  col_ana <- sapply(distances, function(d) colonisation_at("sigmoidal_col_ana", d))
  testthat::expect_true(all(diff(col) < 0))
  testthat::expect_true(all(diff(col_ana) < 0))
  # same colonisation sigmoid in both models
  testthat::expect_equal(col_ana, col)
})

test_that("sign conventions: alpha is a negative exponent in power model only", {
  pars1 <- c(0.5, 0.1, 0.2, 0.1, 10, 0.1, 2, 3, 0.4, 0.1)
  testthat::expect_equal(sign_conventions_MW(pars1, "power")[c(4, 8)], -pars1[c(4, 8)])
  pars1_sigmoidal <- c(pars1, 100)
  testthat::expect_equal(sign_conventions_MW(pars1_sigmoidal, "sigmoidal_col")[c(4, 8)], c(-0.1, 3))
  testthat::expect_equal(sign_conventions_MW(pars1_sigmoidal, "sigmoidal_col_ana")[c(4, 8)], c(-0.1, 3))
  testthat::expect_equal(sign_conventions_MW(pars1_sigmoidal, "sigmoidal_ana")[c(4, 8)], c(-0.1, -3))
})
