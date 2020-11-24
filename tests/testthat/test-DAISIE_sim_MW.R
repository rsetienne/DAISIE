context("DAISIE_sim_MW")

test_that("DAISIE_sim_MW works", {

  skip(message = "Too slow to run")
  archipelago_data <- NULL
  rm(archipelago_data)
  utils::data(archipelago_data)
  replicates <- 2
  M <- 1000

  #M1 model
  pars <- c(0.025,	0.249,	1.873,	0.145,	1.15098E-13,	8.408,	51.37,	0.254,	0.047,	0.416)
  distance_dep <- 'power'
  cladogenesis_dep <- 'NULL'
  sigmoidal_par <- 'NULL'
  M1 <- DAISIE_sim_MW(
    archipelago_data = archipelago_data,
    M = M,
    pars = pars,
    replicates = replicates,
    distance_dep = distance_dep,
    cladogenesis_dep = cladogenesis_dep,
    sigmoidal_par = sigmoidal_par,
    divdepmodel = 'CS'
  )
  testthat::expect(!is.null(M1),'Empty result')

  # M14 model
  pars <- c(0.025849706,	0.25587591,	2.024864154,	0.143435441,	Inf,	0,	66.03008037,	0.285332264,	0.049967897,	0.406500832)
  distance_dep <- 'power'
  cladogenesis_dep <- 'NULL'
  sigmoidal_par <- 'NULL'
  M14 <- DAISIE_sim_MW(
    archipelago_data = archipelago_data,
    M = M,
    pars = pars,
    replicates = replicates,
    distance_dep = distance_dep,
    cladogenesis_dep = cladogenesis_dep,
    sigmoidal_par = sigmoidal_par,
    divdepmodel = 'CS'
  )
  testthat::expect(!is.null(M14),'Empty result')

  # M15 model
  pars <- c(0.008,	0.188,	1.942,	0.150,	Inf, 0, 67.65,	0.295,	0.060,	0.380, 0.231)
  distance_dep <- 'power'
  cladogenesis_dep <- 'additive'
  sigmoidal_par <- 'NULL'
  M15 <- DAISIE_sim_MW(
    archipelago_data = archipelago_data,
    M = M,
    pars = pars,
    replicates = replicates,
    distance_dep = distance_dep,
    cladogenesis_dep = cladogenesis_dep,
    sigmoidal_par = sigmoidal_par,
    divdepmodel = 'CS'
  )
  testthat::expect(!is.null(M15),'Empty result')

  # M16 model
  pars <- c( 0.040,	0.00002,	1.943,	0.150, Inf, 0,	67.50,	0.294,	0.059,	0.384,	0.026)
  distance_dep <- 'power'
  cladogenesis_dep <- 'interactive'
  sigmoidal_par <- 'NULL'
  M16 <- DAISIE_sim_MW(
    archipelago_data = archipelago_data,
    M = M,
    pars = pars,
    replicates = replicates,
    distance_dep = distance_dep,
    cladogenesis_dep = cladogenesis_dep,
    sigmoidal_par = sigmoidal_par,
    divdepmodel = 'CS'
  )
  testthat::expect(!is.null(M16),'Empty result')

  # M19 model
  pars <- c(0.040073803,	0,	1.945656546,	0.150429656,	Inf,	0,	67.25643672,	0.293635061,
            0.059096872,	0.382688527,	0.026510781)
  distance_dep <- 'power'
  cladogenesis_dep <- 'interactive'
  sigmoidal_par <- 'NULL'
  M19 <- DAISIE_sim_MW(
    archipelago_data = archipelago_data,
    M = M,
    pars = pars,
    replicates = replicates,
    distance_dep = distance_dep,
    cladogenesis_dep = cladogenesis_dep,
    sigmoidal_par = sigmoidal_par,
    divdepmodel = 'CS'
  )
  testthat::expect(!is.null(M19),'Empty result')

  # M17 model
  pars <- c(0.049,	0.116, 2.006,	0.156, Inf, 0,	67.29,	0.294,	0.057,	0.389, 42764.54)
  distance_dep <- 'power'
  cladogenesis_dep <- 'interactive1'
  sigmoidal_par <- 'NULL'
  M17 <- DAISIE_sim_MW(
    archipelago_data = archipelago_data,
    M = M,
    pars = pars,
    replicates = replicates,
    distance_dep = distance_dep,
    cladogenesis_dep = cladogenesis_dep,
    sigmoidal_par = sigmoidal_par,
    divdepmodel = 'CS'
  )
  testthat::expect(!is.null(M17),'Empty result')

  # M18 model
  pars <- c(0.044,	0.131,		1.964,	0.152, Inf, 0,	63.77,	0.286,	0.056,	0.391, 42566.78)
  distance_dep <- 'power'
  cladogenesis_dep <- 'interactive2'
  sigmoidal_par <- 'NULL'
  M18 <- DAISIE_sim_MW(
    archipelago_data = archipelago_data,
    M = M,
    pars = pars,
    replicates = replicates,
    distance_dep = distance_dep,
    cladogenesis_dep = cladogenesis_dep,
    sigmoidal_par = sigmoidal_par,
    divdepmodel = 'CS'
  )
  testthat::expect(!is.null(M18),'Empty result')

  # sigmoidal colonisation
  pars <- c(0.02,	0.25,	1.86,	0.14,	0.00000000001,	8.74,			108.35,	0.28,	0.05,	0.42,	0.17)
  distance_dep <- 'sigmoidal'
  sigmoidal_par <- 'colonisation'
  cladogenesis_dep <- 'NULL'
  M_sigm_col <- DAISIE_sim_MW(
    archipelago_data = archipelago_data,
    M = M,
    pars = pars,
    replicates = replicates,
    distance_dep = distance_dep,
    cladogenesis_dep = cladogenesis_dep,
    sigmoidal_par = sigmoidal_par,
    divdepmodel = 'CS'
  )
  testthat::expect(!is.null(M_sigm_col),'Empty result')

  # sigmoidal anagenesis
  pars <- c(0.02,	0.26,	1.85,	0.14,	3.82029371116082E-13,	8.07,	51.51,	0.25,	1.07,	1.70,	293.64)
  distance_dep <- 'sigmoidal'
  sigmoidal_par <- 'anagenesis'
  cladogenesis_dep <- 'NULL'
  M_sigm_ana <- DAISIE_sim_MW(
    archipelago_data = archipelago_data,
    M = M,
    pars = pars,
    replicates = replicates,
    distance_dep = distance_dep,
    cladogenesis_dep = cladogenesis_dep,
    sigmoidal_par = sigmoidal_par,
    divdepmodel = 'CS'
  )
  testthat::expect(!is.null(M_sigm_ana),'Empty result')

  # sigmoidal cladogenesis
  pars <- c(8.68,	0.31,	2.33,	0.17,	0.006,	1.57,	71.48,	0.30,	0.06,	0.37,2.26E+08)
  distance_dep <- 'sigmoidal'
  sigmoidal_par <- 'cladogenesis'
  cladogenesis_dep <- 'NULL'
  M_sigm_clado <- DAISIE_sim_MW(
    archipelago_data = archipelago_data,
    M = M,
    pars = pars,
    replicates = replicates,
    distance_dep = distance_dep,
    cladogenesis_dep = cladogenesis_dep,
    sigmoidal_par = sigmoidal_par,
    divdepmodel = 'CS'
  )
  testthat::expect(!is.null(M_sigm_clado),'Empty result')
})
