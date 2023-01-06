test_that("use", {
  skip_if(Sys.getenv("CI") == "", message = "Run only on CI")
  # THIS FUNCTION DOESN'T WORK CORRECTLY YET! FOR NOW, WE TEST IT THROWS AN
  # APPROPRIATE ERROR

  utils::data(Galapagos_datalist, package = "DAISIE")
  lac0 <- 1.000
  mu0 <- 0.400
  K0 <- 20.000
  gam0 <- 0.009
  laa0 <- 1.010
  d <- 0
  x <- 0
  ka <- 0

  area_pars <- c(
    max_area = 10,
    current_area = 1,
    proportional_peak_t = 0.5,
    total_island_age = 5,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0
  )
  island_ontogeny <- 1
  sea_level <- 0

  pars1_time_dep <- c(
    lac0,
    mu0,
    K0,
    gam0,
    laa0,
    d,
    x,
    ka,
    area_pars
  )

  methode <- "ode45"
  optimmethod <- "simplex"
  ddmodel <- 11
  cond <- 0

  time_dep_mle <- DAISIE_ML3(
    datalist = Galapagos_datalist,
    initparsopt = pars1_time_dep[1:5],
    idparsopt = 1:5,
    parsfix = pars1_time_dep[6:15],
    idparsfix = 6:15,
    island_ontogeny = 1,
    sea_level = 0,
    CS_version = 1,
    cond = cond,
    methode = methode,
    ddmodel = ddmodel,
    optimmethod = optimmethod, verbose = TRUE
  )
  constant_mle <- DAISIE_ML1(
    datalist = Galapagos_datalist,
    initparsopt = pars1_time_dep[1:5],
    idparsopt = 1:5,
    parsfix = NULL,
    idparsfix = NULL,
    CS_version = 1,
    cond = cond,
    methode = methode,
    ddmodel = ddmodel,
    optimmethod = optimmethod, verbose = TRUE
  )

  # All code below refers to future reference test when function is completed
  idpars <- sort(c(5:10, 1:4))
  expected_MLE <- data.frame(
    lambda_c = c(0.0,
                 0.133766934,
                 0.0,
                 0.0),
    mu = c(1.0531518319999997,
           1.0531518319999997,
           1.0531518319999997,
           1.0531518319999997),
    K = c(Inf,
          Inf,
          Inf,
          Inf),
    gamma = c(0.052148978999999998,
              0.152763178999999999,
              0.052148978999999998,
              0.052148978999999998),
    lambda_a = c(0.51293901099999994,
                 0.51293901099999994,
                 0.51293901099999994,
                 0.51293901099999994),
    loglik = c(-454.93478332906614,
               -454.93478332906614,
               -454.93478332906614,
               -454.93478332906614),
    df = c(5L, 5L, 5L, 5L),
    conv = c(0L, 0L, 0L, 0L)
  )

})

test_that("abuse", {
  skip_if(Sys.getenv("CI") == "", message = "Run only on CI")
  expect_error(tested_MLE <- DAISIE:::DAISIE_ML2(
    datalist = "nonsense",
    initparsopt = c(
      1.053151832,
      0.052148979,
      0.512939011,
      0.133766934,
      0.152763179
    ),
    idparsmat = rbind(
      1:5,
      c(6, 2, 3, 7, 5),
      1:5,1:5
    ),
    idparsopt = c(2, 4, 5, 6, 7),
    parsfix = c(0, Inf),
    idparsfix = c(1, 3),
    tol = c(0.01, 0.1, 0.001),
    res = 15,
    tolint = c(0.1, 0.01)
  ))
})
