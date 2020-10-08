context("DAISIE_ML3")
test_that("use", {
  if (Sys.getenv("TRAVIS") != "" || Sys.getenv("APPVEYOR") != "") {
    skip("WIP")
    # This is a rough MLE test, built for fast execution. A more thorough test
    # can be found in the GitHub repository Neves-P/DAISIEtesting

    utils::data(Galapagos_datalist, package = "DAISIE")
    pars1 <- c(0.2, 0.1, 17, 0.001, 0.3)
    pars1_td <- c(
      max_area = 1,
      proportional_peak_t = 0.2,
      peak_sharpness = 1,
      total_island_age = 15,
      lac = pars1[1],
      mu_min = pars1[2],
      mu_max = pars1[2],
      K0 = pars1[3],
      gam = pars1[4],
      laa = pars1[5]
    )
    tested_MLE <- DAISIE:::DAISIE_ML3(
      datalist = Galapagos_datalist,
      initparsopt = pars1_td[5:10],
      idparsopt = 5:10,
      parsfix = pars1_td[1:4],
      idparsfix = 1:4,
      island_ontogeny = 1
    )
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
    expect_equal(tested_MLE, expected_MLE)
  } else {
    testthat::skip("Run only on Travis or AppVeyor")
  }
})

test_that("abuse", {
  if (Sys.getenv("TRAVIS") != "" || Sys.getenv("APPVEYOR") != "") {
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
  } else {
    testthat::skip("Run only on Travis or AppVeyor")
  }
})
