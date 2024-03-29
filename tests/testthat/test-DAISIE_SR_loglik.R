test_that("The SR loglik code works", {
  Biwa_datalist <- NULL
  rm(Biwa_datalist)
  utils::data(Biwa_datalist, package = "DAISIE")
  pars1 <- c(
    0.077,
    0.956,
    Inf,
    0.138,
    0.442,
    0.077,
    0.956,
    Inf,
    0.138,
    0.442,
    0.1951
  )
  pars2 <- c(100, 11, 0, 0)
  loglik_shift <- DAISIE_SR_loglik_CS(pars1 = pars1,
                                   pars2 = pars2,
                                   datalist = Biwa_datalist,
                                   methode = 'odeint::runge_kutta_fehlberg78')
  loglik_no_shift <- DAISIE_loglik_CS(pars1 = pars1[1:5],
                                   pars2 = pars2,
                                   datalist = Biwa_datalist)
  testthat::expect_equal(loglik_shift, loglik_no_shift, tol = 1E-3)
  testthat::expect_equal(loglik_no_shift, -245.6298, tol = 1E-3)

  Macaronesia_datalist <- NULL
  rm(Macaronesia_datalist)
  utils::data(Macaronesia_datalist, package = "DAISIE")
  pars1 <- c(
    0.1,
    1.1,
    10,
    0.6,
    0.05,
    0.1,
    1.1,
    10,
    0.6,
    0.05,
    7
  )
  pars1mat <- matrix(pars1, nrow = 8, ncol = 11, byrow = T)
  expected_loglik <- c(
    -Inf,
    -252.7293,
    -Inf,
    -246.9084,
    -Inf,
    -272.9101,
    -Inf,
    -247.1446
  )
  # No shift in cladogenesis rate older before
  # colonization of diversifying lineages
  pars1mat[1, 1] <- 0.01; pars1mat[1, 6] <- 0.3
  pars1mat[2, 1] <- 0.01; pars1mat[2, 6] <- 0.3; pars1mat[2, 11] <- 2.1
  # No shift in extinction rate older than known ages
  pars1mat[3, 2] <- 0.1; pars1mat[3, 11] <- 12
  pars1mat[4, 2] <- 0.1;
  # No shift in colonization rate older than known ages
  pars1mat[5, 4] <- 0.1; pars1mat[5, 9] <- 1.5; pars1mat[5, 11] <- 10
  pars1mat[6, 4] <- 0.1; pars1mat[6, 9] <- 0.8; pars1mat[6, 11] <- 7
  # No shift in anagenetic rate older than any known non-endemic
  pars1mat[7, 5] <- 0.01; pars1mat[7, 10] <- 0.1;
  pars1mat[8, 5] <- 0.01; pars1mat[8, 10] <- 0.1; pars1mat[8, 11] <- 2.9
  loglik <- rep(0,8)
  for (i in 1:8) {
    loglik[i] <- DAISIE_SR_loglik_CS(
      pars1 = pars1mat[i, -12],
      pars2 = pars2,
      datalist = Macaronesia_datalist[[2]],
      methode = 'odeint::runge_kutta_fehlberg78'
    )
  }
  testthat::expect_equal(loglik, expected_loglik, tol = 1E-3)
})
