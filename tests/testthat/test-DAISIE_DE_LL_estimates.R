test_that("DAISIE_logp0 is correct", {
  data("Galapagos_datalist", package = "DAISIE")
  datalist <- Galapagos_datalist

  parameter <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)

  res1 <- DAISIE:::DAISIE_DE_logp0(
    island_age = datalist[[1]]$island_age,
    pars1 = c(2.546591, 2.678781, 2.678781, 0.009326754, 1.008583),
    methode = "lsodes",
    reltolint = 1e-15,
    abstolint = 1e-15,
    use_rcpp = FALSE
  )

  res2 <- DAISIE:::DAISIE_DE_logp0(island_age = datalist[[1]]$island_age,
                                   parameter,
                                   abstolint = 1e-15,
                                   reltolint = 1e-15,
                                   methode = "ode45",
                                   use_rcpp = FALSE
  )

  testthat::expect_equal(res1, res2, tolerance = 1e-6)

  res3 <- DAISIE:::DAISIE_DE_logp0(island_age = datalist[[1]]$island_age,
                                   parameter,
                                   abstolint = 1e-15,
                                   reltolint = 1e-15,
                                   use_rcpp = TRUE
  )

  testthat::expect_equal(res1, res3, tolerance = 1e-6)
})

test_that("logpEC", {

  data("Galapagos_datalist", package = "DAISIE")
  datalist <- Galapagos_datalist

  i <- 4
  brts <- datalist[[i]]$branching_times


  parameter <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)

  res1 <- DAISIE:::DAISIE_DE_logpEC(brts,
                                    missnumspec = 0,
                                    pars1 = parameter,
                                    abstolint  = 1e-15,
                                    reltolint  = 1e-15,
                                    methode                 = "ode45",
                                    use_rcpp = FALSE)

  res2 <- DAISIE:::DAISIE_loglik_CS_choice(pars1 = c(2.546591, 2.678781, Inf, 0.009326754, 1.008583),
                                           pars2 = c(100, 11, 0, 2),
                                           brts = brts,
                                           stac = 2,
                                           missnumspec = 0,
                                           datalist = datalist)
  testthat::expect_equal(res1, res2)

  res3 <- DAISIE:::DAISIE_DE_logpEC(brts,
                                    missnumspec = 0,
                                    pars1 = parameter,
                                    abstolint  = 1e-15,
                                    reltolint  = 1e-15,
                                    methode                 = "ode45",
                                    use_rcpp = TRUE)
  testthat::expect_equal(res3, res1)
})

test_that("logpES", {

  data("Galapagos_datalist", package = "DAISIE")
  datalist <- Galapagos_datalist

  i <- 6
  brts <- datalist[[i]]$branching_times


  parameter <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)

  res1 <- DAISIE:::DAISIE_DE_logpES(brts,
                                    missnumspec = 0,
                                    pars1 = parameter,
                                    abstolint  = 1e-15,
                                    reltolint  = 1e-15,
                                    methode = "ode45",
                                    use_rcpp = FALSE)

  res2 <- DAISIE:::DAISIE_loglik_CS_choice(pars1 = c(2.546591, 2.678781, Inf, 0.009326754, 1.008583),
                                           pars2 = c(100, 11, 0, 2),
                                           brts = brts,
                                           stac = 2,
                                           missnumspec = 0,
                                           datalist = datalist)

  testthat::expect_equal(res1, res2)

  res3 <- DAISIE:::DAISIE_DE_logpES(brts,
                                    missnumspec = 0,
                                    pars1 = parameter,
                                    abstolint  = 1e-15,
                                    reltolint  = 1e-15,
                                    methode = "ode45",
                                    use_rcpp = TRUE)
  testthat::expect_equal(res3, res1)
})

test_that("logpNE", {

  data("Galapagos_datalist", package = "DAISIE")
  datalist <- Galapagos_datalist

  i <- 3
  brts <- datalist[[i]]$branching_times


  parameter <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)

  res1 <- DAISIE:::DAISIE_DE_logpNE(brts,
                                    pars1 = parameter,
                                    abstolint  = 1e-15,
                                    reltolint  = 1e-15,
                                    methode = "ode45",
                                    use_rcpp = FALSE)

  res2 <- DAISIE:::DAISIE_loglik_CS_choice(pars1 = parameter,
                                           pars2 = c(100, 11, 0, 2),
                                           brts = brts,
                                           stac = 4,
                                           missnumspec = 0,
                                           datalist = datalist)
  # TODO: FAILS
  testthat::expect_equal(res1, res2)

  res2 <- DAISIE:::DAISIE_DE_logpNE(brts,
                                    pars1 = parameter,
                                    abstolint  = 1e-15,
                                    reltolint  = 1e-15,
                                    methode = "ode45",
                                    use_rcpp = TRUE)
  testthat::expect_equal(res1, res2)
})

test_that("logpES_max_min_age_coltime", {

  data("Galapagos_datalist", package = "DAISIE")
  datalist <- Galapagos_datalist


  brts <- c(8, 5, 3)


  parameter <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)

  res1 <- DAISIE:::DAISIE_DE_logpES_max_min_age_coltime(brts,
                                                        pars1 = parameter,
                                                        missnumspec = 0,
                                                        abstolint  = 1e-15,
                                                        reltolint  = 1e-15,
                                                        methode                 = "ode45",
                                                        use_rcpp = FALSE)

  res2 <- DAISIE:::DAISIE_loglik_CS_choice(pars1 = parameter,
                                           pars2 = c(100, 11, 0, 2),
                                           brts = brts,
                                           stac = 9,
                                           missnumspec = 0,
                                           datalist = datalist)
  testthat::expect_equal(res1, res2, tolerance = 1e-2)

  res3 <- DAISIE:::DAISIE_DE_logpES_max_min_age_coltime(brts,
                                                        pars1 = parameter,
                                                        missnumspec = 0,
                                                        abstolint  = 1e-15,
                                                        reltolint  = 1e-15,
                                                        methode                 = "ode45",
                                                        use_rcpp = TRUE)
  testthat::expect_equal(res3, res1, tolerance = 1e-2)
})

test_that("logpNE_max_min_age_coltime", {

  data("Galapagos_datalist", package = "DAISIE")
  datalist <- Galapagos_datalist


  brts <- c(5, 4, 3)


  parameter <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)

  res1 <- DAISIE:::DAISIE_DE_logpNE_max_min_age_coltime(brts,
                                                        pars1 = parameter,
                                                        abstolint  = 1e-15,
                                                        reltolint  = 1e-15,
                                                        methode                 = "ode45",
                                                        use_rcpp = FALSE)

  res2 <- DAISIE:::DAISIE_loglik_CS_choice(pars1 = c(2.546591, 2.678781, Inf, 0.009326754, 1.008583),
                                           pars2 = c(100, 11, 0, 2),
                                           brts = brts,
                                           stac = 8,
                                           missnumspec = 0,
                                           datalist = datalist)

  testthat::expect_equal(res1, res2, tolerance = 1e-2)

  res3 <- DAISIE:::DAISIE_DE_logpNE_max_min_age_coltime(brts,
                                                        pars1 = parameter,
                                                        abstolint  = 1e-15,
                                                        reltolint  = 1e-15,
                                                        methode                 = "ode45",
                                                        use_rcpp = TRUE)
  testthat::expect_equal(res3, res1, tolerance = 1e-2)
})

test_that("logpEC general", {
  data("Galapagos_datalist", package = "DAISIE")
  datalist <- Galapagos_datalist

  i <- 4
  brts <- datalist[[i]]$branching_times


  parameter <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)

  res1 <- DAISIE:::DAISIE_DE_logpEC(brts,
                                    missnumspec = 0,
                                    pars1 = parameter,
                                    abstolint  = 1e-15,
                                    reltolint  = 1e-15,
                                    methode                 = "ode45",
                                    use_rcpp = FALSE)

  res2 <- DAISIE:::DAISIE_DE_logpEC_general(brts = brts,
                                            missnumspec = 0,
                                            mainland = FALSE,
                                            coltime = "not_chosen",
                                            pars1 = parameter,
                                            abstolint  = 1e-15,
                                            reltolint  = 1e-15,
                                            methode                 = "ode45",
                                            use_rcpp = FALSE)
  testthat::expect_equal(res1, res2)

  res1 <- DAISIE:::DAISIE_DE_logpEC_mainland(brts,
                                             missnumspec = 0,
                                             pars1 = parameter,
                                             abstolint  = 1e-15,
                                             reltolint  = 1e-15,
                                             methode                 = "ode45",
                                             use_rcpp = FALSE)
  res2 <- DAISIE:::DAISIE_DE_logpEC_general(brts = brts,
                                            missnumspec = 0,
                                            mainland = TRUE,
                                            coltime = "not_chosen",
                                            pars1 = parameter,
                                            abstolint  = 1e-15,
                                            reltolint  = 1e-15,
                                            methode                 = "ode45",
                                            use_rcpp = FALSE)
  testthat::expect_equal(res1, res2)

  res1 <- DAISIE:::DAISIE_DE_logpEC_unknown_coltime(brts,
                                                    missnumspec = 0,
                                                    pars1 = parameter,
                                                    abstolint  = 1e-15,
                                                    reltolint  = 1e-15,
                                                    methode                 = "ode45",
                                                    use_rcpp = FALSE)
  res2 <- DAISIE:::DAISIE_DE_logpEC_general(brts = brts,
                                            missnumspec = 0,
                                            coltime = "unknown",
                                            pars1 = parameter,
                                            abstolint  = 1e-15,
                                            reltolint  = 1e-15,
                                            methode                 = "ode45",
                                            use_rcpp = FALSE)
  testthat::expect_equal(res1, res2)

  res1 <- DAISIE:::DAISIE_DE_logpEC_max_age_coltime(brts,
                                                    missnumspec = 0,
                                                    pars1 = parameter,
                                                    abstolint  = 1e-15,
                                                    reltolint  = 1e-15,
                                                    methode                 = "ode45",
                                                    use_rcpp = FALSE)
  res2 <- DAISIE:::DAISIE_DE_logpEC_general(brts = brts,
                                            missnumspec = 0,
                                            coltime = "max_age",
                                            pars1 = parameter,
                                            abstolint  = 1e-15,
                                            reltolint  = 1e-15,
                                            methode                 = "ode45",
                                            use_rcpp = FALSE)
  testthat::expect_equal(res1, res2)

  res1 <- DAISIE:::DAISIE_DE_logpEC_max_age_coltime_and_mainland(brts,
                                                                 missnumspec = 0,
                                                                 pars1 = parameter,
                                                                 abstolint  = 1e-15,
                                                                 reltolint  = 1e-15,
                                                                 methode                 = "ode45",
                                                                 use_rcpp = FALSE)
  res2 <- DAISIE:::DAISIE_DE_logpEC_general(brts = brts,
                                            missnumspec = 0,
                                            coltime = "max_age",
                                            mainland = TRUE,
                                            pars1 = parameter,
                                            abstolint  = 1e-15,
                                            reltolint  = 1e-15,
                                            methode                 = "ode45",
                                            use_rcpp = FALSE)
  testthat::expect_equal(res1, res2)
})

test_that("logpNE general", {
  data("Galapagos_datalist", package = "DAISIE")
  datalist <- Galapagos_datalist

  i <- 3
  brts <- datalist[[i]]$branching_times
  parameter <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)

  res1 <- DAISIE:::DAISIE_DE_logpNE(brts,
                                    pars1 = parameter,
                                    abstolint  = 1e-15,
                                    reltolint  = 1e-15,
                                    methode = "ode45",
                                    use_rcpp = FALSE)
  res2 <- DAISIE:::DAISIE_DE_logpNE_general(brts,
                                            pars1 = parameter,
                                            coltime = "known",
                                            abstolint  = 1e-15,
                                            reltolint  = 1e-15,
                                            methode = "ode45",
                                            use_rcpp = FALSE)
  testthat::expect_equal(res1, res2)

  res1 <- DAISIE:::DAISIE_DE_logpNE_max_age_coltime(brts,
                                                    pars1 = parameter,
                                                    abstolint  = 1e-15,
                                                    reltolint  = 1e-15,
                                                    methode = "ode45",
                                                    use_rcpp = FALSE)
  res2 <- DAISIE:::DAISIE_DE_logpNE_general(brts,
                                            pars1 = parameter,
                                            coltime = "max_age",
                                            abstolint  = 1e-15,
                                            reltolint  = 1e-15,
                                            methode = "ode45",
                                            use_rcpp = FALSE)
  testthat::expect_equal(res1, res2)

  res1 <- DAISIE:::DAISIE_DE_logpNE_unknown_coltime(brts,
                                                    pars1 = parameter,
                                                    abstolint  = 1e-15,
                                                    reltolint  = 1e-15,
                                                    methode = "ode45",
                                                    use_rcpp = FALSE)
  res2 <- DAISIE:::DAISIE_DE_logpNE_general(brts,
                                            pars1 = parameter,
                                            coltime = "unknown",
                                            abstolint  = 1e-15,
                                            reltolint  = 1e-15,
                                            methode = "ode45",
                                            use_rcpp = FALSE)
  testthat::expect_equal(res1, res2)

  brts <- c(5, 4, 3)
  parameter <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)

  res1 <- DAISIE:::DAISIE_DE_logpNE_max_min_age_coltime(brts,
                                                        pars1 = parameter,
                                                        abstolint  = 1e-15,
                                                        reltolint  = 1e-15,
                                                        methode = "ode45",
                                                        use_rcpp = FALSE)
  res2 <- DAISIE:::DAISIE_DE_logpNE_general(brts,
                                            pars1 = parameter,
                                            coltime = "max_min_age",
                                            abstolint  = 1e-15,
                                            reltolint  = 1e-15,
                                            methode = "ode45",
                                            use_rcpp = FALSE)
  testthat::expect_equal(res1, res2)
})
