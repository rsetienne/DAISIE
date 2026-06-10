test_that("logpNE_max_min_age_coltime", {
    brts <- c(4, 3, 2.5)
    data("Galapagos_datalist", package = "DAISIE")
    datalist <- Galapagos_datalist
    datalist[[1]]$Mainland_pool_sizes <- c(550, 250)
    datalist[[1]]$M <- 1000

    parameter <- list(2.546591, 2.678781, 0.009326754, 1.008583,
                      matrix(c(0), nrow = 1), 0, NA)

    res1 <-  DAISIE_DE_trait_logpNE_max_min_age_hidden(
                                             datalist               = datalist,
                                              brts                  = brts,
                                              trait                 = 0,
                                              status                = 8,
                                              parameter             = parameter,
                                              trait_mainland_ancestor = NA,
                                              num_observed_states   = 1,
                                              num_hidden_states     = 1,
                                              sampling_fraction     = 1,
                                              atol                  = 1e-15,
                                              rtol                  = 1e-15,
                                              methode               = "ode45"
                                              )

    pars1 <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)
    res2  <-  DAISIE:::DAISIE_loglik_CS_choice(pars1 = pars1,
                                               pars2 = c(100, 11, 0, 2),
                                               brts = brts,
                                               stac = 8,
                                               missnumspec = 0,
                                               datalist = datalist)

    testthat::expect_equal(res1$loglik, res2, tolerance = 0.01)

    res3 <-  DAISIE_DE_trait_logpNE_max_min_age_hidden(
                                              datalist = datalist,
                                              brts                  = brts,
                                              trait                 = 0,
                                              status                = 8,
                                              parameter             = parameter,
                                              trait_mainland_ancestor = NA,
                                              num_observed_states   = 1,
                                              num_hidden_states     = 1,
                                              sampling_fraction     = 1,
                                              atol                  = 1e-15,
                                              rtol                  = 1e-15,
                                              methode               = "ode45",
                                              use_Rcpp              = 2)
    testthat::expect_equal(res1, res3, tolerance = 0.01)
  }
)
