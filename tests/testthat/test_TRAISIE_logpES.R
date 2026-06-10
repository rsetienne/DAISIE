test_that("logpES", {

  if (requireNamespace("DAISIE")) {
    data("Galapagos_datalist", package = "DAISIE")
    datalist <- Galapagos_datalist
    datalist[[1]]$Mainland_pool_sizes <- c(550, 250)
    datalist[[1]]$M <- 1000


    i <- 7
    brts <- datalist[[i]]$branching_times
    trait <- 0
    sf <- 1

    parameter <- list(2.546591, 2.678781, 0.009326754, 1.008583,
                      matrix(c(0), nrow = 1), 0, NA)


    res1 <-  TRAISIE::DAISIE_DE_trait_logpES(
      datalist                = datalist,
      brts                    = brts,
      trait                   = trait,
      status                  = 2,
      parameter               = parameter,
      sampling_fraction       = c(1, 1),
      trait_mainland_ancestor = NA,
      num_observed_states     = 1,
      num_hidden_states       = 1,
      atol                    = 1e-10,
      rtol                    = 1e-10,
      methode                 = "ode45"
    )

    res2 <- DAISIE:::DAISIE_loglik_CS_choice(
                      pars1 = c(2.546591, 2.678781, Inf, 0.009326754, 1.008583),
                      pars2 = c(100, 11, 0, 2),
                      brts = brts,
                      stac = 2,
                      missnumspec = 0,
                      datalist = datalist)

    testthat::expect_equal(res1$loglik, res2)

    res3 <-  TRAISIE::DAISIE_DE_trait_logpES(
      datalist                = datalist,
      brts                    = brts,
      trait                   = trait,
      status                  = 2,
      parameter               = parameter,
      num_observed_states     = 1,
      num_hidden_states       = 1,
      sampling_fraction       = c(1, 1),
      methode                 = "ode45",
      use_Rcpp                = 0)

    testthat::expect_equal(res1, res3)

    res4 <-  TRAISIE::DAISIE_DE_trait_logpES(
      datalist                = datalist,
      brts                    = brts,
      trait                   = trait,
      status                  = 2,
      parameter               = parameter,
      num_observed_states     = 1,
      num_hidden_states       = 1,
      sampling_fraction       = c(1, 1),
      methode                 = "ode45",
      use_Rcpp                = 2)

    testthat::expect_equal(res1, res4)
  }
})
