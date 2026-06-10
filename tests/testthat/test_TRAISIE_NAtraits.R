test_that("NA values", {

  if (requireNamespace("DAISIE")) {
    data("Galapagos_datalist", package = "DAISIE")
    datalist <- Galapagos_datalist

    i <- 7
    brts <- datalist[[i]]$branching_times
    trait <- NA
    sf <- 1

    parameter <- list(2.546591, 2.678781, 0.009326754, 1.008583,
                      matrix(c(0), nrow = 1), 0, NA)


    res1 <-  TRAISIE::DAISIE_DE_trait_logpES(
      datalist                = datalist,
      brts                    = brts,
      trait                   = trait,
      status                  = 2,
      parameter               = parameter,
      sampling_fraction       = rep(1, length(trait)),
      num_observed_states     = 1,
      num_hidden_states       = 1,
      atol                    = 1e-10,
      rtol                    = 1e-10,
      methode                 = "ode45"
    )
    testthat::expect_true(!is.na(res1$loglik))
  }
})
