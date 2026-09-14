test_that("NA values", {
  data("Galapagos_datalist", package = "DAISIE")
  datalist <- Galapagos_datalist

  i <- 7
  brts <- datalist[[i]]$branching_times
  traits <- NA
  sf <- 1

  p <- q <- matrix(c(0), nrow = 1)

  parameter <- list(2.546591, 2.678781, 0.009326754, 1.008583,
                    p, q, NA)


  res1 <-  TRAISIE_logpES(
    datalist                = datalist,
    brts                    = brts,
    traits                  = traits,
    stac                    = 2,
    parameter               = parameter,
    sampling_fraction       = rep(1, length(traits)),
    num_observed_states     = 1,
    num_hidden_states       = 1,
    atol                    = 1e-10,
    rtol                    = 1e-10,
    methode                 = "ode45"
  )
  testthat::expect_true(!is.na(res1$loglik))
})
