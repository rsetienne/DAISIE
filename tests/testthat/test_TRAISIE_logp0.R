test_that("logp0", {

  if (requireNamespace("DAISIE")) {

    if (packageVersion("DAISIE") >= "4.6.0") {

      data("Galapagos_datalist", package = "DAISIE")
      datalist <- Galapagos_datalist

      datalist[[1]]$Mainland_pool_sizes <- c(550, 250)
      datalist[[1]]$M <- 1000

      parameter <- list(2.546591, 2.678781, 0.009326754, 1.008583,
                        matrix(c(0), nrow = 1), 0, 1)

      res1 <-  DAISIE_DE_trait_logp0(
           datalist            = datalist,
           parameter           = parameter,
           num_observed_states = 1,
           num_hidden_states   = 1,
           atol                = 1e-10,
           rtol                = 1e-10,
           methode             = "lsodes"
         )

      res2 <- DAISIE:::DAISIE_DE_logp0(
                island_age = datalist[[1]]$island_age,
                pars1 = c(2.546591, 2.678781, 2.678781, 0.009326754, 1.008583),
                methode = "ode45",
                reltolint = 1e-10,
                abstolint = 1e-10)

      testthat::expect_equal(res1, res2)

      res3 <-  DAISIE_DE_trait_logp0(
        datalist            = datalist,
        parameter           = parameter,
        num_observed_states = 1,
        num_hidden_states   = 1,
        atol                = 1e-10,
        rtol                = 1e-10,
        use_Rcpp = 2
      )
      testthat::expect_equal(res1, res3)
    }
  }
})
