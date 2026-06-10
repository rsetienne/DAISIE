test_that("CS", {

  if (requireNamespace("DAISIE")) {
      data("Galapagos_datalist", package = "DAISIE")
      datalist <- Galapagos_datalist

      parameter <- list(c(2.546591, 2.546591),
                        c(2.678781, 2.678781),
                        c(0.009326754, 0.009326754),
                        c(1.008583, 1.008583),
                        matrix(c(0.001), nrow = 2, ncol = 2),
                        0)

      datalist[[1]]$M0 <- datalist[[1]]$not_present / 2
      datalist[[1]]$M1 <- datalist[[1]]$not_present / 2

      for (i in 2:length(datalist)) {
        datalist[[i]]$phylogeny <-
              DDD::brts2phylo(datalist[[i]]$branching_times[-c(1, 2)])
        datalist[[i]]$phylogeny$root.edge <- 0
        num_tips <- length(datalist[[i]]$phylogeny$tip.label)
        datalist[[i]]$sampling_fraction <- rep(1, 2)
        datalist[[i]]$traits <- sample(c(0, 1), size = num_tips,
                                       replace = TRUE)
        datalist[[i]]$root_state <- c(0.5, 0.5)
      }

      res1 <-  TRAISIE::DAISIE_DE_trait_loglik_CS(
        datalist            = datalist,
        parameter           = parameter,
        num_observed_states = 2,
        num_hidden_states   = 1,
        atol                = 1e-12,
        rtol                = 1e-12,
        methode             = "lsodes",
        use_Rcpp            = 0
      )

      res3 <-  TRAISIE::DAISIE_DE_trait_loglik_CS(
        datalist            = datalist,
        parameter           = parameter,
        num_observed_states = 2,
        num_hidden_states   = 1,
        atol                = 1e-12,
        rtol                = 1e-12,
        methode             = "lsodes",
        use_Rcpp            = 2
      )
      testthat::expect_equal(res1, res3, tol = 0.005)
  }
})
