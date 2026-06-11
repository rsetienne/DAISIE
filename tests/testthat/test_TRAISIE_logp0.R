test_that("logp0", {
  data("Galapagos_datalist", package = "DAISIE")
  datalist <- Galapagos_datalist

  datalist[[1]]$Mainland_pool_sizes <- c(550, 250)
  datalist[[1]]$M <- 1000

  parameter <- list(2.546591, 2.678781, 0.009326754, 1.008583,
                    matrix(c(0), nrow = 1), 0, 1)

  res1 <-  TRAISIE_logp0(
    datalist            = datalist,
    parameter           = parameter,
    num_observed_states = 1,
    num_hidden_states   = 1,
    atol                = 1e-10,
    rtol                = 1e-10,
    methode             = "ode45"
  )

  res2 <- DAISIE:::DAISIE_DE_logp0(
    island_age = datalist[[1]]$island_age,
    pars1 = c(2.546591, 2.678781, 2.678781, 0.009326754, 1.008583),
    methode = "ode45",
    reltolint = 1e-10,
    abstolint = 1e-10)

  testthat::expect_equal(res1, res2)

  res3 <-  TRAISIE_logp0(
    datalist            = datalist,
    parameter           = parameter,
    num_observed_states = 1,
    num_hidden_states   = 1,
    atol                = 1e-10,
    rtol                = 1e-10,
    methode             = "odeint::runge_kutta_cash_karp54"
  )
  testthat::expect_equal(res1, res3)
})
