test_that("all three n-sampling methods agree", {
  skip_if_not_installed("cubature")
  data("datalist_test_TRAISIE_n", package = "DAISIE")

  parameter <- list(
    c(2.546591, 1.2, 1, 0.2),
    c(2.678781, 2, 1.9, 3),
    c(0.009326754, 0.003, 0.002, 0.2),
    c(1.008583, 1, 2, 1.5),
    matrix(c(
      0,     0.001,  0.005, 0,
      0.001, 0,      0,     0.005,
      0.015, 0,      0,     0.005,
      0,     0.0025, 0.005, 0
    ), nrow = 4, byrow = TRUE),
    matrix(1, nrow = 4, ncol = 4)
  )

  loglik <- sapply(
    c("grid", "power_series", "hcubature"),
    function(method) {
      TRAISIE_loglik_CS(
        parameter = parameter,
        datalist = datalist_test_TRAISIE_n,
        methode = "odeint::runge_kutta_fehlberg78",
        atol = 1e-15,
        rtol = 1e-15,
        num_observed_states = 2,
        num_hidden_states = 2,
        cond = 1,
        num_threads = 1,
        verbose = FALSE,
        sampling = "n",
        n_integration = method
      )
    }
  )

  expect_true(all(is.finite(loglik)))
  expect_lt(max(loglik) - min(loglik), 1e-5)
})
