# For strongly peaked island ontogeny (beta function shape parameter ~114) and
# area exponents d = x = 0.1 the extinction rate mu * A^(-x) reaches the cutoff
# of 1e6 while the area is tiny, and the likelihood of the clades is far below
# the resolution of the ODE solution (probabilities < 1e-16). The C++ ODE
# solvers (dopri5, cash_karp54, bulirsch_stoer and rosenbrock4 alike) then stop
# with "cpp_daisie_cs_time: too many steps", and DAISIE_ML used to abort with
# this error instead of treating the parameter set as having likelihood 0.
# Whatever finite value other solvers (lsodes) return for such parameters is
# not stable under changes of the tolerances, i.e. it is numerical noise.

age <- 5
ap_vec <- function(ap) {
  unlist(ap)[c("max_area", "current_area", "proportional_peak_t",
               "total_island_age", "sea_level_amplitude",
               "sea_level_frequency", "island_gradient_angle")]
}
ap_low <- create_area_pars(1000, 100, 0.6, 10, 0, 0, 0)

test_that("a parameter set for which the ODE solver needs too many steps gets loglik -Inf", {
  skip_on_cran()
  set.seed(1)
  dl <- DAISIE_sim_cr(time = age, M = 40, pars = c(0.5, 0.3, Inf, 0.1, 0.4),
                      replicates = 1, plot_sims = FALSE, verbose = FALSE)[[1]]
  pars <- c(0.5, 0.3, Inf, 0.1, 0.4, 0.1, 0.1)   # lac, mu, K, gam, laa, d, x
  trpars <- pars / (1 + pars)
  trpars[3] <- 1
  expect_message(
    ll <- DAISIE_loglik_all_choosepar_time(
      trparsopt = trpars, trparsfix = numeric(0), idparsopt = 1:7,
      idparsfix = numeric(0), pars2 = c(40, 11, 0, 0, 1, 0),
      area_pars = ap_vec(ap_low), datalist = dl,
      methode = "odeint::runge_kutta_dopri5"
    ),
    "too many steps"
  )
  expect_equal(ll, -Inf)
})

test_that("DAISIE_ML does not abort when the start values need too many steps", {
  skip_on_cran()
  set.seed(1)
  dl <- DAISIE_sim_cr(time = age, M = 40, pars = c(0.5, 0.3, Inf, 0.1, 0.4),
                      replicates = 1, plot_sims = FALSE, verbose = FALSE)[[1]]
  expect_warning(
    fit <- suppressMessages(DAISIE_ML(
      datalist = dl, initparsopt = c(0.5, 0.3, 0.1, 0.4, 0.1, 0.1),
      idparsopt = c(1, 2, 4, 5, 6, 7), parsfix = Inf, idparsfix = 3,
      island_ontogeny = 1, sea_level = 0, area_pars = ap_low, ddmodel = 11,
      res = 40, methode = "odeint::runge_kutta_dopri5",
      optimmethod = "simplex", verbose = 0
    )),
    "likelihood that is equal to 0"
  )
  expect_true(is.na(fit$loglik))
})
