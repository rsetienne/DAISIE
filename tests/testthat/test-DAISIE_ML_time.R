context("DAISIE_ML_time")

area <- create_area_pars(
  max_area = 10, current_area = 5, proportional_peak_t = 0.3,
  total_island_age = 10, sea_level_amplitude = 0,
  sea_level_frequency = 0, island_gradient_angle = 0
)

test_that("ML on the time-dependent model improves on the start and is finite", {
  skip_on_cran()
  set.seed(1)
  pars <- c(0.5, 0.3, 0.5, 0.2, 0.3)
  sim <- DAISIE_sim_time_dep(
    time = 5, M = 30, pars = pars, replicates = 1, area_pars = area,
    hyper_pars = create_hyper_pars(d = 0.3, x = 1.2),
    island_ontogeny = "beta", sea_level = "const",
    plot_sims = FALSE, verbose = FALSE
  )
  dl <- sim[[1]]
  ll0 <- DAISIE_loglik_all(
    pars1 = c(0.2, 0.2, 0.5, 0.05, 0.3, 0.3, 1.2, unlist(area)[1:7]),
    pars2 = c(30, 11, 0, 0, 1, 0), datalist = dl,
    methode = "odeint::runge_kutta_dopri5"
  )
  fit <- DAISIE_ML(
    datalist = dl, initparsopt = c(0.2, 0.05), idparsopt = c(1, 4),
    parsfix = c(0.2, 0.5, 0.3, 0.3, 1.2), idparsfix = c(2, 3, 5, 6, 7),
    island_ontogeny = 1, sea_level = 0, area_pars = area, res = 30,
    ddmodel = 11, cond = 0, tol = c(1e-8, 1e-8, 1e-8), optimmethod = 'simplex',
    verbose = 0, maxiter = 300
  )
  expect_true(is.finite(fit$loglik))
  expect_true(fit$loglik >= ll0 - 1e-8)
  expect_true(all(c(fit$lambda_c, fit$gamma) >= 0))
})

test_that("DAISIE_ML_time validates its input", {
  dl <- list(
    list(island_age = 5, not_present = 10),
    list(colonist_name = "a", branching_times = c(5, 2), stac = 2,
         missing_species = 0, type1or2 = 1)
  )
  args <- list(
    datalist = dl, initparsopt = c(0.2, 0.05), idparsopt = c(1, 4),
    parsfix = c(0.2, 0.5, 0.3, 0.3, 1.2), idparsfix = c(2, 3, 5, 6, 7),
    area_pars = area
  )
  expect_error(
    do.call(DAISIE_ML_time, modifyList(args, list(idparsfix = c(2, 3, 5, 6, 6)))),
    "must be specified exactly once"
  )
  expect_error(
    do.call(DAISIE_ML_time, modifyList(args, list(parsfix = c(-1, 0.5, 0.3, 0.3, 1.2)))),
    "non-negative"
  )
  area_young <- area
  area_young$total_island_age <- 4
  expect_error(
    do.call(DAISIE_ML_time, modifyList(args, list(area_pars = area_young))),
    "island age"
  )
  expect_error(
    DAISIE_ML(datalist = dl, initparsopt = 0.2, idparsopt = 1,
              parsfix = c(0.2, 0.5, 0.3, 0.3, 0.1), idparsfix = 2:6,
              island_ontogeny = 1, ddmodel = 11),
    "area_pars"
  )
})
