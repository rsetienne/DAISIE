test_that("use", {
  expect_false(is_simulation_outputs("nonsense"))
  expect_true(is_simulation_outputs(DAISIE_sim_constant_rate(
    time = 0.4,
    M = 10,
    pars = c(2, 2, Inf, 0.001, 1),
    replicates = 2,
    plot_sims = FALSE,
    verbose = FALSE)))

  expect_true(is_simulation_outputs(DAISIE::DAISIE_sim_time_dependent(
    time = 2,
    M = 500,
    pars = c(0.00001, 1, 0.05, 0.001, 1),
    replicates = 1,
    area_pars = create_area_pars(
      max_area = 10000,
      proportional_peak_t = 0.1,
      peak_sharpness = 1,
      total_island_age = 4,
      sea_level_amplitude = 0,
      sea_level_frequency = 0,
      island_gradient_angle = 0),
    ext_pars = c(0.1, 15),
    island_ontogeny = "beta",
    plot_sims = FALSE,
    verbose = FALSE)
  )
  )
})
