context("DAISIE_sumstats_rates")
test_that("use simple ontogeny code", {
  out <- DAISIE_calc_sumstats_pcrates(
    pars = c(2, 2, 40, 0.1, 1),
    area_pars = create_area_pars(max_area = 1000,
                                 proportional_peak_t = 0.1,
                                 peak_sharpness = 1,
                                 total_island_age = 12,
                                 sea_level_amplitude = 0,
                                 sea_level_frequency = 0,
                                 island_gradient_angle = 0),
    ext_pars = c(0.4, 0.6),
    totaltime = 10,
    island_ontogeny = 1,
    extcutoff = 100,
    mainland_n = 1000,
    resol = 100
  )
  expect_true(is.list(out))
})
