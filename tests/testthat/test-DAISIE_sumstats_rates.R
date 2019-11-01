context("DAISIE_sumstats_rates")

test_that("use", {
  out <- DAISIE_calc_sumstats_pcrates(pars = c(2, 2, 40, 0.1, 1),
                                      Apars = create_area_params(max_area = 1000,
                                                                 proportional_peak_t = 0.1,
                                                                 peak_sharpness = 1,
                                                                 total_island_age = 12),
                                      Epars = c(0.4, 0.6),
                                      totaltime = 10,
                                      island_ontogeny = 2,
                                      extcutoff = 1100,
                                      mainland_n = 1000,
                                      resol = 100)
  expect_true(is.list(out))
})
