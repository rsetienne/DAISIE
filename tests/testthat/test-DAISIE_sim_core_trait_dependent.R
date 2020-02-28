context("DAISIE_sim_core_trait_dependent")

test_that("nonontogeny oceanic trait_dependnet island should run silent IW", {
  set.seed(234567890)
  expect_silent(
    DAISIE:::DAISIE_sim_core_trait_dependent(
      time = 10,
      mainland_n = 100,
      pars = c(0.0001, 2.2, 0.005, 0.001, 1),
      trait_pars = create_trait_pars(
        trans_rate = 0,
        immig_rate2 = 0.001,
        ext_rate2 = 0.2,
        ana_rate2 = 1,
        clado_rate2 = 0.004,
        trans_rate2 = 0,
        M2 = 200)
    )
  )
})

test_that("Ontogeny oceanic should run silent CS", {
  set.seed(420)
  expect_silent(
    DAISIE:::DAISIE_sim_core_trait_dependent(
      time = 10,
      mainland_n = 0,
      pars = c(0.0001, 2.2, 0.005, 0.001, 1),
      trait_pars = create_trait_pars(
        trans_rate = 0,
        immig_rate2 = 0.001,
        ext_rate2 = 0.2,
        ana_rate2 = 1,
        clado_rate2 = 0.004,
        trans_rate2 = 0,
        M2 = 200)
    )
  )
})
