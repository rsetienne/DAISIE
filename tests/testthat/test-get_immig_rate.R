context("get_immig_rate")

test_that("immig rate plots", {
  calc_immig <- c()
  timepoints <- seq(0, 10, by = 0.1)
  for (i in 1:length(timepoints)) {
    calc_immig[i] <- get_immig_rate(
      timepoints[i],
      totaltime = 10,
      gam = 0.001,
      ddmodel_sim = 11,
      hyper_pars = NULL,
      area_pars = create_area_pars(5000, 0.5, 1, 15, 0, 0),
      dist_pars = NULL,
      island_ontogeny = 1,
      sea_level = 0,
      num_spec = 5,
      K = 0.05,
      mainland_n = 1000
      )
  }
  expected_immig <- c(0.0000000, 0.8771152,0.9128143, 0.9285714, 0.9379296,
                      0.9442914, 0.9489690, 0.9525895, 0.9554957, 0.9578924,
                      0.9599108, 0.9616392, 0.9631395, 0.9644566, 0.9656238,
                      0.9666667, 0.9676049, 0.9684542, 0.9692271, 0.9699338,
                      0.9705826, 0.9711805, 0.9717333, 0.9722460, 0.9727228,
                      0.9731672, 0.9735824, 0.9739710, 0.9743355, 0.9746779,
                      0.9750000, 0.9753034, 0.9755896, 0.9758598, 0.9761151,
                      0.9763567, 0.9765854, 0.9768020, 0.9770073, 0.9772020,
                      0.9773866, 0.9775619, 0.9777282, 0.9778861, 0.9780360,
                      0.9781782, 0.9783132, 0.9784412, 0.9785627, 0.9786778,
                      0.9787868, 0.9788900, 0.9789876, 0.9790797, 0.9791667,
                      0.9792486, 0.9793256, 0.9793979, 0.9794655, 0.9795287,
                      0.9795876, 0.9796422, 0.9796926, 0.9797390, 0.9797814,
                      0.9798198, 0.9798544, 0.9798852, 0.9799123, 0.9799357,
                      0.9799554, 0.9799715, 0.9799840, 0.9799929, 0.9799982,
                      0.9800000, 0.9799982, 0.9799929, 0.9799840, 0.9799715,
                      0.9799554, 0.9799357, 0.9799123, 0.9798852, 0.9798544,
                      0.9798198, 0.9797814, 0.9797390, 0.9796926, 0.9796422,
                      0.9795876, 0.9795287, 0.9794655, 0.9793979, 0.9793256,
                      0.9792486, 0.9791667, 0.9790797, 0.9789876, 0.9788900,
                      0.9787868)
  expect_equal(calc_immig, expected_immig)
})

test_that("classic behavior", {
  carr_cap <- 10
  ps_imm_rate <- 0.1
  n_island_species <- 5
  n_mainland_species <- 2
  expected <- DAISIE_calc_clade_imm_rate(
    ps_imm_rate = ps_imm_rate,
    n_island_species = n_island_species,
    n_mainland_species = n_mainland_species,
    carr_cap = carr_cap
  )
  created <- get_immig_rate(
    timeval = 1.0,
    totaltime = 10.0,
    gam = ps_imm_rate,
    ddmodel_sim = 11,
    hyper_pars = NULL,
    area_pars =  NULL,
    dist_pars = NULL,
    island_ontogeny = 0,
    sea_level = 0,
    num_spec = n_island_species,
    K = carr_cap,
    mainland_n = n_mainland_species
  )
  expect_equal(expected, created)
})
