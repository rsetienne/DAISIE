context("are_trait_pars")

test_that("minimal use", {
  expect_true(
    are_trait_pars(
      create_trait_pars(
        trans_rate = 0.5,
        immig_rate2 = 0.1,
        ext_rate2 = 0.2,
        ana_rate2 = 0.3,
        clado_rate2 = 0.4,
        trans_rate2 = 0.5,
        M2 = 1000)))
})
