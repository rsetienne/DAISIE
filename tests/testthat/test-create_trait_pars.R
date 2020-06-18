context("create_trait_pars")

test_that("minimal use", {
  expect_silent(
    create_trait_pars(
      trans_rate = 0.5,
      immig_rate2 = 0.1,
      ext_rate2 = 0.2,
      ana_rate2 = 0.3,
      clado_rate2 = 0.4,
      trans_rate2 = 0.5,
      M2 = 1000
    )
  )
  out <- create_trait_pars(
    trans_rate = 0.5,
    immig_rate2 = 0.1,
    ext_rate2 = 0.2,
    ana_rate2 = 0.3,
    clado_rate2 = 0.4,
    trans_rate2 = 0.5,
    M2 = 1000
  )
  reference <- list(
    trans_rate = 0.5,
    immig_rate2 = 0.1,
    ext_rate2 = 0.2,
    ana_rate2 = 0.3,
    clado_rate2 = 0.4,
    trans_rate2 = 0.5,
    M2 = 1000
  )
  expect_equal(out, reference)
})

