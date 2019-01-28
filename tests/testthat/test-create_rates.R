context("create_rates")

test_that("basic use", {
  expect_silent(
    create_rates(
      immig_rate = 0.1,
      ext_rate = 0.2,
      ana_rate = 0.3,
      clado_rate = 0.4,
      ext_rate_max = 0.5,
      immig_rate_max = 0.6,
      clado_rate_max = 0.7
    )
  )
})
