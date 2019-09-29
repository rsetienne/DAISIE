context("are_rates")

test_that("basic use", {
  expect_true(
    are_rates(
      create_rates(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = 0.4,
        ext_rate_max = 0.5,
        immig_rate_max = 0.6,
        clado_rate_max = 0.7,
        Tpars = list(trans_rate = 0.5, 
                     immig_rate2 = 0.1, 
                     ext_rate2 = 0.2, 
                     ana_rate2 = 0.3, 
                     clado_rate2 = 0.4, 
                     trans_rate2 = 0.5, 
                     M2 = 1000)
      )
    )
  )
  expect_true(
    are_rates(
      create_rates(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = 0.4,
        ext_rate_max = 0.5,
        immig_rate_max = 0.6,
        clado_rate_max = 0.7,
        Tpars = NULL
      )
    )
  )
})
