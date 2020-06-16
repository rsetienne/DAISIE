context("are_max_rates")

test_that("basic use", {
  expect_true(
    are_max_rates(
      list(
        immig_max_rate = 0.1,
        ext_max_rate = 0.2,
        ana_max_rate = 0.3,
        clado_max_rate = 0.4
      )
    )
  )
})

test_that("check returns FALSE when wrong", {

  expect_false(
    are_max_rates(
      list(
        immig_max_rate = 0.1,
        ext_max_rate = 0.2,
        ana_max_rate = 0.3
      )
    )
  )

  expect_false(
    are_max_rates(
      list(
        immig_max_rate = "nonsense",
        ext_max_rate = 0.2,
        ana_max_rate = 0.3,
        clado_max_rate = 0.4
      )
    )
  )
  expect_false(
    are_max_rates(
      list(
        clado_max_rate = 0.1,
        ext_max_rate = 0.2,
        ana_max_rate = 0.3
      )
    )
  )

  expect_false(
    are_max_rates(
      list(
        immig_max_rate = 0.1,
        ext_max_rate = 0.2,
        ana_max_rate = 0.3,
        clado_max_rate = "nonsense"
      )
    )
  )
  expect_false(
    are_max_rates(
      list(
        clado_max_rate = 0.1,
        immig_max_rate = 0.2,
        ana_max_rate = 0.3
      )
    )
  )

  expect_false(
    are_max_rates(
      list(
        immig_max_rate = 0.1,
        ext_max_rate = "nonsense",
        ana_max_rate = 0.3,
        clado_max_rate = 0.4
      )
    )
  )
  expect_false(
    are_max_rates(
      list(
        immig_max_rate = 0.1,
        ext_max_rate = 0.2,
        clado_max_rate = 0.3
      )
    )
  )

  expect_false(
    are_max_rates(
      list(
        immig_max_rate = 0.1,
        ext_max_rate = 0.2,
        ana_max_rate = "nonsense",
        clado_max_rate = 0.4
      )
    )
  )
  expect_false(
    are_max_rates(
      list(
        immig_max_rate = 0.1,
        ext_max_rate = 0.2,
        ana_max_rate = -1,
        clado_max_rate = 0.4
      )
    )
  )
  expect_false(
    are_max_rates(
      list(
        immig_max_rate = -1,
        ext_max_rate = 0.2,
        ana_max_rate = 0.2,
        clado_max_rate = 0.4
      )
    )
  )
  expect_false(
    are_max_rates(
      list(
        immig_max_rate = 0.1,
        ext_max_rate = -1,
        ana_max_rate = 0.3,
        clado_max_rate = 0.4
      )
    )
  )
  expect_false(
    are_max_rates(
      list(
        immig_max_rate = 0.1,
        ext_max_rate = 0.2,
        ana_max_rate = 0.5,
        clado_max_rate = -1
      )
    )
  )
})
