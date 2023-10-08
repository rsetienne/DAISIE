test_that("basic use", {
  testthat::expect_true(
    are_rates(
      list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = 0.4
      )
    )
  )
})

test_that("check returns FALSE when wrong", {

  testthat::expect_false(
    are_rates(
      list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3
      )
    )
  )

  testthat::expect_false(
    are_rates(
      list(
        immig_rate = "nonsense",
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = 0.4
      )
    )
  )
  testthat::expect_false(
    are_rates(
      list(
        clado_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3
      )
    )
  )

  testthat::expect_false(
    are_rates(
      list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = "nonsense"
      )
    )
  )
  testthat::expect_false(
    are_rates(
      list(
        clado_rate = 0.1,
        immig_rate = 0.2,
        ana_rate = 0.3
      )
    )
  )

  testthat::expect_false(
    are_rates(
      list(
        immig_rate = 0.1,
        ext_rate = "nonsense",
        ana_rate = 0.3,
        clado_rate = 0.4
      )
    )
  )
  testthat::expect_false(
    are_rates(
      list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        clado_rate = 0.3
      )
    )
  )

  testthat::expect_false(
    are_rates(
      list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = "nonsense",
        clado_rate = 0.4
      )
    )
  )
  testthat::expect_false(
    are_rates(
      list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = -1,
        clado_rate = 0.4
      )
    )
  )
  testthat::expect_false(
    are_rates(
      list(
        immig_rate = -1,
        ext_rate = 0.2,
        ana_rate = 0.2,
        clado_rate = 0.4
      )
    )
  )
  testthat::expect_false(
    are_rates(
      list(
        immig_rate = 0.1,
        ext_rate = -1,
        ana_rate = 0.3,
        clado_rate = 0.4
      )
    )
  )
  testthat::expect_false(
    are_rates(
      list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.5,
        clado_rate = -1
      )
    )
  )
})

test_that("including two trait states", {
  testthat::expect_true(
    are_rates(
      list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = 0.4,
        immig_rate2 = 0.1,
        ext_rate2 = 0.2,
        ana_rate2 = 0.3,
        clado_rate2 = 0.4,
        trans_rate = 0.5,
        trans_rate2 = 0.5
      )
    )
  )
})

test_that("check returns FALSE when wrong", {

  testthat::expect_false(
    are_rates(
      list(
        immig_rate = 0.1,
        clado_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3,
        ext_rate2 = 0.2,
        ana_rate2 = 0.3,
        clado_rate2 = 0.4,
        trans_rate = 0.5,
        trans_rate2 = 0.5
      )
    )
  )

  testthat::expect_false(
    are_rates(
      list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = 0.4,
        immig_rate2 = "nonsense",
        ext_rate2 = 0.2,
        ana_rate2 = 0.3,
        clado_rate2 = 0.4,
        trans_rate = 0.5,
        trans_rate2 = 0.5
      )
    )
  )

  testthat::expect_false(
    are_rates(
      list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = 0.4,
        immig_rate2 = 0.1,
        ana_rate2 = 0.3,
        clado_rate2 = 0.4,
        trans_rate = 0.5,
        trans_rate2 = 0.5
      )
    )
  )
  testthat::expect_false(
    are_rates(
      list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = 0.4,
        immig_rate2 = 0.1,
        ext_rate2 = "nonsense",
        ana_rate2 = 0.3,
        clado_rate2 = 0.4,
        trans_rate = 0.5,
        trans_rate2 = 0.5
      )
    )
  )

  testthat::expect_false(
    are_rates(
      list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = 0.4,
        immig_rate2 = 0.1,
        ext_rate2 = 0.2,
        clado_rate2 = 0.4,
        trans_rate = 0.5,
        trans_rate2 = 0.5
      )
    )
  )
  testthat::expect_false(
    are_rates(
      list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = 0.4,
        immig_rate2 = 0.1,
        ext_rate2 = 0.2,
        ana_rate2 = "nonsense",
        clado_rate2 = 0.4,
        trans_rate = 0.5,
        trans_rate2 = 0.5
      )
    )
  )

  testthat::expect_false(
    are_rates(
      list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = 0.4,
        immig_rate2 = 0.1,
        ext_rate2 = 0.2,
        ana_rate2 = 0.3,
        trans_rate = 0.5,
        trans_rate2 = 0.5
      )
    )
  )
  testthat::expect_false(
    are_rates(
      list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = -1,
        clado_rate = 0.4
      )
    )
  )
  testthat::expect_false(
    are_rates(
      list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = 0.4,
        immig_rate2 = 0.1,
        ext_rate2 = 0.2,
        ana_rate2 = 0.3,
        clado_rate2 = "nonsense",
        trans_rate = 0.5,
        trans_rate2 = 0.5
      )
    )
  )
  testthat::expect_false(
    are_rates(
      list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = 0.4,
        immig_rate2 = 0.1,
        ext_rate2 = 0.2,
        ana_rate2 = 0.3,
        clado_rate2 = 0.4,
        trans_rate2 = 0.5
      )
    )
  )
  testthat::expect_false(
    are_rates(
      list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = 0.4,
        immig_rate2 = 0.1,
        ext_rate2 = 0.2,
        ana_rate2 = 0.3,
        clado_rate2 = 0.4,
        trans_rate = "nonsense",
        trans_rate2 = 0.5
      )
    )
  )
  testthat::expect_false(
    are_rates(
      list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = 0.4,
        immig_rate2 = 0.1,
        ext_rate2 = 0.2,
        ana_rate2 = 0.3,
        clado_rate2 = 0.4,
        trans_rate = 0.5
      )
    )
  )
  testthat::expect_false(
    are_rates(
      list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = 0.4,
        immig_rate2 = 0.1,
        ext_rate2 = 0.2,
        ana_rate2 = 0.3,
        clado_rate2 = 0.4,
        trans_rate = 0.5,
        trans_rate2 = "nonsense"
      )
    )
  )
  testthat::expect_false(
    are_rates(
      list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = 0.4,
        immig_rate2 = -1,
        ext_rate2 = 0.2,
        ana_rate2 = 0.3,
        clado_rate2 = 0.4,
        trans_rate = 0.5,
        trans_rate2 = 0.5
      )
    )
  )
  testthat::expect_false(
    are_rates(
      list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = 0.4,
        immig_rate2 = 0.1,
        ext_rate2 = -1,
        ana_rate2 = 0.3,
        clado_rate2 = 0.4,
        trans_rate = 0.5,
        trans_rate2 = 0.5
      )
    )
  )
  testthat::expect_false(
    are_rates(
      list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = 0.4,
        immig_rate2 = 0.1,
        ext_rate2 = 0.2,
        ana_rate2 = -1,
        clado_rate2 = 0.4,
        trans_rate = 0.5,
        trans_rate2 = 0.5
      )
    )
  )
  testthat::expect_false(
    are_rates(
      list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = 0.4,
        immig_rate2 = 0.1,
        ext_rate2 = 0.2,
        ana_rate2 = 0.3,
        clado_rate2 = -1,
        trans_rate = 0.5,
        trans_rate2 = 0.5
      )
    )
  )
  testthat::expect_false(
    are_rates(
      list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = 0.4,
        immig_rate2 = 0.1,
        ext_rate2 = 0.2,
        ana_rate2 = 0.3,
        clado_rate2 = 0.4,
        trans_rate = -1,
        trans_rate2 = 0.5
      )
    )
  )
  testthat::expect_false(
    are_rates(
      list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.5,
        clado_rate = 1,
        clado_rate = 0.4,
        immig_rate2 = 0.1,
        ext_rate2 = 0.2,
        ana_rate2 = 0.3,
        clado_rate2 = 0.4,
        trans_rate = 0.5,
        trans_rate2 = -1
      )
    )
  )
  testthat::expect_false(
    are_rates(
      list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.5,
        clado_rate = -1,
        clado_rate = 0.4,
        immig_rate2 = 0.1,
        ext_rate2 = 0.2,
        ana_rate2 = 0.3,
        clado_rate2 = 0.4,
        trans_rate = 0.5,
        trans_rate2 = 1
      )
    )
  )
})


