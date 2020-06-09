context("are_max_rates_gt_rates")

test_that("basic use", {
  expect_true(
    DAISIE:::are_max_rates_gt_rates(
      rates = list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = 0.4
      ),
      max_rates = list(
        immig_max_rate = 0.2,
        ext_max_rate = 0.3,
        ana_max_rate = 0.4,
        clado_max_rate = 0.5
      )
    )
  )
  expect_true(
    DAISIE:::are_max_rates_gt_rates(
      rates = list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = 0.4
      ),
      max_rates = list(
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
    DAISIE:::are_max_rates_gt_rates(
      rates = list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = 0.4
      ),
      max_rates = list(
        immig_max_rate = 0.02,
        ext_max_rate = 0.3,
        ana_max_rate = 0.4,
        clado_max_rate = 0.5
      )
    )
  )

  expect_false(
    DAISIE:::are_max_rates_gt_rates(
      rates = list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = 0.4
      ),
      max_rates = list(
        immig_max_rate = 0.2,
        ext_max_rate = 0.03,
        ana_max_rate = 0.4,
        clado_max_rate = 0.5
      )
    )
  )
  expect_false(
    DAISIE:::are_max_rates_gt_rates(
      rates = list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = 0.4
      ),
      max_rates = list(
        immig_max_rate = 0.2,
        ext_max_rate = 0.3,
        ana_max_rate = 0.04,
        clado_max_rate = 0.5
      )
    )
  )

  expect_false(
    DAISIE:::are_max_rates_gt_rates(
      rates = list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = 0.4
      ),
      max_rates = list(
        immig_max_rate = 0.2,
        ext_max_rate = 0.3,
        ana_max_rate = 0.4,
        clado_max_rate = 0.05
      )
    )
  )
  expect_false(
    DAISIE:::are_max_rates_gt_rates(
      rates = list(
        immig_rate = "nonsense",
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = 0.4
      ),
      max_rates = list(
        immig_max_rate = 0.2,
        ext_max_rate = 0.3,
        ana_max_rate = 0.4,
        clado_max_rate = 0.05
      )
    )
  )
  expect_false(
    DAISIE:::are_max_rates_gt_rates(
      rates = list(
        immig_rate = 0.1,
        ext_rate = 0.2,
        ana_rate = 0.3,
        clado_rate = 0.4
      ),
      max_rates = list(
        immig_max_rate = "nonsense",
        ext_max_rate = 0.3,
        ana_max_rate = 0.4,
        clado_max_rate = 0.05
      )
    )
  )
})



