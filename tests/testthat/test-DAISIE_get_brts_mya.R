context("DAISIE_get_brts_mya")

test_that("use", {
  data(Galapagos_datatable)
  branching_times_mya <- DAISIE_get_brts_mya(data_table = Galapagos_datatable)
  expect_true(class(branching_times_mya) == "numeric")
  expect_true(length(branching_times_mya) > 1)
  expect_true(all(branching_times_mya > 0))
})
