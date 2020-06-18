context("DAISIE_dataprep")

test_that("DAISIE_dataprep produces a named list of length 9 for one type", {
  utils::data(Galapagos_datatable)
  output <- DAISIE_dataprep(datatable = Galapagos_datatable,
                  island_age = 4,
                  M = 1000,
                  verbose = FALSE)
  expect_length(output, 9)
  expect_true(is.list(output))
  expect_named(output[[1]], expected = c("island_age", "not_present"))
})

test_that("DAISIE_dataprep produces a named list of length 9 for two types", {
  utils::data(Galapagos_datatable)
  output <- DAISIE_dataprep(datatable = Galapagos_datatable,
                            island_age = 4,
                            M = 1000,
                            number_clade_types = 2,
                            list_type2_clades = "Finches",
                            verbose = FALSE)
  expect_length(output, 9)
  expect_true(is.list(output))
  expect_named(output[[1]], expected = c("island_age",
                                         "not_present_type1",
                                         "not_present_type2"))
})

test_that("DAISIE_dataprep produces a named list of length 9 for two types
          with the proportion of mainland pool which is type two set", {
  utils::data(Galapagos_datatable)
  output <- DAISIE_dataprep(datatable = Galapagos_datatable,
                            island_age = 4,
                            M = 1000,
                            number_clade_types = 2,
                            list_type2_clades = "Finches",
                            prop_type2_pool = 0.163,
                            verbose = FALSE)
  expect_length(output, 9)
  expect_true(is.list(output))
  expect_named(output[[1]], expected = c("island_age",
                                         "not_present_type1",
                                         "not_present_type2"))
})
