test_that("abuse", {
  
  expect_error(
    DAISIE_convert_to_classic_plot("nonsense"),
    "'simulation_outputs' should be a set of simulation outputs"
  )
})
