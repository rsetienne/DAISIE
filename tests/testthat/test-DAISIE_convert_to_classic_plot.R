test_that("abuse", {
  expect_error(
    DAISIE_convert_to_classic_plot("nonsense"),
    "'simulation_outputs' should be a set of simulation outputs"
  )
})


# test_that("use", {
#   utils::data("islands_1type_1000reps", package = "DAISIE")
#   simulation_outputs <- DAISIE::DAISIE_convert_to_classic_plot(
#     islands_1type_1000reps
#   )
#   expect_true(is_simulation_outputs(simulation_outuputs))
# })
