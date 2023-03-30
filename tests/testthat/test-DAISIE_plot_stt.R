test_that("use", {
  skip("Plots: Run and manually inspect output")
  utils::data(islands_1type_1000reps, package = "DAISIE")
  plot_lists <- DAISIE_convert_to_classic_plot(
    simulation_outputs = islands_1type_1000reps
  )
  type <- names(plot_lists)[1]
  expect_silent(
    DAISIE_plot_stt(
      plot_plus_one = FALSE,
      time = 4,
      plot_lists = plot_lists,
      type = type
    )
  )
})
