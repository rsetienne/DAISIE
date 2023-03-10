test_that("Example 1", {
  skip("Plots: Run and manually inspect output")
  data(islands_1type_1000reps)
  expect_silent(
    DAISIE_plot_sims(
      island_replicates = islands_1type_1000reps
    )
  )
})

test_that("Example 2", {
  skip("Plots: Run and manually inspect output")
  data(islands_2types_1000reps)
  expect_silent(
    DAISIE_plot_sims(
      island_replicates = islands_2types_1000reps
    )
  )
})

test_that("Plot plus one", {
  skip("Plots: Run and manually inspect output")
  data(islands_1type_1000reps)
  expect_silent(
    DAISIE_plot_sims(
      island_replicates = islands_1type_1000reps,
      plot_plus_one = TRUE
    )
  )
  expect_silent(
    DAISIE_plot_sims(
      island_replicates = islands_1type_1000reps,
      plot_plus_one = FALSE
    )
  )
})
