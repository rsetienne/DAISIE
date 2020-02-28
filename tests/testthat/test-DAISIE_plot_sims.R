context("DAISIE_plot_sims")

test_that("Example 1", {

  data(islands_1type_1000reps)
  expect_silent(
    DAISIE_plot_sims(
      island_replicates = islands_1type_1000reps
    )
  )
})

test_that("Example 2", {
  data(islands_2types_1000reps)
  expect_silent(
    DAISIE_plot_sims(
     island_replicates = islands_2types_1000reps
    )
  )
})

test_that("Plot plus one", {
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
