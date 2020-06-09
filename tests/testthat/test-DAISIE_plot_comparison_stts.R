test_that("abuse", {
  expect_error(
    DAISIE_plot_comparison_stts(
      time = 5,
      plot_lists_simulations = "?",
      plot_lists_simulations_MLE = "?",
      type = "nonsense"
    ),
    "type should be 'all_species', 'type1_species' or 'type2_species'"
  )
})
