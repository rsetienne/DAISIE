context("DAISIE_extract_stt_median")

test_that("use", {
  data(islands_1type_1000reps)

  created <- as.data.frame(DAISIE_extract_stt_median(islands_1type_1000reps))
  # Trick to extract values:
  # paste0(DAISIE_extract_stt_median(islands_1type_1000reps)[ , "Total" ], ",", collapse = " ") # nolint
  expected <- data.frame(
    Time = c(4.00, 3.84, 3.68, 3.52, 3.36, 3.20, 3.04, 2.88, 2.72, 2.56, 2.40, 2.24, 2.08, 1.92, 1.76, 1.60, 1.44, 1.28, 1.12, 0.96, 0.80, 0.64, 0.48, 0.32, 0.16, 0.00), # nolint keep this one line
    nI = c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), # nolint keep this one line
    nA = c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1, 1, 2, 2, 2, 2), # nolint keep this one line
    nC = c(0, 0, 0, 2, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 12, 14, 14, 16, 17, 18, 19, 20, 20, 21, 22, 23), # nolint keep this one line
    Endemic = c(0, 0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 21, 22, 23, 24, 25), # nolint keep this one line
    Total = c(0, 1, 2, 4, 5, 6, 7, 9, 10, 10, 11.5, 13, 14, 15, 15, 17, 18, 19, 20, 22, 22, 23, 24, 24, 26, 26) # nolint keep this one line
  )
    
  expect_equal(created, expected)
})
