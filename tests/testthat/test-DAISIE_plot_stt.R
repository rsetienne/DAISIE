context("DAISIE_plot_stt")

test_that("use", {
  time <- 4
  n_times <- 10
  times <- seq(time, 0, length.out = n_times)
  df <- data.frame(
    Time = times, 
    nI = rep(0, n_times), 
    Endemic = rep(0, n_times), 
    Total = rep(0, n_times)
  )  
  expect_silent(
    DAISIE:::DAISIE_plot_stt(
      plot_plus_one = FALSE,
      time = 10,
      stt_q0.025 = df,
      stt_q0.25 = df,
      stt_average = df,
      stt_q0.75 = df,
      stt_q0.975 = df
    )
  )
})
