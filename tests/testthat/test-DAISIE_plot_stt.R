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
    DAISIE_plot_stt(
      plot_plus_one = FALSE,
      time = 10,
      stt_q0.025_all = df,
      stt_q0.25_all = df,
      stt_average_all = df,
      stt_q0.75_all = df,
      stt_q0.975_all = df
    )
  )
})
