test_that("area windows bound the island area", {
  ap <- create_area_pars(
    max_area = 1000, current_area = 500, proportional_peak_t = 0.5,
    total_island_age = 10, sea_level_amplitude = 0,
    sea_level_frequency = 0, island_gradient_angle = 0
  )
  tt <- 7
  pk <- calc_peak(tt, ap)
  amax <- get_global_max_area(tt, ap, pk, 1, 0)
  amin <- get_global_min_area(tt, ap, pk, 1, 0)
  w <- get_area_windows(tt, ap, pk, 1, 0, amin, amax)
  expect_equal(w$end[length(w$end)], tt)
  expect_true(all(diff(w$end) > 0))
  expect_true(all(w$Amin <= w$Amax))
  set.seed(5)
  times <- stats::runif(2000, 0, tt)
  for (t in times) {
    k <- which(w$end >= t)[1]
    a <- island_area(t, tt, ap, pk, 1, 0)
    expect_true(a >= w$Amin[k] - 1e-12 && a <= w$Amax[k] + 1e-12)
  }
  # much tighter than the global bounds for most of the simulation
  expect_gt(length(w$end), 1)
  expect_true(stats::median(w$Amin) > 100 * max(amin, 1e-3))
})

test_that("windowed thinning gives the same distribution as global bounds", {
  skip_on_cran()
  ap <- create_area_pars(
    max_area = 1000, current_area = 500, proportional_peak_t = 0.5,
    total_island_age = 10, sea_level_amplitude = 0,
    sea_level_frequency = 0, island_gradient_angle = 0
  )
  tt <- 1
  pk <- calc_peak(tt, ap)
  amax <- get_global_max_area(tt, ap, pk, 1, 0)
  amin <- get_global_min_area(tt, ap, pk, 1, 0)
  global <- list(end = tt, Amin = amin, Amax = amax)
  run <- function(windows, n) {
    vapply(seq_len(n), function(i) {
      out <- DAISIE_sim_core_time_dep(
        time = tt, mainland_n = 1, pars = c(2, 1, 0.01, 1, 0.5),
        nonoceanic_pars = c(0, 0), island_ontogeny = "beta",
        sea_level = "const",
        hyper_pars = create_hyper_pars(d = 0.2, x = 1),
        area_pars = ap, peak = pk, Amax = amax, Amin = amin,
        extcutoff = 20, area_windows = windows
      )
      st <- out$stt_table
      sum(as.numeric(st[nrow(st), 2:4]))
    }, numeric(1))
  }
  set.seed(42)
  n <- 400
  a <- run(global, n)
  b <- run(NULL, n)
  se <- sqrt(stats::var(a) / n + stats::var(b) / n)
  expect_lt(abs(mean(a) - mean(b)), 4.5 * se)
  expect_lt(abs(mean(a > 0) - mean(b > 0)), 4.5 * sqrt(0.25 / n * 2))
})
