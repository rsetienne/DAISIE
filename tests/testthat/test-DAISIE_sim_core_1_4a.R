context("DAISIE_sim_core_1_4a")

test_that("abuse gam = 0", {
  expect_error(
    DAISIE:::DAISIE_sim_core_1_4a(
      time = 5,
      mainland_n = 1000,
      pars = c(1, 1, 20, 0, 1)
    ),
    regexp = "Rate of colonisation is zero. Island cannot be colonised.")
})

test_that("anagenesis sampling on DAISIE_sim_core_1_4a
          multiple immigrants works", {
  set.seed(1)
  expect_silent(
    out <- DAISIE:::DAISIE_sim_core_1_4a(
      time = 5,
      mainland_n = 1000,
      pars = c(0, 1, 20, 5, 5)
    )
  )

  expect_length(out, 2)
  expect_length(out[[1]], 1236)
  expect_length(out[[2]], 20)

})
