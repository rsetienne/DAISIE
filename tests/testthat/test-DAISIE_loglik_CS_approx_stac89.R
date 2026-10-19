context("DAISIE_approx for stac 8 and 9")

# Regression test: under DAISIE_approx the stac 8/9 part was integrated and
# added to the loglik twice, giving a loglik that differed strongly from the
# exact one (about -18.0 versus -12.7 in this example).
test_that("DAISIE_approx agrees with DAISIE for stac 8 and 9", {
  datalist <- list(
    list(island_age = 5, not_present = 20),
    list(colonist_name = "a", branching_times = c(5, 3, 2, 1), stac = 8,
         missing_species = 0, type1or2 = 1),
    list(colonist_name = "b", branching_times = c(5, 3, 2, 1), stac = 9,
         missing_species = 0, type1or2 = 1)
  )
  pars1 <- c(0.5, 0.3, 20, 0.05, 0.4)
  loglik <- function(fun) {
    DAISIE_loglik_all(
      pars1 = pars1, pars2 = c(100, 11, 0, 0), datalist = datalist,
      methode = "lsodes",
      CS_version = list(model = 1, function_to_optimize = fun, sampling = "n")
    )
  }
  ll_exact <- loglik("DAISIE")
  ll_approx <- loglik("DAISIE_approx")
  expect_equal(ll_approx, ll_exact, tolerance = 1e-3)
  expect_equal(ll_exact, -12.69575, tolerance = 1e-6)
})
