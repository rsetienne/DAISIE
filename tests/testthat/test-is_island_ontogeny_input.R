test_that("Returns TRUE when correct", {
  testthat::expect_true(
    is_island_ontogeny_input(island_ontogeny = "const")
  )
  testthat::expect_true(
    is_island_ontogeny_input(island_ontogeny = "beta")
  )
})

test_that("Returns FALSE when not string", {
  testthat::expect_false(
    is_island_ontogeny_input(island_ontogeny = 2)
  )
})

test_that("Returns FALSE when wrong string", {
  testthat::expect_false(
    is_island_ontogeny_input(island_ontogeny = "nonsense")
  )
})
