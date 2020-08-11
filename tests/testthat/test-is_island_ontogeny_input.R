context("is_island_ontogeny_input")

test_that("Returns TRUE when correct", {
  expect_true(
    is_island_ontogeny_input(island_ontogeny = "const")
  )
  expect_true(
    is_island_ontogeny_input(island_ontogeny = "beta")
  )
})

test_that("Returns FALSE when not string", {
  expect_false(
    is_island_ontogeny_input(island_ontogeny = 2)
  )
})

test_that("Returns FALSE when wrong string", {
  expect_false(
    is_island_ontogeny_input(island_ontogeny = "nonsense")
  )
})
