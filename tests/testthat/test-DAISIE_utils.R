context("daisie_utils")
test_that("creates singleton phylogeny", {
  tree <- create_singleton_phylo(age = 10)
  expect_true(class(tree) == "phylo")
})
