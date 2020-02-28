context("DAISIE_utils")

test_that("DAISIE_nonoceanic_spec output is silent", {
  expect_silent(DAISIE_nonoceanic_spec(prob_samp = 0.5,
                                       prob_nonend = 0.5,
                                       mainland_n = 100))
})

test_that("DAISIE_nonoceanic_spec output is a list of three vectors", {
  native_spec <- DAISIE_nonoceanic_spec(prob_samp = 0.1,
                                        prob_nonend = 0.9,
                                        mainland_n = 1000)
  expect_true(class(native_spec) == "list")
})

test_that("DAISIE_nonoceanic_spec samples native species
          when probability of sampling is non-zero", {
  native_spec <- DAISIE_nonoceanic_spec(prob_samp = 0.1,
                                        prob_nonend = 0.9,
                                        mainland_n = 1000)
  expect_true(is.list(native_spec))
  expect_true(is.vector(native_spec[[1]]))
  expect_true(is.numeric(native_spec[[1]]))
  expect_true(is.vector(native_spec[[2]]))
  expect_true(is.numeric(native_spec[[2]]))
  expect_true(is.vector(native_spec[[3]]))
  expect_true(is.numeric(native_spec[[3]]))
  expect_gt(length(native_spec[[1]]), 0)
  expect_gt(length(native_spec[[2]]), 0)
})

test_that("DAISIE_nonoceanic_spec samples no native species
          with zero probability of sampling", {
  prob_samp <- 0.0
  prob_nonend <- 0.9
  mainland_n <- 1000
  nonoceanic_sample <- DAISIE_nonoceanic_spec(prob_samp = prob_samp,
                                        prob_nonend = prob_nonend,
                                        mainland_n = mainland_n)
  expect_true(nonoceanic_sample$init_nonend_spec == 0)
  expect_true(nonoceanic_sample$init_end_spec == 0)
  expect_equal(length(nonoceanic_sample$mainland_spec), mainland_n)
})

test_that("DAISIE_nonoceanic_spec correctly samples number of species
          with seed", {
  set.seed(17)
  nonoceanic_sample <- DAISIE_nonoceanic_spec(prob_samp = 0.5,
                                              prob_nonend = 0.5,
                                              mainland_n = 100)
  expect_equivalent(nonoceanic_sample$init_nonend_spec_vec,
                    c(51, 33, 34,  8, 45,  4,
                      74, 85, 31, 75, 49, 21,
                      55, 92, 39, 81, 61, 41,
                      58, 24, 13, 26, 72, 42))
  expect_equivalent(nonoceanic_sample$init_end_spec_vec,
                    c(80,  5, 44,  3, 12, 25,
                      40, 17, 84,  1, 22, 79,
                      99, 16,  9, 78, 83, 14,
                      50, 18, 64, 20, 70, 69,
                      53, 28, 67, 93, 73,  7,
                      95, 30))
  expect_equivalent(nonoceanic_sample$mainland_spec,
                    c(2, 4, 6, 8, 10, 11, 13,
                      15, 19, 21, 23, 24, 26,
                      27, 29, 31, 32, 33, 34,
                      35, 36, 37, 38, 39, 41,
                      42, 43, 45, 46, 47, 48,
                      49, 51, 52, 54, 55, 56,
                      57, 58, 59, 60, 61, 62,
                      63, 65, 66, 68, 71, 72,
                      74, 75, 76, 77, 81, 82,
                      85, 86, 87, 88, 89, 90,
                      91, 92, 94, 96, 97, 98,
                      100))
})

test_that("DAISIE_spec_tables output is silent", {
  island_spec <- c()
  stt_table <- matrix(ncol = 4)
  colnames(stt_table) <- c("Time", "nI", "nA", "nC")
  totaltime <- 10
  timeval <- 0
  nonoceanic_sample <- list(init_nonend_spec = 4,
                            init_end_spec = 1,
                            init_nonend_spec_vec = c(28, 43, 15, 25),
                            init_end_spec_vec = 31,
                            mainland_spec = c(1:50))
  expect_silent(DAISIE_spec_tables(stt_table,
                                   totaltime,
                                   timeval,
                                   nonoceanic_sample,
                                   island_spec))
})

test_that("DAISIE_spec_tables produces correct output", {
  island_spec <- c()
  stt_table <- matrix(ncol = 4)
  colnames(stt_table) <- c("Time", "nI", "nA", "nC")
  totaltime <- 10
  timeval <- 0
  mainland_spec <- c(1:50)
  mainland_spec <- mainland_spec[-31]
  nonoceanic_sample <- list(init_nonend_spec = 4,
                            init_end_spec = 1,
                            init_nonend_spec_vec = c(28, 43, 15, 25),
                            init_end_spec_vec = 31,
                            mainland_spec = mainland_spec)

  nonoceanic_tables <- DAISIE_spec_tables(stt_table,
                                          totaltime,
                                          timeval,
                                          nonoceanic_sample,
                                          island_spec)
  expected_stt <- stt_table <- matrix(ncol = 4)
  colnames(expected_stt) <- c("Time", "nI", "nA", "nC")
  expected_stt[1, ] <- c(10, 4, 1, 0)
  expected_mainland_spec <- c(1, 2, 3, 4, 5, 6, 7, 8, 9,
                              10, 11, 12, 13, 14, 15, 16,
                              17, 18, 19, 20, 21, 22, 23,
                              24, 25, 26, 27, 28, 29, 30,
                              32, 33, 34, 35, 36, 37, 38,
                              39, 40, 41, 42, 43, 44, 45,
                              46, 47, 48, 49, 50)
  expected_island_spec <- matrix(ncol = 7, nrow = 5)
  expected_island_spec[1, ] <- c("28", "28", "0", "I", NA, NA, NA)
  expected_island_spec[2, ] <- c("43", "43", "0", "I", NA, NA, NA)
  expected_island_spec[3, ] <- c("15", "15", "0", "I", NA, NA, NA)
  expected_island_spec[4, ] <- c("25", "25", "0", "I", NA, NA, NA)
  expected_island_spec[5, ] <- c("31", "31", "0", "A", NA, NA, NA)
  expect_true(length(nonoceanic_tables) == 5)
  expect_true("stt_table" %in% names(nonoceanic_tables))
  expect_true("init_nonend_spec" %in% names(nonoceanic_tables))
  expect_true("init_end_spec" %in% names(nonoceanic_tables))
  expect_true("mainland_spec" %in% names(nonoceanic_tables))
  expect_true("island_spec" %in% names(nonoceanic_tables))
  expect_equal(nonoceanic_tables$stt_table, expected_stt)
  expect_equal(nonoceanic_tables$init_nonend_spec, 4)
  expect_equal(nonoceanic_tables$init_end_spec, 1)
  expect_equal(nonoceanic_tables$mainland_spec, expected_mainland_spec)
  expect_equal(nonoceanic_tables$island_spec, expected_island_spec)
})

test_that("creates singleton phylogeny", {
  tree <- create_singleton_phylo(age = 10)
  expect_true(class(tree) == "phylo")
})

test_that("translate_island_ontogeny", {
  expect_silent(translate_island_ontogeny("const"))
  expect_equal(translate_island_ontogeny("const"), 0)
  expect_equal(translate_island_ontogeny("beta"), 1)
  expect_false(is_island_ontogeny_input("ontogeny"))
})

test_that("translate_sea_level", {
  expect_silent(translate_sea_level("const"))
  expect_equal(translate_sea_level("const"), 0)
  expect_equal(translate_sea_level("sine"), 1)
  expect_false(is_sea_level_input("sea_level"))
})

test_that("counstspecies", {
  utils::data(Galapagos_datalist, package = "DAISIE")
  expect_equal(countspecies(Galapagos_datalist[[2]]), 1)
  expect_error(countspecies("nonsense"))
})

test_that("counttype1", {
  utils::data(Galapagos_datalist, package = "DAISIE")
  expect_equal(counttype1(Galapagos_datalist[[2]]), TRUE)
  expect_error(counttype1("nonsense"))
})

