test_that("range codes and strings", {
  occ <- c(1, 0, 1)
  expect_equal(DAISIE_A_range_code(occ), 5)
  expect_equal(DAISIE_A_range_string(occ), "101")
  expect_equal(DAISIE_A_code_to_occ(5, 3), occ)
  m <- rbind(c(1, 0, 0), c(0, 1, 1))
  expect_equal(DAISIE_A_range_code(m), c(1, 6))
  expect_equal(DAISIE_A_range_string(m), c("100", "011"))
})

test_that("number of allowed ranges under restrictions", {
  base <- function(...) {
    DAISIE_A_create_pars(n_islands = 5, lac = 1, mu = 1, gam = 1, laa = 1,
                         ...)
  }
  expect_length(DAISIE_A_ranges(base())$sets, 31)
  expect_length(DAISIE_A_ranges(base(max_range_size = 2))$sets, 15)
  expect_length(
    DAISIE_A_ranges(base(connectivity = "line", connected_ranges = TRUE))$sets,
    15
  )
  expect_length(
    DAISIE_A_ranges(base(connectivity = "ring", connected_ranges = TRUE))$sets,
    21
  )
})

test_that("connectivity of sets", {
  adj <- DAISIE_A_connectivity(5, "line") > 0
  expect_true(DAISIE_A_is_connected(c(2, 3, 4), adj))
  expect_false(DAISIE_A_is_connected(c(1, 3), adj))
  expect_equal(DAISIE_A_components(c(1, 2, 4, 5), adj), list(c(1, 2), c(4, 5)))
})

test_that("admissible vicariance splits", {
  mk <- function(n, conn, mode, ...) {
    DAISIE_A_create_pars(n_islands = n, lac = 1, mu = 1, gam = 1, laa = 1,
                         connectivity = conn, vic_mode = mode, ...)
  }
  # range 1-2-3 on a line
  expect_length(DAISIE_A_admissible_splits(1:3, mk(3, "line", "connected")), 2)
  expect_length(DAISIE_A_admissible_splits(1:3, mk(3, "line", "all")), 3)
  expect_length(DAISIE_A_admissible_splits(1:3, mk(3, "line", "single")), 3)
  expect_length(DAISIE_A_admissible_splits(1:2, mk(3, "line", "single")), 1)
  # on a full graph "connected" equals "all"
  expect_length(DAISIE_A_admissible_splits(1:4, mk(4, "full", "connected")),
                7)
  # a range of three separate blocks on a line splits along the gaps
  splits <- DAISIE_A_admissible_splits(c(1, 3, 5), mk(5, "line", "connected"))
  expect_length(splits, 3)
  # blocks {1,2} and {4,5}: only one split
  splits <- DAISIE_A_admissible_splits(c(1, 2, 4, 5),
                                       mk(5, "line", "connected"))
  expect_equal(splits, list(c(1, 2)))
  # every split is a proper subset containing the first island
  splits <- DAISIE_A_admissible_splits(1:5, mk(5, "full", "all"))
  expect_length(splits, 15)
  expect_true(all(vapply(splits, function(a) 1 %in% a && length(a) < 5,
                         logical(1))))
  # connected_ranges restricts all modes to connected parts
  expect_length(DAISIE_A_admissible_splits(
    1:4, mk(4, "line", "all", connected_ranges = TRUE)), 3)
})

test_that("sampled vicariance splits are uniform over admissible splits", {
  pars <- DAISIE_A_create_pars(n_islands = 4, lac = 1, mu = 1, gam = 1,
                               laa = 1, connectivity = "line",
                               vic_mode = "connected")
  cache <- new.env()
  set.seed(1)
  draws <- replicate(3000, paste(DAISIE_A_sample_split(1:4, pars, cache),
                                 collapse = ""))
  tab <- table(draws)
  expect_equal(sort(names(tab)), c("1", "12", "123"))
  expect_gt(stats::chisq.test(tab)$p.value, 0.001)
  # rejection sampling for large ranges in "all" mode
  pars13 <- DAISIE_A_create_pars(n_islands = 13, lac = 1, mu = 1, gam = 1,
                                 laa = 1, vic_mode = "all")
  a <- DAISIE_A_sample_split(1:13, pars13, cache)
  expect_true(length(a) >= 1 && length(a) <= 12)
})

test_that("number of admissible splits matches the enumeration", {
  for (conn in c("line", "ring", "full")) {
    for (mode in c("connected", "all", "single")) {
      for (cr in c(FALSE, TRUE)) {
        pars <- DAISIE_A_create_pars(n_islands = 5, lac = 1, mu = 1, gam = 1,
                                     laa = 1, connectivity = conn,
                                     vic_mode = mode, connected_ranges = cr)
        cache <- new.env()
        for (set in DAISIE_A_ranges(pars)$sets) {
          expect_equal(DAISIE_A_n_splits(set, pars, cache),
                       length(DAISIE_A_admissible_splits(set, pars)))
        }
      }
    }
  }
  pars <- DAISIE_A_create_pars(n_islands = 6, lac = 1, mu = 1, gam = 1,
                               laa = 1, connectivity = "line")
  # a contiguous range on a line: one split per link between neighbours
  expect_equal(DAISIE_A_n_splits(2:6, pars), 4)
  expect_equal(DAISIE_A_n_splits(3, pars), 0)
  # closed form for large ranges without enumeration
  pars20 <- DAISIE_A_create_pars(n_islands = 20, lac = 1, mu = 1, gam = 1,
                                 laa = 1, vic_mode = "all")
  expect_equal(DAISIE_A_n_splits(1:20, pars20), 2^19 - 1)
})
