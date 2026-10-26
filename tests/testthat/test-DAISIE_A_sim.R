strip_archipelago <- function(sims) {
  lapply(sims, function(rep) {
    rep[[1]] <- rep[[1]][c("island_age", "not_present", "stt_all")]
    for (i in seq_along(rep)[-1]) {
      rep[[i]] <- rep[[i]][setdiff(names(rep[[i]]),
                                   c("species", "clades", "mixed_clade"))]
    }
    rep
  })
}

test_that("with one island DAISIE_A_sim reproduces DAISIE_sim_cr exactly", {
  parsets <- list(c(0.5, 0.2, Inf, 0.05, 0.5),
                  c(2, 1.5, 3, 1, 0.1),
                  c(1.5, 0.8, 5, 0.3, 0.6))
  for (p in parsets) {
    for (seed in 1:3) {
      for (sf in c(25, Inf)) {
        set.seed(seed)
        daisie <- DAISIE_sim_cr(time = 3, M = 10, pars = p, replicates = 2,
                                sample_freq = sf, plot_sims = FALSE,
                                verbose = FALSE)
        # island-level anagenesis does not apply to a single island
        pars_a <- DAISIE_A_create_pars(n_islands = 1, lac = p[1], mu = p[2],
                                       K = p[3], gam = p[4], laa = p[5],
                                       laa_isl = if (seed == 2) 0.7 else 0)
        set.seed(seed)
        arch <- DAISIE_A_sim(time = 3, M = 10, pars = pars_a, replicates = 2,
                             sample_freq = sf, verbose = FALSE)
        expect_identical(daisie, strip_archipelago(arch))
      }
    }
  }
})

test_that("archipelago output is internally consistent", {
  pars <- DAISIE_A_create_pars(
    n_islands = 3, lac = 0.8, mu = 0.3, K = 6, gam = 0.05, laa = 0.3,
    lav = 0.4, disp = 0.4, laa_isl = 0.2, connectivity = "line"
  )
  set.seed(3)
  sims <- DAISIE_A_sim(time = 4, M = 20, pars = pars, replicates = 2,
                       verbose = FALSE)
  for (rep in sims) {
    expect_equal(rep[[1]]$n_islands, 3)
    expect_length(rep[[1]]$stt_islands, 3)
    if (length(rep) == 1) next
    sp <- do.call(rbind, lapply(rep[-1], function(x) x$species))
    occ <- do.call(rbind, lapply(strsplit(sp$range, ""), as.numeric))
    expect_true(all(rowSums(occ) >= 1))
    expect_equal(rowSums(occ), sp$n_islands)
    # final per-island STT equals the species table
    for (i in 1:3) {
      fin <- rep[[1]]$stt_islands[[i]]
      fin <- fin[nrow(fin), ]
      expect_equal(unname(fin["nI"]), sum(occ[, i] & sp$type == "I"))
      expect_equal(unname(fin["nA"]), sum(occ[, i] & sp$type == "A"))
      expect_equal(unname(fin["nC"]), sum(occ[, i] & sp$type == "C"))
      expect_equal(unname(fin["nSIE"]),
                   sum(occ[, i] & sp$type != "I" & rowSums(occ) == 1))
    }
    # archipelago-level STT equals the species table
    stt_all <- rep[[1]]$stt_all
    expect_equal(unname(stt_all[nrow(stt_all), c("nI", "nA", "nC")]),
                 c(sum(sp$type == "I"), sum(sp$type == "A"),
                   sum(sp$type == "C")))
    for (x in rep[-1]) {
      # per island clade diversity never exceeds K
      o <- do.call(rbind, lapply(strsplit(x$species$range, ""), as.numeric))
      expect_true(all(colSums(o) <= 6))
      # every species belongs to exactly one clade tree
      tips <- unlist(lapply(x$clades, function(cl) cl$tips$species))
      expect_setequal(tips, x$species$species)
      for (cl in x$clades) {
        expect_equal(cl$branching_times[1], 4)
        expect_equal(length(cl$branching_times), nrow(cl$tips) + 1)
        if (requireNamespace("ape", quietly = TRUE) && nrow(cl$tips) > 1) {
          tr <- ape::read.tree(text = cl$newick)
          expect_setequal(tr$tip.label, as.character(cl$tips$species))
          expect_equal(sort(unname(ape::branching.times(tr)),
                            decreasing = TRUE),
                       cl$branching_times[-(1:2)], tolerance = 1e-8)
          expect_equal(max(ape::branching.times(tr)) + tr$root.edge,
                       cl$colonisation_time, tolerance = 1e-8)
        }
      }
    }
  }
})

test_that("range restrictions are respected", {
  pars <- DAISIE_A_create_pars(
    n_islands = 5, lac = 0.5, mu = 0.3, gam = 0.05, laa = 0.3, lav = 0.3,
    disp = 0.8, connectivity = "line", connected_ranges = TRUE,
    max_range_size = 3, vic_mode = "all"
  )
  set.seed(5)
  sims <- DAISIE_A_sim(time = 4, M = 20, pars = pars, replicates = 2,
                       verbose = FALSE)
  adj <- pars$adjacency
  for (rep in sims) {
    for (x in rep[-1]) {
      for (r in x$species$range) {
        set <- which(strsplit(r, "")[[1]] == "1")
        expect_true(length(set) <= 3)
        expect_true(DAISIE_A_is_connected(set, adj))
      }
    }
  }
})

test_that("islands cannot be occupied before they emerge", {
  pars <- DAISIE_A_create_pars(
    n_islands = 3, lac = 0.5, mu = 0.2, gam = 0.1, laa = 0.3, lav = 0.3,
    disp = 0.5, connectivity = "line", island_ages = c(3, 2, 1)
  )
  set.seed(2)
  sims <- DAISIE_A_sim(time = 3, M = 20, pars = pars, replicates = 2,
                       sample_freq = Inf, verbose = FALSE)
  for (rep in sims) {
    for (i in 1:3) {
      m <- rep[[1]]$stt_islands[[i]]
      before <- m[, "Time"] > pars$island_ages[i]
      expect_true(all(m[before, c("nI", "nA", "nC")] == 0))
    }
  }
  expect_error(DAISIE_A_sim(time = 2, M = 2, pars = pars, replicates = 1,
                            verbose = FALSE), "Island ages")
})

test_that("budding keeps the identity of a widespread parent", {
  # no dispersal, vicariance or extinction; a lineage on two islands created
  # by recolonisation can only bud
  pars <- DAISIE_A_create_pars(
    n_islands = 2, lac = 0.5, mu = 0, gam = 1, laa = 0, lav = 0, disp = 0,
    clado_single = "bud"
  )
  set.seed(1)
  sims <- DAISIE_A_sim(time = 2, M = 1, pars = pars, replicates = 3,
                       verbose = FALSE)
  for (rep in sims) {
    sp <- rep[[2]]$species
    # exactly one non-endemic lineage: the mainland species itself
    expect_equal(sum(sp$type == "I"), 1)
    expect_true(all(sp$n_islands[sp$type != "I"] == 1))
  }
})

test_that("island-level anagenesis replaces one island population", {
  pars <- DAISIE_A_create_pars(n_islands = 3, lac = 0, mu = 0, gam = 0.1,
                               laa = 0, laa_isl = c(0, 1, 0.5))
  st <- DAISIE_A_init_state(n_islands = 3, mainland_n = 1)
  st$n_clades <- 1L
  st <- DAISIE_A_add_node(st, parent = 0, time = 0.5, kind = 1, side = 0,
                          clade = 1, species = 1)
  # the non-endemic mainland species on islands 1 and 2
  st$spec <- rbind(st$spec, c(1, 1, 0.5, 1, 0, st$n_gen, 0, 1, 1, 1, 0))
  # a widespread endemic on islands 2 and 3 (no island-level anagenesis:
  # for endemics this event is vicariance)
  st$spec <- rbind(st$spec, c(2, 1, 0.5, 2, 1, st$n_gen, 0, 1, 0, 1, 1))
  st$maxspecID <- 2
  rates <- DAISIE_A_update_rates(st, timeval = 1, pars = pars,
                                 emerge = rep(0, 3), mainland_spec = 1,
                                 split_cache = new.env())
  # only island 2 of the range of the non-endemic has a positive rate
  expect_equal(rates$isl_ana, 1)
  expect_equal(sum(rates$Wi[2, ]), 0)
  st2 <- DAISIE_A_update_state(st, possible_event = 7, rates = rates,
                               timeval = 1, pars = pars,
                               split_cache = new.env())
  occ <- st2$spec[, .DA_NFIX + 1:3]
  expect_equal(nrow(st2$spec), 3)
  # the parent stays non-endemic, keeps its ID and loses island 2
  expect_equal(st2$spec[1, c(.DA_ID, .DA_TYPE)], c(1, 1))
  expect_equal(unname(occ[1, ]), c(1, 0, 0))
  # the endemic is unchanged
  expect_equal(unname(occ[2, ]), c(0, 1, 1))
  # the new endemic occupies island 2 only
  expect_equal(st2$spec[3, .DA_TYPE], 3)
  expect_equal(unname(occ[3, ]), c(0, 1, 0))
  # the number of species on every island is unchanged
  expect_equal(unname(colSums(occ)), c(1, 2, 1))
  # a lineage on one island has no island-level anagenesis
  rates2 <- DAISIE_A_update_rates(st2, timeval = 1, pars = pars,
                                  emerge = rep(0, 3), mainland_spec = 1,
                                  split_cache = new.env())
  expect_equal(rates2$isl_ana, 0)
  # in a simulation, endemics are only formed on single islands and the
  # mainland species is never replaced
  pars <- DAISIE_A_create_pars(n_islands = 2, lac = 0, mu = 0, gam = 1,
                               laa = 0, laa_isl = 0.5)
  set.seed(1)
  sims <- DAISIE_A_sim(time = 2, M = 1, pars = pars, replicates = 3,
                       verbose = FALSE)
  for (rep in sims) {
    sp <- rep[[2]]$species
    expect_equal(sum(sp$type == "I"), 1)
    expect_true(all(sp$n_islands[sp$type != "I"] == 1))
  }
})

test_that("a submerged island loses all its species", {
  pars <- DAISIE_A_create_pars(3, lac = 0.6, mu = 0.2, gam = 0.1, laa = 0.3,
                               lav = 0.3, disp = 0.5, laa_isl = 0.1,
                               connectivity = "line",
                               submergence_ages = c(NA, NA, 0.7))
  set.seed(9)
  sims <- DAISIE_A_sim(time = 2.5, M = 20, pars = pars, replicates = 3,
                       sample_freq = Inf, verbose = FALSE)
  for (rep in sims) {
    m <- rep[[1]]$stt_islands[[3]]
    after <- m[, "Time"] <= 0.7
    expect_true(all(m[after, c("nI", "nA", "nC")] == 0))
    expect_true(any(m[!after, "nI"] + m[!after, "nA"] + m[!after, "nC"] > 0))
    for (x in rep[-1]) {
      expect_true(all(substr(x$species$range, 3, 3) == "0"))
    }
  }
  # all-NA submergence ages give the same simulation as none
  p0 <- DAISIE_A_create_pars(2, lac = 0.6, mu = 0.2, gam = 0.1, laa = 0.3,
                             lav = 0.3, disp = 0.5)
  p1 <- DAISIE_A_create_pars(2, lac = 0.6, mu = 0.2, gam = 0.1, laa = 0.3,
                             lav = 0.3, disp = 0.5,
                             submergence_ages = c(NA, NA))
  set.seed(1)
  a <- DAISIE_A_sim(time = 2, M = 10, pars = p0, replicates = 1,
                    verbose = FALSE)
  set.seed(1)
  b <- DAISIE_A_sim(time = 2, M = 10, pars = p1, replicates = 1,
                    verbose = FALSE)
  expect_identical(a, b)
  expect_error(DAISIE_A_sim(time = 0.5, M = 2, pars = pars, replicates = 1,
                            verbose = FALSE), "Submergence")
})

test_that("DAISIE_A_plot_stt_islands runs", {
  pars <- DAISIE_A_create_pars(n_islands = 2, lac = 0.5, mu = 0.2,
                               gam = 0.05, laa = 0.5, lav = 0.2, disp = 0.3)
  set.seed(1)
  sims <- DAISIE_A_sim(time = 2, M = 20, pars = pars, replicates = 3,
                       verbose = FALSE)
  grDevices::pdf(NULL)
  expect_silent(q <- DAISIE_A_plot_stt_islands(sims, what = "endemic"))
  grDevices::dev.off()
  expect_length(q, 2)
})
