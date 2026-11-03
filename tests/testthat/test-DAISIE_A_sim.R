sim_configs <- function() {
  star <- matrix(0, 4, 4)
  star[1, 2:4] <- 1
  star[2:4, 1] <- 1
  list(
    list(n = 1, conn = "full", k = 1, gam = 0.5, K = 3, cs = "split"),
    list(n = 3, conn = "line", k = 2, gam = 0.3, K = Inf, cs = "bud"),
    list(n = 4, conn = "full", k = 1, gam = 0.3, K = 3, cs = "split"),
    list(n = 4, conn = "ring", k = 3, gam = c(0.3, 0, 0.2, 0), K = Inf,
         cs = "bud", ages = c(3, 2, 3, 1.5), sub = c(NA, 0.5, NA, NA)),
    list(n = 4, conn = star, k = 2, gam = c(0.4, 0.1, 0, 0.1), K = Inf,
         cs = "bud"),
    list(n = 5, conn = "full", k = 4, gam = 0.3, K = 5, cs = "bud")
  )
}
sim_pars <- function(x) {
  DAISIE_A_create_pars(n_islands = x$n, lac = 0.3, mu = 0.6, K = x$K,
                       gam = x$gam, alpha = 1.2, k = x$k, disp = 0.8,
                       connectivity = x$conn, clado_single = x$cs,
                       island_ages = x$ages, submergence_ages = x$sub)
}

test_that("the C++ and R simulations are identical for the same seed", {
  for (x in sim_configs()) {
    p <- sim_pars(x)
    for (prot in c(FALSE, TRUE)) {
      for (mn in c(1, 2)) {
        for (s in 1:3) {
          set.seed(s)
          a <- DAISIE_A_sim_core_R(time = 3, mainland_n = mn, pars = p,
                                   protracted = prot)
          set.seed(s)
          b <- DAISIE_A_sim_core_cpp(time = 3, mainland_n = mn, pars = p,
                                     protracted = prot)
          expect_identical(a$spec, b$spec)
          expect_identical(unname(a$gen), unname(b$gen))
          expect_identical(unname(a$stt_table), unname(b$stt_table))
          expect_identical(unname(a$stt_islands), unname(b$stt_islands))
          expect_equal(a$n_events, b$n_events)
        }
      }
    }
  }
  p <- sim_pars(sim_configs()[[2]])
  set.seed(4)
  a <- DAISIE_A_sim(time = 3, M = 10, pars = p, replicates = 2,
                    verbose = FALSE, use_cpp = FALSE)
  set.seed(4)
  b <- DAISIE_A_sim(time = 3, M = 10, pars = p, replicates = 2,
                    verbose = FALSE)
  expect_identical(a, b)
})

test_that("the state of the simulation stays consistent", {
  withr::local_options(DAISIE_A.debug = TRUE)
  for (x in sim_configs()) {
    p <- sim_pars(x)
    for (prot in c(FALSE, TRUE)) {
      for (s in 1:4) {
        set.seed(s)
        expect_silent(DAISIE_A_sim_core_R(time = 3, mainland_n = 1, pars = p,
                                          protracted = prot))
      }
    }
  }
})

test_that("archipelago output is internally consistent", {
  pars <- DAISIE_A_create_pars(
    n_islands = 3, lac = 0.8, mu = 0.3, K = 6, gam = 0.05, alpha = 0.6,
    k = 2, disp = 0.4, connectivity = "line"
  )
  for (prot in c(FALSE, TRUE)) {
    set.seed(3)
    sims <- DAISIE_A_sim(time = 4, M = 20, pars = pars, replicates = 2,
                         verbose = FALSE, protracted = prot)
    for (rep in sims) {
      expect_equal(rep[[1]]$n_islands, 3)
      expect_length(rep[[1]]$stt_islands, 3)
      if (length(rep) == 1) next
      sp <- do.call(rbind, lapply(rep[-1], function(x) x$species))
      occ <- do.call(rbind, lapply(strsplit(sp$range, ""), as.numeric))
      expect_true(all(rowSums(occ) >= 1))
      expect_equal(rowSums(occ), sp$n_islands)
      for (i in 1:3) {
        fin <- rep[[1]]$stt_islands[[i]]
        fin <- fin[nrow(fin), ]
        expect_equal(unname(fin["nI"]), sum(occ[, i] & sp$type == "I"))
        expect_equal(unname(fin["nA"]), sum(occ[, i] & sp$type == "A"))
        expect_equal(unname(fin["nC"]), sum(occ[, i] & sp$type == "C"))
        expect_equal(unname(fin["nSIE"]),
                     sum(occ[, i] & sp$type != "I" & rowSums(occ) == 1))
      }
      stt_all <- rep[[1]]$stt_all
      expect_equal(unname(stt_all[nrow(stt_all), c("nI", "nA", "nC")]),
                   c(sum(sp$type == "I"), sum(sp$type == "A"),
                     sum(sp$type == "C")))
      for (x in rep[-1]) {
        o <- do.call(rbind, lapply(strsplit(x$species$range, ""), as.numeric))
        expect_true(all(colSums(o) <= 6))
        tips <- unlist(lapply(x$clades, function(cl) cl$tips$species))
        expect_setequal(tips, x$species$species)
        # at most one non-endemic lineage per mainland species
        expect_true(sum(x$species$type == "I") <= 1)
        for (cl in x$clades) {
          expect_equal(cl$branching_times[1], 4)
          expect_equal(length(cl$branching_times), nrow(cl$tips) + 1)
          nd <- cl$nodes
          # nodes are not older than their parents, the crown is not older
          # than the colonisation, and with the protracted option nodes are
          # dated at or before the split
          int <- nd[is.na(nd$species), ]
          par_age <- nd$age[match(int$parent, nd$node)]
          expect_true(all(is.na(par_age) | int$age <= par_age + 1e-12))
          expect_true(all(nd$age <= cl$colonisation_time + 1e-12))
          expect_true(all(nd$age >= nd$age_split - 1e-12))
          if (!prot) expect_equal(nd$age, nd$age_split)
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
  }
})

test_that("islands cannot be occupied before they emerge", {
  pars <- DAISIE_A_create_pars(
    n_islands = 3, lac = 0.5, mu = 0.2, gam = 0.1, alpha = 0.3, k = 2,
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
  # no dispersal, divergence or extinction: the non-endemic can only spread
  # by colonisation and bud
  pars <- DAISIE_A_create_pars(n_islands = 2, lac = 0.5, mu = 0, gam = 1,
                               alpha = 0, k = 1, disp = 0)
  set.seed(1)
  sims <- DAISIE_A_sim(time = 2, M = 1, pars = pars, replicates = 3,
                       verbose = FALSE)
  for (rep in sims) {
    if (length(rep) < 2) next
    sp <- rep[[2]]$species
    expect_equal(sum(sp$type == "I"), 1)
    expect_true(all(sp$n_islands[sp$type != "I"] == 1))
  }
})

test_that("a submerged island loses all its species", {
  pars <- DAISIE_A_create_pars(3, lac = 0.6, mu = 0.2, gam = 0.1, alpha = 0.3,
                               k = 2, disp = 0.5, connectivity = "line",
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
  p0 <- DAISIE_A_create_pars(2, lac = 0.6, mu = 0.2, gam = 0.1, alpha = 0.3,
                             disp = 0.5)
  p1 <- DAISIE_A_create_pars(2, lac = 0.6, mu = 0.2, gam = 0.1, alpha = 0.3,
                             disp = 0.5, submergence_ages = c(NA, NA))
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

test_that("with one island and k = 1 the model is DAISIE", {
  skip_on_cran()
  # anagenesis at rate alpha, DAISIE's cladogenesis with clado_single = "split"
  summ <- function(sims) {
    vapply(sims, function(r) {
      s <- r[[1]]$stt_all
      last <- s[nrow(s), ]
      paste(last[["nI"]], last[["nA"]], min(last[["nC"]], 4),
            if (length(r) > 1) r[[2]]$stac else 0)
    }, "")
  }
  for (K in c(3, Inf)) {
    set.seed(1)
    d <- DAISIE_sim_cr(time = 2.5, M = 1, pars = c(0.8, 0.5, K, 0.4, 0.6),
                       replicates = 4000, divdepmodel = "CS",
                       plot_sims = FALSE, verbose = FALSE)
    pA <- DAISIE_A_create_pars(n_islands = 1, lac = 0.8, mu = 0.5, K = K,
                               gam = 0.4, alpha = 0.6, k = 1,
                               clado_single = "split")
    set.seed(2)
    a <- DAISIE_A_sim(time = 2.5, M = 1, pars = pA, replicates = 4000,
                      verbose = FALSE)
    kd <- summ(d)
    ka <- summ(a)
    lev <- union(kd, ka)
    td <- table(factor(kd, lev))
    ta <- table(factor(ka, lev))
    keep <- (td + ta) >= 20
    tab <- rbind(c(td[keep], sum(td[!keep])), c(ta[keep], sum(ta[!keep])))
    expect_gt(suppressWarnings(stats::chisq.test(tab))$p.value, 1e-3)
  }
})

test_that("DAISIE_A_plot_stt_islands runs", {
  pars <- DAISIE_A_create_pars(n_islands = 2, lac = 0.5, mu = 0.2,
                               gam = 0.05, alpha = 0.5, k = 2, disp = 0.3)
  set.seed(1)
  sims <- DAISIE_A_sim(time = 2, M = 20, pars = pars, replicates = 3,
                       verbose = FALSE)
  grDevices::pdf(NULL)
  expect_silent(q <- DAISIE_A_plot_stt_islands(sims, what = "endemic"))
  grDevices::dev.off()
  expect_length(q, 2)
})
