# Independent check of the simulation: the distribution of the final state of
# one mainland species, computed exactly from the master equation of the
# model (helper-DAISIE_A_exact.R), against the frequencies in simulations.

test_that("the exact chain of one island with K = 1 has 3 states", {
  # empty -> I (gam), I -> A (alpha: the mainland pair reaches k = 1),
  # I -> empty and A -> empty (mu); recolonization of I changes no state and
  # nothing else can happen at K = 1
  gam <- 0.7
  alpha <- 0.4
  mu <- 0.3
  pars <- DAISIE_A_create_pars(n_islands = 1, lac = 0.5, mu = mu, K = 1,
                               gam = gam, alpha = alpha, k = 1)
  chain <- exact_chain(pars)
  expect_equal(nrow(chain$Q), 3)
  expect_equal(unname(Matrix::rowSums(chain$Q)), rep(0, 3))
  summ <- vapply(chain$states, exact_summary, "")
  expect_setequal(summ, c("", "I1", "A1"))
  Q <- matrix(c(-gam, gam, 0,
                mu, -(mu + alpha), alpha,
                mu, 0, -mu), 3, 3, byrow = TRUE)
  t <- 1.7
  ref <- as.numeric(as.matrix(Matrix::expm(Matrix::Matrix(Q * t)))[1, ])
  p <- exact_distribution(chain, t)
  expect_equal(p[match(c("", "I1", "A1"), summ)], ref, tolerance = 1e-10)
})

test_that("simulated final states agree with the exact distribution", {
  skip_on_cran()
  configs <- list(
    # two adjacent islands, both reached from the mainland
    list(t = 2.5, pars = DAISIE_A_create_pars(
      n_islands = 2, lac = 0.5, mu = 0.3, K = 2, gam = 0.4, alpha = 0.7,
      k = 2, disp = 0.5, clado_single = "split")),
    # the mainland reaches only island 1: gap pairs with the mainland
    list(t = 3, pars = DAISIE_A_create_pars(
      n_islands = 2, lac = 0.5, mu = 0.3, K = 2, gam = c(0.5, 0),
      alpha = 0.7, k = 2, disp = 0.6, clado_single = "bud")),
    # three islands on a line: gap pairs between islands 1 and 3
    list(t = 3, pars = DAISIE_A_create_pars(
      n_islands = 3, lac = 0.6, mu = 0.3, K = c(1, 2, 1),
      gam = c(0.4, 0, 0.4), alpha = 0.8, k = 2, disp = 0.6,
      connectivity = "line", clado_single = "split")),
    # k = 1: every step completes a divergence
    list(t = 2, pars = DAISIE_A_create_pars(
      n_islands = 2, lac = 0.5, mu = 0.3, K = 3, gam = 0.3, alpha = 0.6,
      k = 1, disp = 0.5, clado_single = "split"))
  )
  n_sims <- 30000
  for (q in seq_along(configs)) {
    x <- configs[[q]]
    chain <- exact_chain(x$pars)
    p <- exact_distribution(chain, x$t)
    expect_equal(sum(p), 1, tolerance = 1e-10)
    expected <- tapply(p, vapply(chain$states, exact_summary, ""), sum)
    cfg <- DAISIE_A_core_cfg(x$pars, x$t, protracted = FALSE)
    cfg$mainland_n <- 1L
    set.seed(q)
    sims <- vapply(seq_len(n_sims), function(i) {
      exact_sim_summary(DAISIE_A_call_cpp(cfg), x$pars$n_islands)
    }, "")
    # every simulated configuration is possible in the exact model
    expect_true(all(sims %in% names(expected)))
    observed <- table(factor(sims, levels = names(expected)))
    e <- n_sims * expected
    keep <- e >= 5
    o <- c(observed[keep], sum(observed[!keep]))
    e <- c(e[keep], sum(e[!keep]))
    chi2 <- sum((o - e)^2 / e)
    expect_gt(stats::pchisq(chi2, df = length(o) - 1, lower.tail = FALSE),
              1e-3)
  }
})
