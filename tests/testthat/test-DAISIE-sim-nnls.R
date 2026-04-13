library(nnls)
library(testthat)

test_that("NNLS properly recovers the true P_n distribution from P_rho", {
  # --- 1. Definieer de bekende parameters ---
  S <- 40

  # Meerdere waarden testen (bijv. n=5 en n=10)
  test_max_n_waarden <- c(5, 10, 20, 30, 40)

  for (max_n in test_max_n_waarden) {
    cat("\n=== Testing voor maximaal", max_n, "missende soorten (n) ===\n")

    # --- 2. 'True' distributie (geen normalisatie) ---
    P_n_true <- dpois(0:max_n, lambda = 5)

    # --- 3. Genereren van P_rho (Simulatie uit ODE-vervanger) ---
    rho_waarden <- seq(0.01, 1.0, length.out = max_n + 10)
    P_rho_gesimuleerd <- numeric(length(rho_waarden))

    for (i in seq_along(rho_waarden)) {
      rho <- rho_waarden[i]
      kans_som <- 0
      for (n in 0:max_n) {
        # De originele formule
        kans_som <- kans_som + P_n_true[n + 1] * choose(S + n, n) * (rho^S) * ((1 - rho)^n)
      }
      P_rho_gesimuleerd[i] <- kans_som
    }

    # --- 4. Bouwen van matrix A en de NNLS fit ---
    A <- matrix(0, nrow = length(rho_waarden), ncol = max_n + 1)
    for (i in seq_along(rho_waarden)) {
      rho <- rho_waarden[i]
      for (n in 0:max_n) {
        kolom <- n + 1
        A[i, kolom] <- choose(S + n, n) * (rho^S) * ((1 - rho)^n)
      }
    }

    fit <- nnls(A, P_rho_gesimuleerd)

    # Gebruik de ruwe resultaten direct binnen de loop
    P_n_geschat <- fit$x

    # --- 5. Validatie ---
    mse <- mean((P_n_geschat - P_n_true)^2)
    cat("Mean Squared Error voor n =", max_n, ":", mse, "\n")

    expect_equal(P_n_geschat, P_n_true, tolerance = 0.1)
  }
})
