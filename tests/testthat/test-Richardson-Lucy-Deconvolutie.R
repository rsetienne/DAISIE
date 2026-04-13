library(testthat)

test_that("RL properly recovers the true P_n distribution from P_rho", {
  richardson_lucy <- function(A, observed, iterations = 1000, tolerance = 1e-8) {
    loglik_old <- -Inf
    eps <- 1e-12
    observed <- pmax(observed, eps)
    n_bins <- ncol(A)

    # Pre-calculate de getransponeerde matrix en kolom-sommen voor snelheid
    At <- t(A)
    col_sums <- colSums(A)
    col_sums <- pmax(col_sums, eps)

    # x <- rep(1 / n_bins, n_bins) is eigenlijk ook normalisatie
    x <- rep(sum(observed) / sum(col_sums), n_bins)

    for (i in 1:iterations) {
      x_old <- x

      # 2. Projecteer de huidige gok naar de geobserveerde ruimte
      reprojection <- as.vector(A %*% x)

      # 3. Bereken de ratio (vermijd delen door nul)
      reprojection <- pmax(reprojection, eps)
      ratio <- observed / reprojection

      # 4. Bereken de correctiefactor (back-projection)
      correction <- as.vector(At %*% ratio)

      # We corrigeren voor de schaal van A via col_sums en maken die niet 0
      x <- x * (correction / col_sums)

      # Stop als de verbetering minimaal is
      # if (mean((x - x_old)^2) < tolerance) break werkt wel
      # maar neemt niet de algemene fluctuatie mee.
      if (any(!is.finite(x))) break
      loglik <- sum(observed * log(reprojection) - reprojection)

      if (abs(loglik - loglik_old) < tolerance) break
      loglik_old <- loglik
      # if (sum((x - x_old)^2) / sum(x_old^2) < tol) break kan natuurlijk ook
    }

    return(x)
  }

  # --- 1. Definieer de bekende parameters ---
  S <- 40

  # Meerdere waarden testen (bijv. n=5 en n=10)
  test_max_n_waarden <- c(1:10, 20, 30, 40)

  for (max_n in test_max_n_waarden) {
    cat("\n=== Testing voor maximaal", max_n, "missende soorten (n) ===\n")

    # --- 2. 'True' distributie (geen normalisatie) ---
    P_n_true <- dpois(0:max_n, lambda = 5)

    # --- 3. Genereren van P_rho (Simulatie uit ODE-vervanger) ---
    rho_waarden <- sort(c(
      seq(0.001, 0.1, length.out = max_n),
      seq(0.1, 0.9, length.out = max_n),
      seq(0.9, 0.999, length.out = max_n)
    ))
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

    P_n_geschat <- richardson_lucy(A, P_rho_gesimuleerd, iterations = 2000)

    # --- 5. Validatie ---
    mse <- mean((P_n_geschat - P_n_true)^2)
    cat("Mean Squared Error voor n =", max_n, ":", mse, "\n")

    expect_equal(P_n_geschat, P_n_true, tolerance = 1e-1)
  }
})
