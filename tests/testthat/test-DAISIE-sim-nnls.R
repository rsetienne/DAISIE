library(nnls)
library(testthat)

test_that("NNLS properly recovers the true P_n distribution from P_rho", {
  # --- 1. Definieer de bekende parameters ---
  S <- 24        # Het aantal geobserveerde soorten in de clade
  
  # DEFINIEER HIER JE VARIABLE OM MEERDERE WAARDEN VAN n (max_n) TE TESTEN:
  test_max_n_waarden <- c(10, 20, 40, 80, 100)
  
  for (max_n in test_max_n_waarden) {
    cat("\n=== Testing voor maximaal", max_n, "missende soorten (n) ===\n")
    
    # --- 2. We definiëren een 'True' distributie voor P_n (de grondwaarheid!) ---
    # Stel dat de echte verdeling van missende soorten een Poisson verdeling is met gemiddelde 5:
    P_n_true <- dpois(0:max_n, lambda = 5)
    P_n_true <- P_n_true / sum(P_n_true) # Normaliseer zodat het exact 1 is
    
    # --- 3. We genereren wat we 'normaal' uit de ODE zouden krijgen (P_rho) ---
    rho_waarden <- seq(0.01, 1.0, length.out = max_n + 10) # Meer parameters vereisen meer datapunten!
    P_rho_gesimuleerd <- numeric(length(rho_waarden))
    
    for (i in seq_along(rho_waarden)) {
      rho <- rho_waarden[i]
      kans_som <- 0
      for (n in 0:max_n) {
        kans_som <- kans_som + P_n_true[n + 1] * choose(S + n, n) * (rho^S) * ((1 - rho)^n)
      }
      P_rho_gesimuleerd[i] <- kans_som
    }
    
    # --- 4. Nu doen we alsof we P_n_true NIET weten, en we P_n gaan 'schatten' via NNLS ---
    A <- matrix(0, nrow = length(rho_waarden), ncol = max_n + 1)
    for (i in seq_along(rho_waarden)) {
      rho <- rho_waarden[i]
      for (n in 0:max_n) {
        kolom <- n + 1
        A[i, kolom] <- choose(S + n, n) * (rho^S) * ((1 - rho)^n)
      }
    }
    
    fit <- nnls(A, P_rho_gesimuleerd)
    P_n_geschat <- fit$x
    
    # Normaliseren
    P_n_genormaliseerd <- P_n_geschat / sum(P_n_geschat)
    
    # Laten we testen of NNLS succesvol was! We stellen dynamisch de tolerantie in
    mse <- mean((P_n_genormaliseerd - P_n_true)^2)
    cat("Mean Squared Error van NNLS ten opzichte van ware P_n:", mse, "\n")
    
    expect_equal(P_n_genormaliseerd, P_n_true, tolerance = 0.01)
  }
})

