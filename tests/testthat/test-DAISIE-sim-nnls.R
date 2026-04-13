library(nnls)
library(testthat)

test_that("NNLS properly recovers the true P_n distribution for different missnumspec cases", {
  # --- 1. Definieer de bekende parameters ---
  S <- 10        # Het aantal geobserveerde soorten in de clade
  
  # DEFINIEER HIER JE VARIABLE OM VERSCHILLENDE MISSNUMSPEC (WARE MISSENDE SOORTEN) TE TESTEN:
  test_missnumspec_waarden <- c(2, 5, 12, 20)
  
  for (missnumspec in test_missnumspec_waarden) {
    cat("\n=== Testing voor missnumspec =", missnumspec, "===\n")
    
    # Zorg dat onze integratie limiet (max_n) altijd ruim groter is dan het echte aantal
    max_n <- missnumspec * 2 + 10
    
    # --- 2. We definiëren een 'True' distributie voor P_n (de grondwaarheid!) ---
    # In deze test bootsen we de situatie na waarin er EXACT 'missnumspec' soorten missen
    P_n_true <- numeric(max_n + 1)
    P_n_true[missnumspec + 1] <- 1 # Alleen de kans op n=missnumspec is 1, de rest 0
    
    # --- 3. We genereren wat we 'normaal' uit de ODE zouden krijgen (P_rho) ---
    rho_waarden <- seq(0.01, 1.0, length.out = max_n + 15) # Flexibele dichtheid
    P_rho_gesimuleerd <- numeric(length(rho_waarden))
    
    for (i in seq_along(rho_waarden)) {
      rho <- rho_waarden[i]
      # Aangezien we weten dat P_n_true alleen op 'missnumspec' staat, is in principe
      # de P_rho simpelweg deze vergelijking. Maar de loop houdt het netjes wiskundig.
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
    
    # Laten we testen of NNLS succesvol was! 
    mse <- mean((P_n_genormaliseerd - P_n_true)^2)
    cat("Mean Squared Error van NNLS ten opzichte van ware missnumspec:", mse, "\n")
    
    expect_equal(P_n_genormaliseerd, P_n_true, tolerance = 0.01)
  }
})


