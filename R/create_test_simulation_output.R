create_test_simulation_output <- function() {
  pars_equal = c(2.550687345,2.683454548,Inf,0.00933207,1.010073119)
  DAISIE::DAISIE_sim(
    time = 4,
    M = 1000,
    pars = pars_equal,
    replicates = 1
  )
  
}
