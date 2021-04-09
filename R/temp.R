## -----------------------------------------------------------------------------
n_replicates <- 1000

## -----------------------------------------------------------------------------
set.seed(42)
island_age <- 5
n_mainland_species <- 1000
sample_counts <- c()
start_time <- Sys.time()
for (i in 1:100) {
  test_pars <- DAISIEutils::setup_model("cr_dd")
  pars <- test_pars$initparsopt

  sample_counts[i] <- DAISIE_sim_constant_rate(
    time = island_age,
    M = n_mainland_species,
    pars = pars,
    replicates = n_replicates,
    plot_sims = FALSE,
    verbose = FALSE
  )$sample_count
  message(sample_counts[i])
  message("iteration: ", i)
}

end_time <- Sys.time()

end_time - start_time
