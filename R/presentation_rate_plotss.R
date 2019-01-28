# time <- 10
# M <- 1000
# pars <- c(0.0001, 2, 0.05, 0.001, 0.1)
# replicates <- 10
# mainland_params <- NULL
# divdepmodel <- "CS"
# island_ontogeny <- "quadratic"
# Apars <- create_area_params(5000, 0.1, 1, 11)
# Epars <- c(0.05, 100)
# rng_value <- 42
# set.seed(rng_value)
# ontogeny_CS_younger <- DAISIE_sim(time = time, M = M, pars = pars, 
#                                   replicates = replicates, divdepmodel = divdepmodel, Apars = Apars,
#                                   Epars = Epars, island_ontogeny = island_ontogeny, sample_freq = 50, plot_sims = FALSE, verbose = T)
# DAISIE_plot_sims(ontogeny_CS_younger)

# parameter plots
totaltime <- 15
lac <- 0.0001
laa <- 0.1
Epars <- c(0.01, 100)
K <- 0.05
gam <- 0.001
mainland_n <- 1000
Apars <- create_area_params(
  max_area = 5000,
  proportional_peak_t = 0.1,
  peak_sharpness = 1,
  total_island_age = 20
)
resolution <- 0.1
removed_timepoints <- 1
island_ontogeny <- "quadratic"

DAISIE_plot_area(
  totaltime = 15, 
  Apars = Apars,
  island_ontogeny = "quadratic",
  resolution = 0.100
)

DAISIE_plot_cladogenesis(
  totaltime = totaltime,
  K = K,
  lac = lac,
  Apars = Apars,
  island_ontogeny = island_ontogeny,
  removed_timepoints = removed_timepoints,
  resolution = resolution
)

DAISIE_plot_extinction(
  totaltime = totaltime,
  K = K,
  Apars = Apars,
  Epars = Epars,
  island_ontogeny = island_ontogeny,
  removed_timepoints = removed_timepoints,
  resolution = resolution
)

DAISIE_plot_immigration(totaltime = totaltime,
                        K = K,
                        Apars = Apars,
                        gam = gam,
                        mainland_n = mainland_n,
                        island_ontogeny = island_ontogeny,
                        removed_timepoints = removed_timepoints,
                        resolution = resolution)
