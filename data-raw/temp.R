# Use robustness parameters instead
# Try to find d and x with ML

lac0 <- 0.02
mu0 <- 0.975
K0 <- 0.001
gam0 <- 0.03363
laa0 <- 0.0295
d <- 0.1108
x <- 0.075
area_pars <- c()
peak <- c()
pars1_time_dep <- list()
total_time <- 2.55
loglik_out <- c()

island_ontogeny <- 1
sea_level <- 0


pars2 <- c(1.0e+02, 1.1e+01, 0.0e+00, 0.0e+00, NA, 0.0e+00, 1.0e-04,
           1.0e-05, 1.0e-07, 3.0e+03, 9.5e-01, 9.8e-01)
brts <- c(4.0000, 3.0282, 1.3227, 0.8223, 0.4286, 0.3462, 0.2450,
          0.0808, 0.0527, 0.0327, 0.0221, 0.1180, 0.0756, 0.0525,
          0.0322, 0.0118)


stac <- 2
missnumspec <- 0
CS_version <- 0

it <- 1
for (i in seq(1, 1000, by = 10)) {

  area_pars <- c(
    max_area = i,
    current_area = 0.99,
    proportional_peak_t = 0.50,
    total_island_age = 2.864,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0
  )

  area_pars_list <- DAISIE::create_area_pars(
    max_area = i,
    current_area = 0.99,
    proportional_peak_t = 0.50,
    total_island_age = 2.864,
    sea_level_amplitude = 0,
    sea_level_frequency = 0,
    island_gradient_angle = 0
  )
  peak <- DAISIE:::calc_peak(total_time = total_time, area_pars = area_pars_list)
  pars1_time_dep[[i]] <- c(
    lac0,
    mu0,
    K0,
    gam0,
    laa0,
    d,
    x,
    area_pars,
    island_ontogeny,
    sea_level,
    total_time,
    peak
  )
  loglik_out[it] <- DAISIE:::DAISIE_loglik_CS_choice(pars1 = pars1_time_dep[[i]],
                                           pars2 = pars2,
                                           brts = brts,
                                           stac = stac,
                                           missnumspec = missnumspec,
                                           CS_version = CS_version

  )
  print(it)
  it <- it + 1
}

loglik_out == sort(loglik_out) # Identical vectors, vector is sorted

# deSolve lsodes time dep function with Hawaii area
pars1_const_rate <- c(
  lac0,
  mu0,
  K0,
  gam0,
  laa0
)
# deSolve lsodes constant rate function
loglik3 <-
  DAISIE:::DAISIE_loglik_CS_choice(
    pars1 = pars1_const_rate,
    pars2 = pars2,
    brts = brts,
    stac = stac,
    missnumspec = missnumspec,
    CS_version = CS_version
  )

out <- loglik_out - loglik3
plot(out[which(!is.na(out))], type = "l", ylab = "loglik difference", xlab = "Area",)


area_pars_list <- DAISIE::create_area_pars(
  max_area = i,
  current_area = 0.99,
  proportional_peak_t = 0.50,
  total_island_age = 2.864,
  sea_level_amplitude = 0,
  sea_level_frequency = 0,
  island_gradient_angle = 0
)
peak <- DAISIE:::calc_peak(total_time = total_time, area_pars = area_pars_list)
DAISIE::DAISIE_plot_area(2.864, area_pars_list, peak, "beta", 0.0001, "const")
DAISIE_plot_extinction(
  total_time = 2.5,
  area_pars = area_pars_list,
  peak = peak,
  mu = mu0,
  hyper_pars =  DAISIE::create_hyper_pars(d = d, x = x),
  island_ontogeny = "beta",
  resolution = 0.0001,
  sea_level = "const",
  extcutoff = 1100,
  removed_timepoints = 1
)

DAISIE_plot_cladogenesis(
  total_time = total_time,
  K = 0.05,
  area_pars = area_pars_list,
  peak = peak,
  lac = lac0,
  island_ontogeny = "beta",
  sea_level = "const",
  hyper_pars = create_hyper_pars(d = d, x = x),
  1,
  0.001
)

timeval_vec <- seq(0, 2.864, by = 0.001)
area <- c()
for (i in 1:1000) {
  area[i] <- island_area(
    timeval = timeval_vec[i],
    total_time = total_time,
    area_pars = area_pars_list,
    peak = peak,
    island_ontogeny = island_ontogeny,
    sea_level = sea_level
  )
}
plot(area)
