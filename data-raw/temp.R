lac0 <- 2.000
mu0 <- 2.700
K0 <- 20.000
gam0 <- 0.009
laa0 <- 1.010
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


for (i in 1:100) {

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
  loglik_out[i] <- DAISIE:::DAISIE_loglik_CS_choice(pars1 = pars1_time_dep[[i]],
                                           pars2 = pars2,
                                           brts = brts,
                                           stac = stac,
                                           missnumspec = missnumspec,
                                           CS_version = CS_version

  )
  print(i)
}

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

out <- c(loglik_out - loglik3)
plot(out, type = "l", ylab = "loglik difference", xlab = "Area")
