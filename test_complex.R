data(Galapagos_datalist)
datalist <- Galapagos_datalist
brts <- datalist[[5]]$branching_times
missnumspec <- datalist[[5]]$missing_species

pars1 <- c(2.546591, 2.678781, 2.678781, 0.009326754, 1.008583)

stac = 2
methode = "odeint::runge_kutta_cash_karp54"
reltolint = 1e-15
abstolint = 1e-15

t0 <- brts[1]
t1 <- brts[2]
t2 <- brts[3]
tp <- 0
ti <- sort(brts)
ti <- ti[1:(length(ti) - 2)]

# Initial conditions
number_of_species <- length(brts) - 1
rho <- number_of_species / (missnumspec + number_of_species)

init_D <- 1 # originally rho

initial_conditions1   <- c(DE = init_D, DM3 = 0, E = 1 - rho, DA3 = 1)
if (stac == 3 || stac == 7) {
  initial_conditions1 <- c(DE = init_D, DM3 = 1, E = 1 - rho, DA3 = 0)
}

solution0 <- DAISIE:::DAISIE_DE_solve_branch(interval_func = interval2_EC,
                                    initial_conditions = initial_conditions1,
                                    time = c(0, ti),
                                    parameter = pars1,
                                    methode = methode,
                                    rtol = reltolint,
                                    atol = abstolint)

# and now with complex numbers
initial_conditions_complex <- c()
for (i in 1:length(initial_conditions1)) {
  initial_conditions_complex[i] <- as.complex(initial_conditions1[i])
}

names(initial_conditions_complex) <- names(initial_conditions1)

pars1_c <- c()
for (i in 1:length(pars1)) {
  pars1_c[i] <- as.complex(pars1[i])
}


solution0_c <- DAISIE:::DAISIE_DE_solve_branch(interval_func = interval2_EC,
                                             initial_conditions = initial_conditions_complex,
                                             time = c(0, ti),
                                             parameter = pars1,
                                             methode = methode,
                                             rtol = reltolint,
                                             atol = abstolint)

testthat::expect_equal(dim(solution0), dim(solution0_c))
for (i in 1:nrow(solution0)) {
  for (j in 1:ncol(solution0)) {
     a <- solution0[i, j]
     b <- solution0_c[i, j]
     testthat::expect_equal(a, Re(b))
  }
}



