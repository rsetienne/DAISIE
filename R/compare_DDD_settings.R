compare_DDD_settings <- function(
    filename = 'simplex_odeint__runge_kutta_cash_karp54_timeout_1a9e5143686a8.rds',
    methode = 'odeint::runge_kutta_cash_karp54',
    optimmethod = 'simplex') {
  ff <- readRDS(filename)
  results <- dd_ML(brts = ff$brts,
                   initparsopt = ff$initpars,
                   methode = methode,
                   optimmethod = optimmethod)
  write.csv(c(filename, optimmethod, methode, initpars, results), file = 'output.csv', append = TRUE)
}
