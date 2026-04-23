#' @keywords internal
DAISIE_DE_solve_branch <- function(interval_func,
                                   initial_conditions,
                                   time,
                                   parameter,
                                   methode = "odeint::runge_kutta_cash_karp54",
                                   atol,
                                   rtol) {
  solution <- c()
  if (startsWith(methode, "odeint::")) {
    interval_name <- as.character(substitute(interval_func))
    if (interval_name == "interval_func") {
      interval_name <- interval_func # got passed as string
    }
    solution <- solve_branch_cpp(interval_name,
                                 initial_conditions,
                                 time,
                                 parameter,
                                 methode,
                                 atol,
                                 rtol)
  } else {
    solution <- deSolve::ode(
      y = initial_conditions,
      times = time,
      func = interval_func,
      parms = parameter,
      method = methode,
      atol = atol,
      rtol = rtol
    )
  }
  return(solution)
}

#' @keywords internal
solve_branch_cpp <- function(chosen_func,
                             initial_conditions,
                             time,
                             parameter,
                             methode = "odeint::runge_kutta_cash_karp54",
                             atol = 1e-15,
                             rtol = 1e-15) {
  lambda_c <- parameter[[1]]
  mu_E     <- parameter[[2]]
  mu_NE    <- parameter[[3]]
  gamma    <- parameter[[4]]
  lambda_a <- parameter[[5]]

  solution <- .Call("DAISIE_DE_cpp_solve",
                    lambda_c,
                    lambda_a,
                    mu_E,
                    mu_NE,
                    gamma,
                    chosen_func,
                    methode,
                    initial_conditions,
                    time,
                    atol,
                    rtol)
  if (length(time) == 2) {
    res <- matrix(data = NA, nrow = 2, ncol = length(solution$states))
    res[1, ] <- initial_conditions
    res[2, ] <- solution$states
    colnames(res) <- names(initial_conditions)
    return(res)
  } else if (length(time) > 2) {
    res <- matrix(data = NA, nrow = length(time), ncol = length(initial_conditions))
    for (i in 1:length(solution$states)) {
      res[i, ] <- solution$states[[i]]
    }
    res <- cbind(time, res)
    colnames(res) <- c("time", names(initial_conditions))
    return(res)
  }
}
