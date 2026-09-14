#' branch solving for TRAISIE
#' @description
#' solve along branch
#' @inheritParams default_params_doc
#' @param interval_func chosen function for interval, can also be string if
#' using Rcpp
#' @param initial_conditions vector of initial conditions
#' @param time vector with two time points
#' @export
TRAISIE_solve_branch <- function(interval_func,
                                 initial_conditions,
                                 time,
                                 parameter,
                                 trait_mainland_ancestor,
                                 methode =  "odeint::runge_kutta_cash_karp54",
                                 atol,
                                 rtol) {
  solution <- c()
  if (startsWith(methode, "odeint::")) {
    interval_name <- as.character(substitute(interval_func))
    if (interval_name == "interval_func") {
      interval_name <- interval_func # got passed as string
    }
    solution <- TRAISIE_solve_branch_cpp(interval_name,
                                         initial_conditions,
                                         time,
                                         parameter,
                                         trait_mainland_ancestor,
                                         methode,
                                         atol,
                                         rtol)
  } else {
    parameter[[7]] <- trait_mainland_ancestor
    solution <- deSolve::ode(
      y = initial_conditions,
      times = time,
      func = interval_func,
      parms = parameter,
      method = methode,
      atol = atol,
      rtol = rtol
    )
    solution <- matrix(solution[, -1], nrow = 2)
  }

  return(solution)
}

#' @keywords internal
TRAISIE_solve_branch_cpp <- function(chosen_func,
                                     initial_conditions,
                                     time,
                                     parameter,
                                     trait_mainland_ancestor,
                                     methode = "odeint::runge_kutta_cash_karp54",
                                     atol = 1e-15,
                                     rtol = 1e-15) {

  lambda_c <- parameter[[1]]
  mus      <- parameter[[2]]
  gammas   <- parameter[[3]]
  lambda_a <- parameter[[4]]
  q_matrix       <- parameter[[5]]
  p_value       <- parameter[[6]]

  solution <- .Call("TRAISIE_branch_cpp",
                    lambda_c,
                    lambda_a,
                    mus,
                    gammas,
                    q_matrix,
                    p_value,
                    trait_mainland_ancestor,
                    chosen_func,
                    methode,
                    initial_conditions,
                    time,
                    atol,
                    rtol)

  res <- matrix(data = NA, nrow = 2, ncol = length(solution$states))
  res[1, ] <- initial_conditions
  res[2, ] <- solution$states

  return(res)
}

