#' branch solving
#' @description
#' solve along branch
#' @inheritParams default_params_doc
#' @param interval_func chosen function for interval, can also be string if using Rcpp
#' @param initial_conditions vector of initial conditions
#' @param time vector with two time points
#' @export
solve_branch <- function(interval_func,
                         initial_conditions,
                         time,
                         parameter,
                         methode = "ode45",
                         rcpp_methode = "odeint::bulirsch_stoer",
                         atol,
                         rtol,
                         use_Rcpp = TRUE) {
  solution <- c()
  if (use_Rcpp == FALSE) {
    solution <- deSolve::ode(
      y = initial_conditions,
      times = time,
      func = interval_func,
      parms = parameter,
      method = methode,
      atol = atol,
      rtol = rtol
    )
  } else {
    interval_name <- as.character(substitute(interval_func))
    solution <- solve_branch_cpp(interval_name,
                                 initial_conditions,
                                 time,
                                 parameter,
                                 rcpp_methode,
                                 atol,
                                 rtol)
  }
  return(solution)
}

#' @keywords internal
solve_branch_cpp <- function(chosen_func,
                             initial_conditions,
                             time,
                             parameter,
                             methode = "odeint::bulirsch_stoer",
                             atol = 1e-15,
                             rtol = 1e-15) {

  lambda_c <- parameter[[1]]
  mu      <- parameter[[2]]
  gamma   <- parameter[[4]]
  lambda_a <- parameter[[5]]

  solution <- cpp_solve(lambda_c,
                        lambda_a,
                        mu,
                        gamma,
                        chosen_func,
                        methode,
                        initial_conditions,
                        time,
                        atol,
                        rtol)

  res <- matrix(data = NA, nrow = 2, ncol = length(solution$states))
  res[1, ] <- initial_conditions
  res[2, ] <- solution$states
  colnames(res) <- names(initial_conditions)
  return(res)
}
