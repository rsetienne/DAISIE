#' Maximum Likelihood Estimation for the DAISIE Trait Model
#'
#' @inheritParams default_params_doc
#' @description
#' This function estimates the maximum likelihood (ML) for a given dataset
#' using a trait-dependent model. The model includes parameters for
#' cladogenesis, extinction, colonization, anagenesis, and state transitions.
#' The optimization is done using the `TRAISIE::calc_ml` function.
#' @param datalist A list containing the data used for likelihood calculation.
#' @param idparslist A list containing the model parameters required for the
#' trait-dependent diversification process.
#' Each element in the list corresponds to a specific set of parameters:
#'
#' \describe{
#'   \item{\code{idparslist[[1]]}}{
#'     A numeric vector specifying the \strong{cladogenesis rates}
#'     (\eqn{\lambda_c}) for each trait state. Each element in the vector
#'     corresponds to the cladogenesis rate for a different state.
#'   }
#'
#'   \item{\code{idparslist[[2]]}}{
#'     A numeric vector specifying the \strong{extinction rates}
#'     (\eqn{\mu}) for each trait state. Each element in the vector
#'     corresponds to the extinction rate for a different state.
#'   }
#'
#'   \item{\code{idparslist[[3]]}}{
#'     A numeric vector specifying the \strong{immigration rates}
#'     (\eqn{\gamma}) for each trait state. Each element in the vector
#'     corresponds to the immigration rate for a different state.
#'   }
#'
#'   \item{\code{idparslist[[4]]}}{
#'     A numeric vector specifying the \strong{anagenesis rates}
#'     (\eqn{\lambda_a}) for each trait state. Each element in the vector
#'     corresponds to the anagenesis rate for a different state.
#'   }
#'
#'   \item{\code{idparslist[[5]]}}{
#'     A transition matrix describing the transition rates between all
#'     trait states. The matrix is square, and its entries represent
#'     transition rates between states:
#'
#'     \describe{
#'       \item{\code{idparslist[[5]][i, j]}}{
#'         The transition rate from state \eqn{i} to state \eqn{j},
#'         where \eqn{i} and \eqn{j} range from 1 to the total number
#'         of states.
#'       }
#'
#'       \item{Diagonal elements}{
#'         \code{idparslist[[5]][i, i]} represents the self-transition
#'         rate and is equal to 0.
#'       }
#'     }
#'   }
#'
#'   \item{\code{idparslist[[6]]}}{
#'     A fixed value for the parameter \eqn{p}, which specifies the
#'     probability that a trait transition results in a new species.
#'     If \eqn{p = 1}, every transition results in a new species.
#'     If \eqn{p = 0}, the transition does not lead to the creation
#'     of a new species.
#'   }
#' }

#' @param idparsopt A vector of integers specifying which parameters are to be
#' optimized.
#' @param idparsfix A vector of integers specifying which parameters are to be
#' kept fixed.
#' @param parsfix A vector of values corresponding to the fixed parameters
#' specified in `idparsfix`.
#' @param verbose sets verbose output; default is TRUE when optimmethod is
#' "simplex". If optimmethod is set to "simplex", then even if set to FALSE,
#' optimizer output will be shown.
#'
#' @details
#' This function runs the maximum likelihood estimation for a dataset based on
#' the trait-dependent model. The optimization process uses the
#' `TRAISIE::calc_ml` function to adjust the parameters for the best fit to the
#' data. The function assumes that the dataset is prepared correctly, and the
#' parameter IDs are provided for both optimization and fixing.
#'
#' @returns
#' A list containing the following elements:
#' \describe{
#'   \item{\code{MLpars}}{The optimized parameter values.}
#'   \item{\code{ML}}{The calculated log-likelihood value for the model.}
#'   \item{\code{conv}}{The convergence status of the optimization process.}
#' }
#'
#' @examples
#' \dontrun{
#' ############## Parameter Estimation
#' ### Binary State Model without Hidden States
#'
#' # Define the parameter list for the model
#' idparslist <- list()
#' idparslist[[1]] <- c(1, 1)  # lambda_c
#' idparslist[[2]] <- c(2, 2)  # mu
#' idparslist[[3]] <- c(3, 3)  # gamma
#' idparslist[[4]] <- c(4, 4)  # lambda_a
#'
#' # Transition matrix for the model (since there are two trait states)
#' idparslist[[5]] <- matrix(0, 2, 2)
#' colnames(idparslist[[5]]) <- c("0", "1")
#' rownames(idparslist[[5]]) <- colnames(idparslist[[5]])
#'
#' # Hidden state transitions (irrelevant here but kept for consistency)
#' idparslist[[5]][1, 2] <- 5  # 0 -> 1
#' idparslist[[5]][2, 1] <- 6  # 1 -> 0
#' idparslist[[5]][1, 1] <- 7  # 0 -> 0
#' idparslist[[5]][2, 2] <- 7  # 1 -> 1
#'
#' # Set p (parameter for anagenesis transition) not to be estimated
#' idparslist[[6]] <- 8
#'
#' # Parameters to optimize
#' idparsopt <- 1:6
#' idparsopt <- idparsopt[idparsopt > 0]  # We will optimize parameters 1 to 6
#'
#' # Example of starting values for the optimization
#' initvals2 <- c(1.07, 1.0102, 0.0035, 0.174, 1, 1)
#'
#' # Preparing the parameter optimization values
#' initparsopt <- initvals2
#' idparsfix = c(0, 7, 8) # We will not optimize parameters 1 to 6
#' parsfix <- c(0, 0, 0)  # Fixed parameter, qs and p not to be estimated
#' trparsopt <- initparsopt / (1 + initparsopt)
#' trparsopt[which(initparsopt == Inf)] <- 1
#' trparsfix <- parsfix / (1 + parsfix)
#' trparsfix[which(parsfix == Inf)] <- 1
#'
#' # Run the ML estimation with the provided dataset
#' ml_estimates <- TRAISIE::calc_ml(datalist,
#'                                 num_observed_states = 2,
#'                                 num_hidden_states = 1,
#'                                 idparslist = idparslist,
#'                                 idparsopt = idparsopt,
#'                                 initparsopt = initvals2,
#'                                 idparsfix = idparsfix,
#'                                 parsfix = parsfix,
#'                                 atol = 1e-15,
#'                                 rtol = 1e-15,
#'                                 num_threads = 8,
#'                                 verbose = TRUE)
#'
#' # View the results
#' print(ml_estimates)
#' }
#'
#'
#' @export
TRAISIE_calc_ml <- function(datalist,
                            num_observed_states,
                            num_hidden_states,
                            idparslist,
                            idparsopt,
                            initparsopt,
                            idparsfix,
                            parsfix,
                            cond = "proper_cond",
                            tol = c(1e-04, 1e-05, 1e-07),
                            maxiter = 2000 * round((1.25) ^ length(idparsopt)),
                            optimmethod = "simplex",
                            methode = "odeint::runge_kutta_cash_karp54",
                            num_cycles = 1,
                            verbose = FALSE,
                            num_threads = 1,
                            atol = 1e-15,
                            rtol = 1e-15
) {
  if (identical(as.numeric(sort(c(idparsopt, idparsfix))),
                as.numeric(sort(unique(unlist(idparslist))))) == FALSE) {
    stop("All elements in idparslist must be included in either
             idparsopt or idparsfix ")
  }

  trparsopt <- initparsopt / (1 + initparsopt)
  trparsopt[which(initparsopt == Inf)] <- 1
  trparsfix <- parsfix / (1 + parsfix)
  trparsfix[which(parsfix == Inf)] <- 1

  optimpars <- c(tol, maxiter, verbose)

  initloglik <- TRAISIE_loglik_choosepar(trparsopt = trparsopt,
                                         trparsfix = trparsfix,
                                         idparsopt = idparsopt,
                                         idparsfix = idparsfix,
                                         idparslist = idparslist,
                                         datalist = datalist,
                                         num_observed_states = num_observed_states,
                                         num_hidden_states = num_hidden_states,
                                         cond = cond,
                                         atol = atol,
                                         rtol = rtol,
                                         methode = methode,
                                         verbose = verbose,
                                         num_threads = num_threads)

  print_init_ll(initloglik = initloglik, verbose = verbose)

  if (initloglik == -Inf) {
    stop("The initial parameter values have a likelihood that is
             equal to 0 or below machine precision.
             Try again with different initial values.")
  } else {
    out <- DDD::optimizer(optimmethod = optimmethod,
                          optimpars = optimpars,
                          fun = TRAISIE_loglik_choosepar,
                          trparsopt = trparsopt,
                          num_cycles = num_cycles,
                          idparsopt = idparsopt,
                          trparsfix = trparsfix,
                          idparsfix = idparsfix,
                          idparslist = idparslist,
                          datalist = datalist,
                          num_observed_states = num_observed_states,
                          num_hidden_states = num_hidden_states,
                          cond = cond,
                          atol = atol,
                          rtol = rtol,
                          methode = methode,
                          verbose = verbose,
                          num_threads = num_threads)
    if (out$conv != 0) {
      stop("Optimization has not converged.
                 Try again with different initial values.")
    } else {
      ml_pars1 <- transform_parameters(as.numeric(unlist(out$par)),
                                       trparsfix,
                                       idparsopt,
                                       idparsfix,
                                       idparslist)
      names(ml_pars1) <- c(
        "ML_lambda_c",
        "ML_mu",
        "ML_gamma",
        "ML_lambda_a",
        "ML_q",
        "p"
      )
      out2 <- list(MLpars = ml_pars1,
                   ML = as.numeric(unlist(out$fvalues)),
                   conv = out$conv)
    }
  }
  return(out2)
}

#' @keywords internal
TRAISIE_loglik_choosepar <- function(trparsopt,
                                     trparsfix,
                                     idparsopt,
                                     idparsfix,
                                     idparslist,
                                     datalist,
                                     num_observed_states,
                                     num_hidden_states,
                                     cond = cond,
                                     atol,
                                     rtol,
                                     methode,
                                     verbose,
                                     num_threads) {
  alltrpars <- c(trparsopt, trparsfix)

  loglik <- NA

  if (max(alltrpars) > 1 || min(alltrpars) < 0) {
    loglik <- -Inf
  } else {
    pars1 <- transform_parameters(trparsopt, trparsfix,
                                  idparsopt, idparsfix,
                                  idparslist)

    loglik <- TRAISIE_loglik_CS(parameter = pars1,
                                datalist = datalist,
                                methode = methode,
                                atol = atol,
                                rtol = rtol,
                                num_observed_states =
                                  num_observed_states,
                                num_hidden_states = num_hidden_states,
                                cond = cond,
                                verbose = verbose,
                                num_threads = num_threads)

    if (is.nan(loglik) || is.na(loglik)) {
      warning("There are parameter values used which cause
                numerical problems.")
      loglik <- -Inf
    }
  }
  return(loglik)
}
