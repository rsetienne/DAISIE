#' Maximum likelihood estimation for the archipelago model (K = Inf)
#'
#' @description Maximises \code{\link{DAISIE_A_loglik}()} over a vector of
#'   free parameters. The mapping from the free parameters to the full
#'   parameter object is given by the user as a function, so that any
#'   constraint can be expressed: equal rates on all islands, rates that
#'   depend on island area, a distance-based connectivity matrix, fixed
#'   parameters, and so on.
#'
#' @param datalist Archipelago datalist, see \code{\link{DAISIE_A_loglik}()}.
#' @param parsfun Function that takes a vector of non-negative free
#'   parameters and returns a \code{DAISIE_A_pars} object (see the example).
#' @param initparsopt Initial values of the free parameters (all
#'   non-negative).
#' @param parnames Optional names of the free parameters.
#' @param sampling,cond Passed to \code{\link{DAISIE_A_loglik}()}.
#' @param optimmethod Optimisation method passed to \code{DDD::optimizer()}
#'   (\code{"simplex"} or \code{"subplex"}).
#' @param optimpars Optimisation settings passed to \code{DDD::optimizer()}.
#' @param num_cycles Number of optimisation cycles.
#' @param verbose Logical; print progress.
#' @param ... Further arguments passed to \code{\link{DAISIE_A_loglik}()}
#'   (e.g. ODE tolerances).
#'
#' @return A data frame with the maximum likelihood estimates of the free
#'   parameters, the maximum log-likelihood (\code{loglik}), the number of
#'   free parameters (\code{df}) and the convergence code (\code{conv}).
#' @family DAISIE_A
#' @export
#' @examples
#' \donttest{
#' pars <- DAISIE_A_create_pars(n_islands = 2, lac = 0.5, mu = 0.3,
#'   gam = 0.05, laa = 0.5, lav = 0.3, disp = 0.2)
#' set.seed(1)
#' sim <- DAISIE_A_sim(time = 2, M = 50, pars = pars, replicates = 1,
#'   verbose = FALSE)
#' datalist <- DAISIE_A_sim_to_datalist(sim[[1]])
#' # estimate lac, mu, gam and disp; laa and lav fixed at their true values
#' parsfun <- function(x) {
#'   DAISIE_A_create_pars(n_islands = 2, lac = x[1], mu = x[2], gam = x[3],
#'     laa = 0.5, lav = 0.3, disp = x[4])
#' }
#' DAISIE_A_ML(datalist, parsfun, initparsopt = c(0.4, 0.2, 0.04, 0.3),
#'   parnames = c("lac", "mu", "gam", "disp"))
#' }
DAISIE_A_ML <- function(datalist,
                        parsfun,
                        initparsopt,
                        parnames = NULL,
                        sampling = "n",
                        cond = 0,
                        optimmethod = "simplex",
                        optimpars = c(1e-4, 1e-4, 1e-6, 1000),
                        num_cycles = 1,
                        verbose = FALSE,
                        ...) {
  if (any(initparsopt < 0)) stop("Initial parameter values must be >= 0.")
  if (is.null(parnames)) parnames <- paste0("par", seq_along(initparsopt))
  fun <- function(trparsopt) {
    if (any(trparsopt < 0) || any(trparsopt >= 1)) return(-Inf)
    x <- trparsopt / (1 - trparsopt)
    ll <- tryCatch({
      pars <- parsfun(x)
      DAISIE_A_loglik(pars = pars, datalist = datalist, sampling = sampling,
                      cond = cond, ...)
    }, error = function(e) -Inf)
    if (!is.finite(ll)) ll <- -Inf
    if (verbose) {
      message(paste(signif(x, 5), collapse = " "), " : ", signif(ll, 8))
    }
    ll
  }
  trparsopt <- initparsopt / (1 + initparsopt)
  initloglik <- fun(trparsopt)
  if (!is.finite(initloglik)) {
    stop("The initial parameter values give a log-likelihood of -Inf.")
  }
  out <- DDD::optimizer(optimmethod = optimmethod,
                        optimpars = c(optimpars[1:4], as.numeric(verbose)),
                        num_cycles = num_cycles,
                        fun = fun,
                        trparsopt = trparsopt)
  est <- out$par / (1 - out$par)
  res <- as.data.frame(t(est))
  names(res) <- parnames
  res$loglik <- out$fvalues
  res$df <- length(initparsopt)
  res$conv <- out$conv
  return(res)
}
