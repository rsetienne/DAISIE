#' The expected number of endemics and non-endemics under the DAISIE model with
#' no diversity-dependence
#'
#' This function calculates the expected number of endemics, non-endemics and
#' the sum of these for a given set of parameter values, a given mainland
#' species pool size and a given time, assuming no diversity-dependence
#'
#' @inheritParams default_params_doc
#'
#' @return \item{out}{The output is a list with three elements: \cr \cr
#' \code{ExpE} The number of endemic species \cr
#' \code{ExpI} The number of non-endemic species \cr
#' \code{ExpN} The sum of the number of endemics and non-endemics }
#' @author Rampal S. Etienne
#' @references Valente, L.M., A.B. Phillimore and R.S. Etienne (2015).
#' Equilibrium and non-equilibrium dynamics simultaneously operate in the
#' Galapagos islands. Ecology Letters 18: 844-852.
#' @keywords models
#' @examples
#'
#' ### Compute the expected values at t = 4, for a mainland pool size of 1000 potential
#' # colonists and a vector of 5 parameters (cladogenesis, extinction, clade-level carrying
#' # capacity, immigration, anagenesis)
#'
#' DAISIE_ExpEIN(
#'    tvec = c(2,4),
#'    pars = c(0.5,0.1,Inf,0.01,0.4),
#'    M = 1000
#'    )
#'
#' @export DAISIE_ExpEIN
DAISIE_ExpEIN <- function(tvec, pars, M, initEI = c(0, 0)) {
   if(pars[3] != Inf) warning('K must be infinite; proceeding with infinite K')
   pars1 <- pars
   lac <- pars1[1]
   mu <- pars1[2]
   ga <- pars1[4]
   laa <- pars1[5]
   if (!is.na(pars1[11])) {
       M2 <- M - DDD::roundn(pars1[11] * M)
   } else {
       M2 <- M
   }
   A <- mu - lac
   B <- lac + mu + ga + laa
   C <- laa + 2 * lac + ga
   DD <- laa + 2 * lac
   E0 <- initEI[1]
   I0 <- initEI[2]
   if (tvec[1] == Inf) {
      Imm <- ga * M2 / B
      End <- DD / A * Imm
   } else {
      #Imm = M2 * ga / B * (1 - exp(-B * t))
      #End = M2 * ga * (laa + 2 * lac) * (1/(A * B) -
      #exp(-A*tvec) / (A * C) + exp(-B*tvec)/(B * C))
      Imm <- M2 * ga / B - (M2 * ga / B - I0) * exp(-B * tvec)
      End <- DD / C * (M2 * ga / A - M2 * ga / B +
                        (C / DD * E0 - M2 * ga / A + I0) *
                        exp(-A * tvec) + (M2 * ga / B - I0) * exp(-B * tvec))
   }
   All <- End + Imm
   expEIN <- list(End, Imm, All)
   names(expEIN) <- c("ExpE", "ExpI", "ExpN")
   return(expEIN)
}

#' The expected number of endemics and non-endemics under the DAISIE model
#'
#' This function calculates the expected number of endemics, non-endemics and
#' the sum of these for a given set of parameter values, a given mainland
#' species pool size and a given time, where there can be diversity-dependence
#'
#' @inheritParams default_params_doc
#' @param pars2 list of settings
#' \itemize{
#' \item{\code{res}: the number of equations}
#' \item{\code{ddep}: the model of diversity-dependence}
#' \item{\code{methode}: the method used to integrate the ODE system}
#' \item{\code{reltolint}: the relative tolerance in integration}
#' \item{\code{abstolint}: the absolute tolerance in integration}
#' }
#'
#' @return \item{out}{The output is a list with three elements: \cr \cr
#' \code{ExpE} The number of endemic species \cr \code{ExpI} The number of
#' non-endemic species \cr \code{ExpN} The sum of the number of endemics and
#' non-endemics }
#' @author Rampal S. Etienne
#' @export DAISIE_ExpEIN2
DAISIE_ExpEIN2 <- function(tvec,
                           pars,
                           M,
                           initEI = c(0, 0),
                           res = 1000,
                           ddmodel = 11,
                           methode = 'odeint::runge_kutta_fehlberg78',
                           reltolint = 1E-16,
                           abstolint = 1E-16) {
  pars1 <- pars
  lac <- pars1[1]
  mu <- pars1[2]
  K <- pars1[3]
  ga <- pars1[4]
  laa <- pars1[5]
  if (!is.na(pars1[11])) {
    M2 <- M - DDD::roundn(pars1[11] * M)
  } else {
    M2 <- M
  }
  Kprime <- lac/(lac - mu) * K
  res <- ceiling(min(Kprime,res))
  initprobs <- rep(0,2 * res + 1)
  initprobs[initEI[1] + 1] <- 1
  initprobs[initEI[2]] <- 1
  tvec <- sort(abs(tvec))
  if(tvec[1] != 0) tvec <- c(0,tvec)
  probs <- DAISIE_integrate(initprobs = initprobs,
                            tvec = tvec,
                            rhs_func = DAISIE_loglik_rhs,
                            pars = c(pars1,0,ddmodel),
                            rtol = reltolint,
                            atol = abstolint,
                            method = methode)
  dp <- dim(probs)
  if(is.null(dp)) {
    probs <- matrix(probs, nrow = 1, byrow = TRUE)
    dp <- dim(probs)
  }
  probs1 <- matrix(probs[,1:res], nrow = dp[1], byrow = FALSE)
  probs2 <- matrix(probs[,(res + 1):(2 * res)], nrow = dp[1], byrow = FALSE)
  nil2resmin1 <- matrix(rep(0:(res - 1), each = dp[1]), nrow = dp[1], byrow = FALSE)
  End <- M2 * rowSums((probs1 + probs2) * nil2resmin1)
  Imm <- M2 * rowSums(probs2)
  All <- End + Imm
  if(any(All > 0.5 * res) & res < Kprime) warning('Result is probably not accurate.
                              Increase the number of equations (pars2$res')
  expEIN <- list(End, Imm, All)
  names(expEIN) <- c("ExpE", "ExpI", "ExpN")
  return(expEIN)
}
