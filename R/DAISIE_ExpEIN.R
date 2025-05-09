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
#' @param initEI matrix where each row represents the initial number of endemic
#' and non-endemic species per colonizing lineage.
#' \itemize{
#' \item{\code{res}: the number of equations}
#' \item{\code{ddep}: the model of diversity-dependence}
#' \item{\code{methode}: the method used to integrate the ODE system}
#' \item{\code{reltolint}: the relative tolerance in integration}
#' \item{\code{abstolint}: the absolute tolerance in integration}
#' }
#'
#' @return \item{tot_expEIN}{The output is a list with three elements: \cr \cr
#' \code{ExpE} The number of endemic species at the times in tvec\cr
#' \code{ExpI} The number of non-endemic species at the times in tvec\cr
#' \code{ExpN} The sum of the number of endemics and non-endemics at the times
#' in tvec}
#' @author Rampal S. Etienne
#' @examples DAISIE_ExpEIN2(tvec = c(0.000001,0.5,0.75,1),
#'                          pars = c(0.3,0.1,10,1,0.1),
#'                          M = 1000,
#'                          initEI = rbind(c(1,0),c(2,0),c(0,1)))
#' @export DAISIE_ExpEIN2
DAISIE_ExpEIN2 <- function(tvec,
                           pars,
                           M,
                           initEI = NULL,
                           res = 1000,
                           ddmodel = 11,
                           methode = 'ode45',
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
  res <- ceiling(min(K,res))
  tvec <- sort(abs(tvec))
  if(tvec[1] != 0) tvec <- c(0,tvec)
  if(is.null(initEI) | all(initEI == c(0,0))) {
    initEI <- t(c(0,0))
  } else {
    initEI <- rbind(c(0,0),initEI)
  }
  num_of_lin <- nrow(initEI)
  if(M < num_of_lin - 1) warning('M should be a positive integer.')
  expEIN <- list()
  tot_expEIN <- list()
  for(i in 1:num_of_lin) {
    initprobs <- rep(0,2 * res + 1)
    if(initEI[i,2] > 0) {
      initprobs[initEI[i,2] + initEI[i,1] + res] <- 1
    } else {
      initprobs[initEI[i,1] + 1] <- 1
    }
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
    fac <- ifelse(i == 1,M2 - num_of_lin + 1, 1)
    End <- fac * rowSums((probs1 + probs2) * nil2resmin1)
    Imm <- fac * rowSums(probs2)
    All <- End + Imm
    if(any(All > 0.5 * res) & res < K) warning('Result is probably not accurate.
                              Increase the number of equations (pars2$res')
    expEIN[[i]] <- list(End, Imm, All)
    if(i == 1) {
      tot_expEIN <- expEIN[[1]]
      names(tot_expEIN) <- c("ExpE", "ExpI", "ExpN")
    } else {
      tot_expEIN <- list(ExpE = tot_expEIN$ExpE + End,
                         ExpI = tot_expEIN$ExpI + Imm,
                         ExpN = tot_expEIN$ExpN + All)
    }
    names(expEIN[[i]]) <- c("ExpE", "ExpI", "ExpN")
  }
  return(tot_expEIN)
}
#' The probability distribution of the number of endemics and non-endemics under
#' the DAISIE model
#'
#' This function calculates the probability distribution of the number of
#' endemics, non-endemics and the sum of these for a given set of parameter
#' values, a given mainland species pool size and a given time, where there can
#' be diversity-dependence
#'
#' @inheritParams default_params_doc
#' @param pars2 list of settings
#' @param initEI matrix where each row represents the initial number of endemic
#' and non-endemic species per colonizing lineage.
#' \itemize{
#' \item{\code{res}: the number of equations}
#' \item{\code{ddep}: the model of diversity-dependence}
#' \item{\code{methode}: the method used to integrate the ODE system}
#' \item{\code{reltolint}: the relative tolerance in integration}
#' \item{\code{abstolint}: the absolute tolerance in integration}
#' }
#'
#' @return \item{probsEIN}{The output is a list with three elements: \cr \cr
#' \code{probsE} The number of endemic species at the times in tvec\cr
#' \code{probsI} The number of non-endemic species at the times in tvec \cr
#' \code{probsN} The sum of the number of endemics and non-endemics at the times
#' in tvec}
#' @author Rampal S. Etienne
#' @examples DAISIE_margprobdist2(tvec = c(0.000001,0.5,0.75,1),
#'                                pars = c(0.3,0.1,10,1,0.1),
#'                                M = 1000,
#'                                initEI = rbind(c(1,0),c(2,0),c(0,1)))
#' @export DAISIE_margprobdist2
DAISIE_margprobdist2 <- function(tvec,
                                 pars,
                                 M,
                                 initEI = NULL,
                                 res = 1000,
                                 ddmodel = 11,
                                 methode = 'ode45',
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
  res <- ceiling(min(K,res))
  tvec <- sort(abs(tvec))
  if(tvec[1] != 0) tvec <- c(0,tvec)
  if(is.null(initEI) | all(initEI == c(0,0))) {
    initEI <- t(c(0,0))
  } else {
    initEI <- rbind(c(0,0),initEI)
  }
  num_of_lin <- nrow(initEI)
  if(M < num_of_lin - 1) warning('M should be a positive integer.')
  for(i in 1:num_of_lin) {
    initprobs <- rep(0,2 * res + 1)
    if(initEI[i,2] > 0) {
      initprobs[initEI[i,2] + initEI[i,1] + res] <- 1
    } else {
      initprobs[initEI[i,1] + 1] <- 1
    }
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
    probs3 <- probs1 + probs2
    probs4 <- cbind(1 - rowSums(probs2),rowSums(probs2))
    probs5 <- cbind(probs1[,1],probs1[,2:res] + probs2[,1:(res - 1)])
    if(i == 1) {
      probsE <- convolve_rows_fft(mat = probs3, n = M2 - num_of_lin + 1)
      probsI <- convolve_rows_fft(mat = probs4, n = M2 - num_of_lin + 1)
      probsN <- convolve_rows_fft(mat = probs5, n = M2 - num_of_lin + 1)
    } else {
      probsE <- convolve_rows_fft2(mat1 = probsE, mat2 = probs3)
      probsI <- convolve_rows_fft2(mat1 = probsI, mat2 = probs4)
      probsN <- convolve_rows_fft2(mat1 = probsN, mat2 = probs5)
    }
    probsEIN <- list(probsE = probsE, probsI = probsI, probsN = probsN)
    nil2E <- matrix(rep(0:(ncol(probsE) - 1), each = dp[1]), nrow = dp[1], byrow = FALSE)
    nil2I <- matrix(rep(0:(ncol(probsI) - 1), each = dp[1]), nrow = dp[1], byrow = FALSE)
    nil2N <- matrix(rep(0:(ncol(probsN) - 1), each = dp[1]), nrow = dp[1], byrow = FALSE)
    End <- rowSums(probsE * nil2E)
    Imm <- rowSums(probsI * nil2I)
    All <- rowSums(probsN * nil2N)
    if(any(All > 0.5 * res) & res < K) warning('Result is probably not accurate.
                              Increase the number of equations (pars2$res')
    names(probsEIN) <- c("probsE", "probsI", "probsN")
  }
  return(probsEIN)
}

convolve_rows_fft <- function(mat, n) {
  num_rows <- nrow(mat)
  input_len <- ncol(mat)
  output_len <- (input_len - 1) * n + 1
  results <- matrix(0, nrow = num_rows, ncol = output_len)
  for (i in 1:num_rows) {
    v <- mat[i, ]
    v_padded <- c(v, rep(0, output_len - length(v)))
    V <- stats::fft(v_padded)
    V_power <- V^n
    result <- pmax(0,Re(stats::fft(V_power, inverse = TRUE)) / output_len)
    results[i, ] <- result
  }
  return(results)
}

convolve_rows_fft2 <- function(mat1, mat2) {
  output_len <- ncol(mat1) + ncol(mat2) - 1
  fft_len <- 2^ceiling(log2(output_len))
  mat3 <- matrix(0, nrow = nrow(mat1), ncol = output_len)
  for(i in 1:nrow(mat1)) {
    mat1_padded <- c(mat1[i,], rep(0, fft_len - length(mat1[i,])))
    mat2_padded <- c(mat2[i,], rep(0, fft_len - length(mat2[i,])))
    tr_mat1 <- stats::fft(mat1_padded)
    tr_mat2 <- stats::fft(mat2_padded)
    tr_mat3 <- tr_mat1 * tr_mat2
    mat3[i,] <- pmax(0,Re(stats::fft(tr_mat3, inverse = TRUE)) / fft_len)[1:output_len]
  }
  return(mat3)
}
