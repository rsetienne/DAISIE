#' IW concurrency control
#'
#' Sets or retrieves the number of threads used by the odeint solver.
#'
#' @param num_threads \code{num_threads < 0 or omitted}: retrieves number of threads. \cr
#' \code{num_threads = 0}: sets the number of threads to the number of available cores. \cr
#' \code{num_threads = 1}: single-threaded execution. \cr
#' \code{num_threads > 1}: sets the number of threads to \code{num_threads}.
#' @return number of threads
#' @note The maximum number of threads is limited to the value of the C++
#' standard library function \code{std::thread::hardware_concurrency()}.
#' This is also the default number of threads upon library load.
#' Multithreading incurs some overhead. Therefore, single-threaded execution
#' might be faster for small systems.
#'
#' @export DAISIE_IW_num_threads
DAISIE_IW_num_threads <- function(num_threads) {
  if (missing(num_threads)) {
    # retrieve only
    return(.Call("daisie_odeint_iw_num_threads", -1))
  }
  return(.Call("daisie_odeint_iw_num_threads", num_threads))
}


dec2bin <- function(y, ly) {
  stopifnot(length(y) == 1, mode(y) == "numeric")
  q1 <- (y / 2) %/% 1
  r <- y - q1 * 2
  res <- c(r)
  while (q1 >= 1) {
    q2 <- (q1 / 2) %/% 1
    r <- q1 - q2 * 2
    q1 <- q2
    res <- c(r, res)
  }
  res <- c(rep(0, ly - length(res)), res)
  return(res)
}

dec2binmat <- function(y) {
  numrows <- 2 ^ y
  res <- matrix(0, numrows, y)
  for (i in 0:(numrows - 1)) {
    res[i + 1, ] <- dec2bin(i, y)
  }
  return(res)
}

bin2dec <- function(y) {
  res <- y %*% 2 ^ ((length(y) - 1):0)
  return(as.numeric(res))
}

kimat <- function(dec2binmatk) {
  ki <- matrix(0, dim(dec2binmatk)[1], dim(dec2binmatk)[1])
  for (i in 2:dim(dec2binmatk)[1]) {
    locationones <- which(dec2binmatk[i, ] == 1)
    for (j in 1:length(locationones)) {
      dec2binmatki <- dec2binmatk[i, ]
      dec2binmatki[locationones[j]] <- 0
      j2 <- 1 + bin2dec(dec2binmatki)
      ki[i,j2] <- 1
    }
  }
  return(ki)
}

create_l0ki <- function(dec2binmatk, lxm, lxe, sysdim = dim(dec2binmatk)[1]) {
  posc <- Matrix::rowSums(dec2binmatk)
  negc <- log2(sysdim) - posc
  l0 <- rep(negc, each = lxm * lxe)
  dim(l0) <- c(lxm, lxe, sysdim)
  ki <- kimat(dec2binmatk)
  res <- list(l0 = l0, ki = ki)
  return(res)
}

nndivdep <- function(lxm, lxe, sysdim, Kprime, k, M, l0) {
  nnm <- c(0, 0:(lxm + 1))
  nne <- c(0, 0, 0:(lxe + 1))
  lnnm <- length(nnm)
  lnne <- length(nne)
  nil2lxm <- 2:(lxm + 1)
  nil2lxe <- 3:(lxe + 2)
  nn <- rowSums(expand.grid(n1 = nnm, n2 = nne))
  dim(nn) <- c(lnnm, lnne)
  nn <- replicate(sysdim, nn)
  nilm <- rep(1, lxm)
  nile <- rep(1, lxe)
  allc <- 1:sysdim
  divdepfac <- pmax(array(0, dim = c(lxm + 3, lxe + 4, sysdim)),
                    1 - (nn + k) / Kprime)
  divdepfacmin1 <- pmax(array(0, dim = c(lxm + 3, lxe + 4, sysdim)),
                        1 - (nn + k - 1) / Kprime)
  divdepfacplus1 <- pmax(array(0, dim = c(lxm + 3, lxe + 4, sysdim)),
                         1 - (nn + k + 1) / Kprime)
  divdepfac <- divdepfac[nil2lxm, nil2lxe, allc]
  divdepfacmin1 <- divdepfacmin1[nil2lxm, nil2lxe, allc]
  divdepfacplus1 <- divdepfacplus1[nil2lxm, nil2lxe, allc]
  mfac <- (nn[nil2lxm,nile,allc] + 1)/(M - l0)
  oneminmfac <- (M - nn[nil2lxm,nile,allc] - l0)/(M - l0)
  res <- list(nn = nn,
              divdepfac = divdepfac,
              divdepfacmin1 = divdepfacmin1,
              divdepfacplus1 = divdepfacplus1,
              mfac = mfac,
              oneminmfac = oneminmfac)
  return(res)
}

selectrows <- function(sysdim, order) {
  mat <- NULL
  for (i in seq(1, sysdim / order, by = 2)) {
    group1 <- (i - 1) * order + (1:order)
    group2 <- i * order + (1:order)
    mat <- rbind(mat, cbind(group1, group2))
  }
  return(mat)
}

DAISIE_IW_pars <- function(parslist) {
  lac <- parslist$pars[1]
  mu <- parslist$pars[2]
  Kprime <- parslist$pars[3]
  gam <- parslist$pars[4]
  laa <- parslist$pars[5]
  M <- parslist$pars[6]
  k <- parslist$k
  ddep <- parslist$ddep
  lxm <- parslist$dime$lxm
  lxe <- parslist$dime$lxe
  sysdim <- parslist$dime$sysdim
  l0 <- parslist$l0ki$l0
  ki <- parslist$l0ki$ki
  nn <- parslist$nndd$nn
  divdepfac <- parslist$nndd$divdepfac
  divdepfacmin1 <- parslist$nndd$divdepfacmin1
  nil2lxm <- 2:(lxm + 1)
  nil2lxe <- 3:(lxe + 2)
  allc <- 1:sysdim
  nilm <- rep(1, lxm)
  nile <- rep(1, lxe)
  cp <- list(
    laa = laa,
    ki = ki,
    lxm = lxm,
    lxe = lxe,
    sysdim = sysdim,
    c1 = gam * divdepfacmin1 * pmax(0,M - l0 - nn[nil2lxm - 1, nile, allc]),
    c2 = mu * nn[nil2lxm + 1, nile, allc],
    c3 = mu * nn[nilm, nil2lxe + 1, allc],
    c4 = laa * nn[nil2lxm + 1, nile, allc],
    c5 = lac * divdepfacmin1 * nn[nil2lxm + 1, nile, allc],
    c6 = lac * divdepfacmin1 * (2 * (k - l0) + nn[nilm, nil2lxe - 1, allc]),
    c7 = gam * divdepfac * pmax(0,M - nn[nil2lxm, nile, allc]) +
         mu * (k + nn[nil2lxm, nil2lxe, allc]) +
         laa * (l0 + nn[nil2lxm, nile, allc]) +
         lac * divdepfac * (k + nn[nil2lxm, nil2lxe, allc]),
    c8 = 2 * lac * divdepfacmin1
  )
  return(cp)
}

DAISIE_loglik_rhs_IW <- function(t,x,cp)
{
  lxm <- cp$lxm
  lxe <- cp$lxe
  sysdim <- cp$sysdim
  dim(x) <- c(lxm ,lxe, sysdim)
  xx <- array(0,dim = c(lxm + 2, lxe + 3, sysdim))
  nil2lxm <- 2:(lxm + 1)
  nil2lxe <- 3:(lxe + 2)
  allc <- 1:sysdim
  xx[nil2lxm, nil2lxe, allc] <- x
  dx <-
    cp$c1 * xx[nil2lxm - 1, nil2lxe, allc] +
    cp$c2 * xx[nil2lxm + 1, nil2lxe, allc] +
    cp$c3 * xx[nil2lxm, nil2lxe + 1, allc]  +
    cp$c4 * xx[nil2lxm + 1, nil2lxe - 1, allc] +
    cp$c5 * xx[nil2lxm + 1, nil2lxe - 2, allc] +
    cp$c6 * xx[nil2lxm, nil2lxe - 1, allc] -
    cp$c7 * xx[nil2lxm, nil2lxe, allc]
  if (sysdim > 1) {
    dx <- dx +
      cp$laa * tensor::tensor(xx[nil2lxm, nil2lxe, allc], cp$ki, 3, 2) +
      cp$c8 * tensor::tensor(xx[nil2lxm, nil2lxe - 1, allc], cp$ki, 3, 2)
  }
  dim(dx) <- c(lxm * lxe * sysdim, 1)
  return(list(dx))
}

#' Computes the loglikelihood of the DAISIE model with island-wide
#' diversity-dependence given data and a set of model parameters
#'
#' Computes the loglikelihood of the DAISIE model given colonization and
#' branching times for lineages on an island, and a set of model parameters for
#' the DAISIE model with island-wide diversity-dependence
#'
#' The output is a loglikelihood value
#'
#' @param pars1 Contains the model parameters: \cr \cr
#' \code{pars1[1]} corresponds to lambda^c (cladogenesis rate) \cr
#' \code{pars1[2]} corresponds to mu (extinction rate) \cr
#' \code{pars1[3]} corresponds to K (clade-level carrying capacity) \cr
#' \code{pars1[4]} corresponds to gamma (immigration rate) \cr
#' \code{pars1[5]} corresponds to lambda^a (anagenesis rate) \cr
#' \code{pars1[6]} is optional; it may contain M, the total number of species
#' on the mainland \cr \cr
#' @param pars2 Contains the model settings \cr \cr
#' \code{pars2[1]} corresponds to lx = length of ODE variable x \cr
#' \code{pars2[2]} corresponds to ddmodel = diversity-dependent model, model of diversity-dependence, which can be one
#' of\cr \cr
#' ddmodel = 0 : no diversity dependence \cr
#' ddmodel = 1 : linear dependence in speciation rate \cr
#' ddmodel = 11: linear dependence in speciation rate and in immigration rate \cr
#' ddmodel = 2 : exponential dependence in speciation rate\cr
#' ddmodel = 21: exponential dependence in speciation rate and in immigration rate\cr
#' Only ddmodel = 11 is currently implemented \cr \cr
#' \code{pars2[3]} corresponds to cond = setting of conditioning\cr \cr
#' cond = 0 : conditioning on island age \cr
#' cond = 1 : conditioning on island age and non-extinction of the island biota \cr \cr
#' \code{pars2[4]} Specifies whether intermediate output should be provided,
#' because computation may take long. Default is 0, no output. A value of 1
#' means the parameters and loglikelihood are printed. A value of 2 means also
#' intermediate progress during loglikelihood computation is shown.
#' @param datalist Data object containing information on colonisation and
#' branching times. This object can be generated using the DAISIE_dataprep
#' function, which converts a user-specified data table into a data object, but
#' the object can of course also be entered directly. It is an R list object
#' with the following elements.\cr
#' The first element of the list has two or
#' three components: \cr \cr
#' \code{$island_age} - the island age \cr
#' Then, depending on whether a distinction between types is made, we have:\cr
#' \code{$not_present} - the number of mainland lineages that are not present
#' on the island \cr The remaining elements of the list each contains
#' information on a single colonist lineage on the island and has 5
#' components:\cr \cr
#' \code{$colonist_name} - the name of the species or clade
#' that colonized the island \cr
#' \code{$branching_times} - island age and stem
#' age of the population/species in the case of Non-endemic, Non-endemic_MaxAge
#' and Endemic anagenetic species. For cladogenetic species these should be
#' island age and branching times of the radiation including the stem age of
#' the radiation.\cr
#' \code{$stac} - the status of the colonist \cr \cr
#' * Non_endemic_MaxAge: 1 \cr
#' * Endemic: 2 \cr
#' * Endemic&Non_Endemic: 3 \cr
#' * Non_endemic: 4 \cr
#' * Endemic_MaxAge: 5 \cr \cr
#' \code{$missing_species} -
#' number of island species that were not sampled for particular clade (only
#' applicable for endemic clades) \cr
#' @param methode Method of the ODE-solver. Supported Boost \code{ODEINT}
#'   solvers (steppers) are:
#'   \code{'odeint::runge_kutta_cash_karp54'}
#'   \code{'odeint::runge_kutta_fehlberg78'} [default]
#'   \code{'odeint::runge_kutta_dopri5'}
#'   \code{'odeint::bulirsch_stoer'}
#'   \code{'odeint::adams_bashforth_[1|2|3|4|5|6|7|8]}
#'   \code{'odeint::adams_bashforth_moulton_[1|2|3|4|5|6|7|8]}
#'   without \code{odeint::}-prefix, \code{\link[deSolve]{ode}} method is
#'   assumed.
#' @param abstolint Absolute tolerance of the integration
#' @param reltolint Relative tolerance of the integration
#' @param verbose Logical controling if progress is printed to console.
#' @return The loglikelihood
#' @author Rampal S. Etienne & Bart Haegeman
#' @seealso \code{\link{DAISIE_ML_IW}}, \code{\link{DAISIE_loglik_CS}},
#' \code{\link{DAISIE_sim_cr}}
#' @references Valente, L.M., A.B. Phillimore and R.S. Etienne (2015).
#' Equilibrium and non-equilibrium dynamics simultaneously operate in the
#' Galapagos islands. Ecology Letters 18: 844-852.
#' @export DAISIE_loglik_IW
DAISIE_loglik_IW <- function(
  pars1,
  pars2,
  datalist,
  methode = "lsodes",
  abstolint = 1E-12,
  reltolint = 1E-10,
  verbose = FALSE
  )
{
  if(is.na(pars2[4]))
  {
    pars2[4] <- 0
  }
  if (is.null(datalist[[1]]$brts_table)) {
    datalist <- add_brt_table(datalist)
  }
  brts <- c(-abs(datalist[[1]]$brts_table[,'brt']),0)
  clade <- datalist[[1]]$brts_table[,'clade']
  event <- datalist[[1]]$brts_table[,'event']
  col <- datalist[[1]]$brts_table[,'col']
  pars1 <- as.numeric(pars1)
  if(length(pars1) == 5)
  {
     np <- datalist[[1]]$not_present
     if(is.null(np))
     {
        np <- datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2
     }
     if(is.null(np))
     {
        warning('Number of species not present is misspecified.')
        loglik <- NA
        return(loglik)
     }
     M <- length(datalist) - 1 + np
     pars1[6] <- M
  } else
  if(length(pars1) == 6) {
    M <- pars1[6]
  } else {
    warning("pars1 should contain 5 or 6 values.")
    loglik <- NA
    return(loglik)
  }

  ddep <- pars2[2]
  cond <- pars2[3]
  if (cond > 1) {
    stop('cond > 1 has not been implemented for the island-wide model.')
  }

  lac <- pars1[1]
  mu <- pars1[2]
  Kprime <- pars1[3]
  if(ddep == 0)
  {
    Kprime <- Inf
  }
  gam <- pars1[4]
  laa <- pars1[5]
  #l <- 0

  if (min(pars1) < 0)
  {
    message("One or more parameters are negative.")
    loglik <- -Inf
    return(loglik)
  }

  if(length(datalist) > 1) for(i in 2:length(datalist))
  {
    if(!datalist[[i]]$stac %in% c(0,2,4) & is.null(datalist[[i]]$all_colonisations)) {
      stop('IW does not work on data with unknown colonization times.')
    }
  }

  if((ddep == 1 | ddep == 11) & (ceiling(Kprime) < length(brts) - 2))
  {
    message('The proposed value of K is incompatible with the number of species
          in the clade. Likelihood for this parameter set will be set to -Inf. \n')
    loglik <- -Inf
    return(loglik)
  }

  if (ddep == 1 | ddep == 11)
  {
    lx <- min(1 + ceiling(Kprime), DDD::roundn(pars2[1]) )
  } else {
    lx <- DDD::roundn(pars2[1])
  }
  lxm <- min(lx,M + 1)
  lxe <- lx

  if(M * (1 - exp((min(brts) * gam))) > 0.2 * lxm) {
    message('With this colonization rate and system size setting, results may not be accurate.')
  }

  sysdimchange <- 1
  sysdim <- 1
  totdim <- lxm * lxe * sysdim
  probs <- rep(0, totdim)
  probs[1] <- 1
  loglik <- 0
  expandvec <- NULL
  for (k in 0:(length(brts) - 2))
  {
    if (isTRUE(identical(pars2[4], 3)))
    {
      message(paste('k = ',k ,', sysdim = ', sysdim, sep = ''))
    }
    dime <- list(lxm = lxm, lxe = lxe, sysdim = sysdim)
    if (sysdimchange == 1) {
      if (sysdim > 1) {
        dec2binmatk <- dec2binmat(log2(sysdim))
        l0ki <- create_l0ki(dec2binmatk, lxm, lxe, sysdim)
      } else if (sysdim == 1) {
        l0ki <- list(l0 = 0, ki = NULL)
      }
      sysdimchange <- 0
    }
    nndd <- nndivdep(lxm = lxm, lxe = lxe, sysdim = sysdim, Kprime = Kprime, k = k, M = M, l0 = l0ki$l0)
    parslist <- list(pars = pars1,k = k,ddep = ddep,dime = dime,l0ki = l0ki,nndd = nndd)
    iw_parms = DAISIE_IW_pars(parslist)
    if (startsWith(methode, "odeint::")) {
      probs <- .Call("daisie_odeint_iw", probs, brts[(k + 1):(k + 2)], iw_parms, methode, abstolint, reltolint)
    } else {
      y <- deSolve::ode(y = probs,
                        times = brts[(k + 1):(k + 2)],
                        func = DAISIE_loglik_rhs_IW,
                        parms = iw_parms,
                        rtol = reltolint,
                        atol = abstolint,
                        method = methode)
      probs <- y[2,2:(totdim + 1)]
    }
    cp <- checkprobs2(NA, loglik, probs, verbose); loglik = cp[[1]]; probs = cp[[2]]
    dim(probs) <- c(lxm, lxe, sysdim)

    if(k < (length(brts) - 2))
    {
      divdepfac <- nndd$divdepfac
      divdepfacplus1 <- nndd$divdepfacplus1
      mfac <- nndd$mfac
      oneminmfac <- nndd$oneminmfac
      if(event[k + 2] == 1) #colonization
      {
        #l <- l + 1
        if(is.na(col[k + 2]))
        {
          test_for_colonization <- TRUE
        } else
        {
          #test_for_colonization <- (max(event[which(clade == col[k + 2])]) > 1)
          # this tests whether the original clade has diversified, but this does not need to be the case
          # it could also be a case of anagenesis followed by colonization, so this test seems inappropriate
          # also, one should expect that colonizations in the data are all real colonizations; a colonization
          # which is not followed by speciation but by recolonization should only occur once (the last colonization)
          test_for_colonization <- TRUE
        }
        if(test_for_colonization) # new colonization or recolonization after speciation
        {
          probs2 <- array(0,dim = dim(probs))
          probs2[1:(lxm - 1),,] <- probs[2:lxm,,]
          #probs <- gam * divdepfac * probs[,,1:sysdim]
          probs <- gam * divdepfac * oneminmfac * probs[,,1:sysdim] +
            gam * divdepfacplus1 * mfac * probs2[,,1:sysdim]
          probs <- c(probs,rep(0,totdim))
          sysdim <- sysdim * 2
        } else # recolonization without speciation
        {
          tocollapse <- which(expandvec == col[k + 2])
          sr <- selectrows(sysdim,2^(tocollapse - 1))
          probs[,,sr[,1]] <- 0
          probs <- gam * divdepfac * probs[,,1:sysdim] #
          probs <- probs[,,sr[,2]]
          expandvec <- expandvec[-tocollapse]
          probs <- c(probs,rep(0,totdim/2))
        }
        expandvec <- c(expandvec,clade[k + 2])
        sysdimchange <- 1
      } else  # speciation
      {
        probs <- lac * divdepfac * probs[,,1:sysdim]
        if(event[k + 2] == 2) # first speciation in clade
        {
          tocollapse <- which(expandvec == clade[k + 2])
          sr <- selectrows(sysdim,2^(tocollapse - 1))
          probs <- probs[,,sr[,1]] + probs[,,sr[,2]]
          sysdim <- sysdim / 2
          dim(probs) <- c(lxm,lxe,sysdim)
          expandvec <- expandvec[-tocollapse]
          sysdimchange <- 1
        }
      }
      cp <- checkprobs2(NA, loglik, probs, verbose); loglik <- cp[[1]]; probs <- cp[[2]]
      totdim <- lxm * lxe * sysdim
      dim(probs) <- c(totdim,1)
    }
  }
  dim(probs) <- c(lxm, lxe, sysdim)
  expandedclades <- which(pracma::histc(clade, 1:length(clade))$cnt == 1)
  lexpandedclades <- length(expandedclades)
  status <- rep(0, lexpandedclades)
  if (lexpandedclades > 0) {
    for (i in lexpandedclades:1) {
      if (datalist[[1 + expandedclades[i]]]$stac == 2) {
        status[i] <- 1
      }
    }
  }
  if (length(status) > 0) {
    decstatus <- bin2dec(rev(status))
  } else {
    decstatus <- 0
  }
  loglik <- loglik + log(probs[1,1,1 + decstatus])

  if(cond > 0)
  {
    sysdim <- 1
    totdim <- lxm * lxe * sysdim
    dime <- list(lxm = lxm,lxe = lxe,sysdim = sysdim)
    probs <- rep(0,totdim)
    probs[1] <- 1
    l0ki <- list(l0 = 0,ki = NULL)
    nndd <- nndivdep(lxm = lxm,lxe = lxe,sysdim = sysdim,Kprime = Kprime,M = M,k = 0,l0 = l0ki$l0)
    parslist <- list(pars = pars1,k = 0,ddep = ddep,dime = dime,l0ki = l0ki,nndd = nndd)
    iw_parms <- DAISIE_IW_pars(parslist)
    if (startsWith(methode, "odeint::")) {
      probs <- .Call("daisie_odeint_iw", probs, c(min(brts),0), iw_parms, methode, abstolint, reltolint)
    } else {
      y <- deSolve::ode(y = probs,
                        times = c(min(brts), 0),
                        func = DAISIE_loglik_rhs_IW,
                        parms = iw_parms,
                        rtol = reltolint,
                        atol = abstolint,
                        method = methode)
      probs <- y[2,2:(totdim + 1)]
    }
    dim(probs) <- c(lxm, lxe, sysdim)
    logcond <- log1p(-probs[1,1,1])
    if(logcond == -Inf)
    {
        message('Parameters lead to probability of extinction of 1. Loglik is set to -Inf')
        loglik <- -Inf
    } else
    {
        loglik <- loglik - logcond
    }
  }
  print_parameters_and_loglik(pars = pars1,
                              loglik = loglik,
                              verbose = pars2[4],
                              type = 'island_loglik')
  return(as.numeric(loglik))
}
