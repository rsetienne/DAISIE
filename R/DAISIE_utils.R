#' Count the number of species
#'
#' @param datalistelement something
#' @keywords internal
#' @return A numeric value
countspecies = function(datalistelement)
{
  N = length(datalistelement$branching_times) - 1 + datalistelement$missing_species
}

counttype1 = function(datalistelement)
{
  N1 = 0
  if(length(datalistelement$type1or2) > 0)
  {
    N1 = (datalistelement$type1or2 == 1)
  }
}

countspeciestype1 = function(datalistelement)
{
  N1 = 0
  if(length(datalistelement$type1or2) > 0)
  {
    if(datalistelement$type1or2 == 1)
    {
      N1 = length(datalistelement$branching_times) - 1 + datalistelement$missing_species
    }
  }
}

countimmi = function(datalistelement)
{
  datalistelement$stac != 2
}

countstac = function(datalistelement,stac)
{
  return(datalistelement$stac == stac)
}

fconstr13 = function(x,pars1,x_E,age)
{
  lac = pars1[1]
  laa = pars1[5]
  ga = pars1[4]
  A = x - lac
  C = ga + laa + 2 * lac
  ff = (1 + A/C * (1 - exp(-C * age))) * exp(-A * age) - (1 - x_E)
  return(ff)
}

fconstr15 = function(x,pars1,x_E,x_I,age)
{
  lac = pars1[1]
  laa = pars1[5]
  A = x - lac
  B_c = -1/age * log(1 - x_I)
  ga = B_c - x - laa - lac
  C = ga + laa + 2 * lac
  ff = (1 + A/C * (1 - exp(-C * age))) * exp(-A * age) - (1 - x_E)
  return(ff)
}

calcMN = function(datalist,pars1)
{
  N = sum(unlist(lapply(datalist,countspecies)))
  if(is.null(datalist[[1]]$not_present))
  {
    M = datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2 + length(datalist) - 1
    if(!is.na(pars1[6]))
    {
      if(is.na(pars1[11]))
      {
        M = datalist[[1]]$not_present_type1 + sum(unlist(lapply(datalist,counttype1)))
      } else {
        M = M - max(0,DDD::roundn(pars1[11] * M))
      }
      N = sum(unlist(lapply(datalist,countspeciestype1)))
    }
  } else {
    M = datalist[[1]]$not_present + length(datalist) - 1
  }
  return(c(M,N))
}

DAISIE_eq = function(datalist,pars1,pars2)
{
  eqmodel = pars2[5]
  ddep = pars2[2]
  MN = calcMN(datalist,pars1)
  M = MN[1]
  N = MN[2]
  I = sum(unlist(lapply(datalist,countimmi)))
  rNM = N/M
  rIM = I/(M - I)
  rIN = I/(N - I)
  clado = pars1[1] * ((1 - N/pars1[3])^(ddep == 1 || ddep == 11)) * (exp(-N/pars1[3]))^(ddep == 2 || ddep == 21)
  ana = pars1[5]
  # Equilibrium based on deterministic model in terms of N
  if(eqmodel == 1)
  {
    immi = pars1[4] * ((1 - N/pars1[3])^(ddep == 11)) * (exp(-N/pars1[3]))^(ddep == 21)
    ext = clado + immi * (1/rNM - 1)
    pars1[2] = ext
  }
  # Equilibrium model based on deterministic model in terms of E and I
  if(eqmodel == 2) # Only eq for N
  {
    ext = pars1[2]
    immitot = 1/(1/rNM * 1/(ext - clado) - 1/(ana + clado + ext))
    immi = immitot / ((1 - N/pars1[3])^(ddep == 11) * (exp(-N/pars1[3]))^(ddep == 21))
    pars1[4] = immi
  }
  if(eqmodel == 3) # Only eq for E
  {
    immi = pars1[4] * ((1 - N/pars1[3])^(ddep == 11)) * (exp(-N/pars1[3]))^(ddep == 21)
    ext = clado + (ana + 2 * clado) * rIN
    pars1[2] = ext
  }
  if(eqmodel == 4) # Only eq for I
  {
    ext = pars1[2]
    immitot = (ext + ana + clado) * rIM
    immi = immitot / ((1 - N/pars1[3])^(ddep == 11) * (exp(-N/pars1[3]))^(ddep == 21))
    pars1[4] = immi
  }
  if(eqmodel == 5) # Eq for E and I
  {
    ext = clado + (ana + 2 * clado) * rIN
    immitot = (ext + ana + clado) * rIM
    immi = immitot / ((1 - N/pars1[3])^(ddep == 11) * (exp(-N/pars1[3]))^(ddep == 21))
    pars1[2] = ext
    pars1[4] = immi
  }
  if(eqmodel == 13) # Within x_E of equilibrium for E - diversity-dependence not implemented
  {
    x_E = pars2[10]
    x_I = pars2[11]
    age = datalist[[1]]$island_age
    pars1[2] = stats::uniroot(f = fconstr13,interval = c(pars1[1] + 1E-6, pars1[1] + 10),pars1 = pars1,x_E = x_E, age = age)$root
    ga_c = -1/age * log(1 - x_I) - pars1[1] - pars1[2] - pars1[5]
    if(pars1[4] < ga_c)
    {
      cat("The non-endemics do not satisfy the equilibrium criterion for these parameters.\n")
    }
  }
  if(eqmodel == 15) # Within x_E and x_I of equilibrium for both E and I - diversity-dependence not implemented
  {
    x_E = pars2[10]
    x_I = pars2[11]
    age = datalist[[1]]$island_age
    pars1[2] = stats::uniroot(f = fconstr15,interval = c(pars1[1] + 1E-6, pars1[1] + 10),pars1 = pars1,x_E = x_E, x_I = x_I, age = age)$root
    pars1[4] = -1/age * log(1 - x_I) - pars1[1] - pars1[2] - pars1[5]
  }
  return(pars1)
}

quantiles = function(probdist,probs)
{
  result = NULL
  cdf = cumsum(probdist[2,])
  for(i in 1:length(probs))
  {
    n = max(which(cdf <= probs[i]))
    x = probdist[1,n]
    if(cdf[n] == probs[i])
    {
      result[i] = x
    } else
      if(n < length(cdf))
      {
        result[i] = ((x + 1) * (probs[i] - cdf[n]) + x * (cdf[n + 1] - probs[i]))/(cdf[n + 1] - cdf[n])
      } else
      {
        result[i] = x
      }
  }
  names(result) = probs
  return(result)
}

antidiagSums = function(mat)
{
  dime = dim(mat)
  out = rep(0,sum(dime) - 1)
  nr = nrow(mat)
  nc = ncol(mat)
  for(i in 1:(nr + nc - 1))
  {
    rownums = min(i,nr):max(1,i - nc + 1)
    colnums = max(1,i - nr + 1):min(i,nc)
    for(j in 1:length(rownums))
    {
      out[i] = out[i] + mat[rownums[j],colnums[j]]
    }
  }
  return(out)
}

order_pars1 <- function(pars1)
{
  np <- names(pars1)
  correct_order <- c('max_area','proportional_peak_t','peak_sharpness','total_island_age','lac','mu_min','mu_max','K0','gam','laa')
  if(!is.null(np))
  {
    pars1ff <- pars1
    pars1ff[1] <- pars1[which(names(pars1) == 'max_area')]
    pars1ff[2] <- pars1[which(names(pars1) == 'proportional_peak_t')]
    pars1ff[3] <- pars1[which(names(pars1) == 'peak_sharpness')]
    pars1ff[4] <- pars1[which(names(pars1) == 'total_island_age')]
    pars1ff[5] <- pars1[which(names(pars1) == 'lac')]
    pars1ff[6] <- pars1[which(names(pars1) == 'mu_min')]
    pars1ff[7] <- pars1[which(names(pars1) == 'mu_max')]
    pars1ff[8] <- pars1[which(names(pars1) == 'K0')]
    pars1ff[9] <- pars1[which(names(pars1) == 'gam')]
    pars1ff[10] <- pars1[which(names(pars1) == 'laa')]
    pars1 <- pars1ff
    names(pars1) <- correct_order
  }
  return(pars1)
}

#' Samples from distribution when parameter is relaxed
#'
#' @param pars A vector of 5 elements with the model parameters
#' @param relaxed_par A string determining which parameter is relaxed
#' @param relaxed_rate_pars A list of two numbers, element one is the
#' distribution mean, element two is the distribution standard deviation (sd)
#'
#' @return A vector of parameters.
#' @keywords internal
sample_relaxed_rate <- function(pars,
                                relaxed_par) {

  mean <- which(c("cladogenesis",
                  "extinction",
                  "carrying_capacity",
                  "immigration",
                  "anagenesis") == relaxed_par)
  mean <- pars[mean]
  sd <- pars[6]
  if (relaxed_par == "cladogenesis") {
    pars[1] <- stats::rgamma(
      n = 1,
      shape = mean^2 / sd^2,
      scale = sd^2 / mean)
  }
  if (relaxed_par == "extinction") {
    pars[2] <- stats::rgamma(
      n = 1,
      shape = mean^2 / sd^2,
      scale = sd^2 / mean)
  }
  if (relaxed_par == "carrying_capacity") {
    pars[3] <- stats::rgamma(
      n = 1,
      shape = mean^2 / sd^2,
      scale = sd^2 / mean)
  }
  if (relaxed_par == "immigration") {
    pars[4] <- stats::rgamma(
      n = 1,
      shape = mean^2 / sd^2,
      scale = sd^2 / mean)
  }
  if (relaxed_par == "anagenesis") {
    pars[5] <- stats::rgamma(
      n = 1,
      shape = mean^2 / sd^2,
      scale = sd^2 / mean)
  }
  pars <- pars[1:5]
  return(pars)
}

#' Translate user-friendly ontogeny codes to numerics
#'
#' @inheritParams default_params_doc
#'
#' @return Numeric, 0 for null-ontogeny, 1 for beta function
#' @keywords internal
#' @examples translated_ontogeny <- DAISIE:::translate_island_ontogeny("const")
translate_island_ontogeny <- function(island_ontogeny) {

  if (island_ontogeny == "const" || island_ontogeny == 0) {
    island_ontogeny <- 0
  }
  if (island_ontogeny == "beta" || island_ontogeny == 1) {
    island_ontogeny <- 1
  }
  return(island_ontogeny)
}

#' Translate user-friendly sea-level codes to numerics
#'
#' @inheritParams default_params_doc
#'
#' @return Numeric, 0 for null-sea-level, 1 for sine function
#' @keywords internal
#' @examples translated_sea_level <- DAISIE:::translate_sea_level("const")
translate_sea_level <- function(sea_level) {

  if (sea_level == "const" || sea_level == 0) {
    sea_level <- 0
  }

  if (sea_level == "sine" || sea_level == 1) {
    sea_level <- 1
  }
  return(sea_level)
}

#' Calculates the species on the island initially when \code{nonoceanic_pars[1]
#' != 0}
#'
#' @param prob_samp probability of a species being sampled
#' from the mainland pool
#' @param prob_nonend probability of a species sampled being non-endemic
#' @param mainland_n number of species in the mainland pool
#'
#' @return A list of non-endemic species, endemic species and the new
#' mainland species pool
#' @examples DAISIE:::DAISIE_nonoceanic_spec(
#' prob_samp = 0.1,
#' prob_nonend = 0.9,
#' mainland_n = 1000)
#' @keywords internal
DAISIE_nonoceanic_spec <- function(prob_samp, prob_nonend, mainland_n) {
  testit::assert(prob_samp <= 1)
  testit::assert(prob_samp >= 0)
  testit::assert(prob_nonend <= 1)
  testit::assert(prob_nonend  >= 0)
  testit::assert(length(mainland_n) > 0)
  if (prob_samp != 0) {
    prob_not_samp <- 1 - prob_samp
    prob_nonend <- prob_samp * prob_nonend
    prob_end <- 1 - (prob_not_samp + prob_nonend)
    num_native_spec <- sample(1:3, length(1:mainland_n),
                              replace = TRUE,
                              prob = c(prob_not_samp, prob_nonend, prob_end))
    init_nonend_spec_vec <- sample(1:mainland_n,
                                   length(which(num_native_spec == 2)),
                                   replace = FALSE)
    new_source_pool <- setdiff(1:mainland_n, init_nonend_spec_vec)
    init_end_spec_vec <- sample(new_source_pool,
                                length(which(num_native_spec == 3)),
                                replace = FALSE)
    mainland_spec <- setdiff(1:mainland_n, init_end_spec_vec)
    testit::assert(sum(length(which(num_native_spec == 1)),
                       length(which(num_native_spec == 2)),
                       length(which(num_native_spec == 3)))
                   == sum(mainland_n))
    init_nonend_spec <- length(init_nonend_spec_vec)
    init_end_spec <- length(init_end_spec_vec)
    if (length(mainland_spec) == 0) {
      mainland_spec <- 0
    }
  } else {
    init_nonend_spec <- 0
    init_end_spec <- 0
    init_nonend_spec_vec <- integer(0)
    init_end_spec_vec <- integer(0)
    if(mainland_n != 0){
      mainland_spec <- seq(1, mainland_n, 1)
    } else {
      mainland_spec = c()
    }
  }
  return(list(init_nonend_spec = init_nonend_spec,
              init_end_spec = init_end_spec,
              init_nonend_spec_vec = init_nonend_spec_vec,
              init_end_spec_vec = init_end_spec_vec,
              mainland_spec = mainland_spec))
}

#' Update internal Gillespie bookeeping objects
#'
#' @param stt_table A species=through-time table.
#' @param total_time Simulated amount of time.
#' @param timeval Current time of simulation.
#' @param mainland_spec A vector with the numeric IDs of the mainland species
#' (i.e. potential colonizers).
#' @param island_spec A matrix with the species on the island (state of the
#' system at each time point).
#'
#' @return A named list with the updated input arguments at time of simulation.
#'
#' @noRd
DAISIE_spec_tables <- function(stt_table,
                               total_time,
                               timeval,
                               nonoceanic_sample,
                               island_spec,
                               maxspecID) {
  init_nonend_spec <- nonoceanic_sample$init_nonend_spec
  init_end_spec <- nonoceanic_sample$init_end_spec
  mainland_spec <- nonoceanic_sample$mainland_spec
  stt_table[1, ] <- c(total_time,
                      init_nonend_spec,
                      init_end_spec,
                      0)
  if (init_nonend_spec != 0) {
    for (i in seq_along(1:init_nonend_spec)) {
      island_spec <- rbind(island_spec,
                           c(nonoceanic_sample$init_nonend_spec_vec[i],
                             nonoceanic_sample$init_nonend_spec_vec[i],
                             timeval,
                             "I",
                             NA,
                             NA,
                             NA))
    }
  }
  if (init_end_spec != 0) {
    for (j in seq_along(1:init_end_spec)) {
      maxspecID <- maxspecID + 1
      island_spec <- rbind(island_spec,
                           c(maxspecID,
                             nonoceanic_sample$init_end_spec_vec[j],
                             timeval,
                             "A",
                             NA,
                             NA,
                             NA))
    }
  }
  return(list(stt_table = stt_table,
              init_nonend_spec = init_nonend_spec,
              init_end_spec = init_end_spec,
              mainland_spec = mainland_spec,
              island_spec = island_spec,
              maxspecID = maxspecID))
}

#' Add a column to a data frame
#'
#' @param df data frame to add the column to
#' @param position location in data frame where to insert the column.
#' Position can also be a name of a column
#' @param column_to_insert the elements of the column to insert. If
#' the column has a name, this name will be copied into the data frame.
#' Id is does not have a name, it will get the name "nc".
#'
#' @return A data frame with the column inserted
add_column_to_dataframe <- function(df, position, column_to_insert) {
  if(is.character(position)) {
    position <- which(names(df) == position)
  }
  df <- data.frame(df[1:position],
                   nc = column_to_insert,
                   df[(position + 1):ncol(df)])
  names(df)[names(df) == 'nc'] <- names(column_to_insert)
  return(df)
}
