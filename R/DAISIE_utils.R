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

is.odd <- function(x)
{ 
  res <- x %% 2 != 0
  return(res)
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

#' Translate user-friendly ontogeny codes to numerics
#'
#' @inherit DAISIE_sim
#'
#' @return Numeric, 0 for null-ontogeny, 1 for linear decrease and 
#' 2 for beta function
#' @export
#' @examples translate_island_ontogeny("const")
translate_island_ontogeny <- function(island_ontogeny) {
 
  if (island_ontogeny == "const" || island_ontogeny == 0) {
    island_ontogeny <- 0
  }
   
  if (island_ontogeny == "linear" || island_ontogeny == 1) {
    island_ontogeny <- 1
  }
   
  if (island_ontogeny == "beta" || island_ontogeny == 2) {
    island_ontogeny <- 2 
  }
  return(island_ontogeny)
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


#' Determine if list has only numerical values.
#' 
#'
#' @param x Object to determine
#'
#' @return Boolean indicating if object is list with only numerical values
#' @note do not forget: NAs are removed from a list!
#' @examples 
#'   testit::assert(
#'     DAISIE:::is_numeric_list(
#'       x = list(char = "character", numerical = 1)
#'     ) == FALSE
#'   )
#'   
#'   testit::assert(
#'     DAISIE:::is_numeric_list(
#'       x = list(numerical_1 = 1, numerical_2 = 2)
#'     ) == TRUE
#'   )
is_numeric_list <- function(x) {
  is.list(x) && is.numeric(unlist(x))
}
