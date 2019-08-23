#' Simulates an island using DAISIE_sim_core and then calculates the half-life
#' and repeats this across a parameter space
#'
#' @inheritParams default_params_doc
#' @param lac numeric vector with all values of cladogenesis giving the 
#' parameter space of cladogenesis.
#' @param mu numeric vector with all values of extinction giving the parameter
#' space of extinction.
#' @param K numeric vector with all values of the carrying capacity giving the 
#' parameter space of the carrying capacity.
#' @param gam numeric vector with all values of immigration giving the parameter
#' space of immigration.
#' @param laa numeric vector with all values of anagenesis giving the parameter
#' space of anagenesis.
#' @param ssr numeric vector with all the supersaturation ratios (ssr) giving
#' the parameter space of the supersaturation ratio.
#' @param ddmodel_sim numeric vector specifying which diversity-depedent models
#' should be used. \code{ddmodel_sim == 0} is a diversity-independent (constant
#' rate) model; \code{ddmodel_sim == 1} is a linear diversity-depedent 
#' cladogenesis model; \code{ddmodel_sim == 11} is a linear diversity-
#' dependent cladogenesis and immigration model; \code{ddmodel_sim == 2} is an
#' exponential cladogenesis model; \code{ddmodel_sim == 3} is a linear 
#' diversity-dependent extinction model.
#' @param x_s numeric vector with all the values of non-oceanic species sampling
#' giving the paramter space of non-oceanic species sampling (x_s is the
#' first element of the nonoceanic parameter vector from DAISIE_sim).
#' @param x_nonend numeric vector with all the values of the non-endemic
#' sampling giving the parameter space of the non-endemic sampling (x_nonend
#' is the second element of the nonoceanic parameter vector from DAISIE_sim).
#' @param verbose boolean determining whether to print
#'
#' @return A list of numerics which are the relaxation half-lives.
#' @export
#'
#' @examples
DAISIE_sim_half_life <- function(time,
                                 replicates,
                                 lac,
                                 mu,
                                 K,
                                 gam, 
                                 laa,
                                 ssr,
                                 divdepmodel,
                                 ddmodel_sim,
                                 island_type,
                                 x_s,
                                 x_nonend,
                                 verbose) {
  testit::assert(is.numeric(time))
  testit::assert(is.numeric(replicates))
  testit::assert(is.numeric(lac))
  testit::assert(is.numeric(mu))
  testit::assert(is.numeric(K))
  testit::assert(is.numeric(gam))
  testit::assert(is.numeric(laa))
  testit::assert(is.numeric(ssr))
  testit::assert(divdepmodel == "CS" || divdepmodel == "IW")
  testit::assert(is.numeric(ddmodel_sim))
  testit::assert(island_type == "oceanic" || island_type == "nonoceanic")
  testit::assert(is.vector(x_s) || is.null(x_s))
  testit::assert(is.vector(x_nonend) || is.null(x_nonend))
  if (is.null(x_s) && is.null(x_nonend)) {
    x_s <- 0
    x_nonend <- 0
  }
  if (divdepmodel == "IW") {
    grid <- expand.grid(c(lac), c(mu), c(K), c(gam), c(laa), c(ssr),
                        c(ddmodel_sim), c(x_s), c(x_nonend))
    names(grid) <- c("lac", "mu", "K", "gam", "laa", "ssr",
                     "ddmodel_sim", "x_s", "x_nonend")
    if (island_type == "nonoceanic") {
      mainland_n <- (grid["K"] * grid["ssr"] / grid["x_s"])
    } else {
      mainland_n <- 1000 #change this to a proper mainland_n
    }
    names(mainland_n) <- "mainland_n"
    grid <- cbind(grid, mainland_n)
    island_replicates <- list()
    half_life <- list()
    for (i in 1:nrow(grid)) {
      island_replicates[[i]] <- list()
      half_life[[i]] <- list()
      if (verbose == TRUE) {
        print(paste("Parameter set", i)) 
      }
      for (rep in 1:replicates) {
        if (verbose == TRUE) {
          print(paste("Replicate", rep))
          
        }
        island_replicates[[i]][[rep]] <- DAISIE_sim_core(
          time = time,
          mainland_n = unname(unlist(grid[i, "mainland_n"])),
          pars = unname(unlist(grid[i, c("lac", "mu", "K", "gam", "laa")])),
          ddmodel_sim = unname(unlist((grid[i, "ddmodel_sim"]))),
          island_type = island_type,
          nonoceanic = unname(unlist(grid[i, c("x_s", "x_nonend")])))
      }
      half_life[[i]] <- DAISIE_calc_half_life(
        island_replicates = island_replicates[[i]],
        mainland_n = unlist(grid[i, "mainland_n"]),
        pars = unname(unlist(grid[i, c("lac", "mu", "K", "gam", "laa")])),
        island_type = island_type,
        divdepmodel = divdepmodel)
    }
  }
  if (divdepmodel == "CS") {
    grid <- expand.grid(c(lac), c(mu), c(K), c(gam), c(laa), c(ssr),
                        c(ddmodel_sim), c(x_s), c(x_nonend))
    names(grid) <- c("lac", "mu", "K", "gam", "laa", "ssr",
                     "ddmodel_sim", "x_s", "x_nonend")
    if (island_type == "nonoceanic") {
      mainland_n <- (grid["K"] * grid["ssr"] / grid["x_s"])
    } else {
      mainland_n <- 1000 #change this to a proper mainland_n
    }
    names(mainland_n) <- "mainland_n"
    grid <- cbind(grid, mainland_n)
    island_replicates <- list()
    half_life <- list()
    for (i in 1:nrow(grid)) {
      island_replicates[[i]] <- list()
      half_life[[i]] <- list()
      if (verbose == TRUE) {
        print(paste("Parameter set", i))
        
      }
      for (rep in 1:replicates) {
        if (verbose == TRUE) {
          print(paste("Replicate", rep))
          
        }
        island_replicates[[i]][[rep]] <- list()
        full_list <- list()
        for (m_spec in 1:grid[i, "mainland_n"]) {
          full_list[[m_spec]] <- DAISIE_sim_core(
            time = time,
            mainland_n = 1,
            pars = unname(unlist(grid[i, c("lac", "mu", "K", "gam", "laa")])),
            ddmodel_sim = unname(unlist((grid[i, "ddmodel_sim"]))),
                                                 island_type = island_type,
                                                 nonoceanic = unname(unlist(grid[i, c("x_s", "x_nonend")])))
          }  
      island_replicates[[i]][[rep]] <- full_list
      }
    }
    half_life[[i]] <- DAISIE_calc_half_life(
      island_replicates = island_replicates[[i]],
      mainland_n = unlist(grid[i, "mainland_n"]),
      pars = unname(unlist(grid[i, c("lac", "mu", "K", "gam", "laa")])),
      island_type = island_type,
      divdepmodel = divdepmodel)
    }
  return(half_life)
}