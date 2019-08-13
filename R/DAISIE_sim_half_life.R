DAISIE_sim_half_life <- function(time, 
                                 replicates, 
                                 lac, 
                                 mu, 
                                 K, 
                                 gam, 
                                 laa, 
                                 ssr, 
                                 ddlac, 
                                 ddmu, 
                                 ddgam, 
                                 island_type,
                                 x_s,
                                 x_nonend) {
  testit::assert(is.numeric(time))
  testit::assert(is.numeric(replicates))
  testit::assert(is.numeric(lac))
  testit::assert(is.numeric(mu))
  testit::assert(is.numeric(K))
  testit::assert(is.numeric(gam))
  testit::assert(is.numeric(laa))
  testit::assert(is.numeric(ssr))
  testit::assert(is.numeric(ddlac))
  testit::assert(ddlac == 0 || ddlac == 1)
  testit::assert(ddmu == 0 || ddmu == 1 || ddmu == 2)
  testit::assert(length(ddmu) == 1 || length(ddmu) == 2 || length(ddmu) == 3)
  testit::assert(length(ddgam) == 1 || length(ddgam) == 2)
  testit::assert(island_type == "oceanic" || island_type == "nonoceanic")
  testit::assert(length(x_s) == 1 ||
                   length(x_s) == 3 || 
                   is.null(x_s))
  testit::assert(length(x_nonend) == 1 ||
                   length(x_nonend) == 3 || 
                   is.null(x_nonend))
  grid <- expand.grid(c(lac), c(mu), c(K), c(gam), c(laa), c(ssr), c(ddlac), c(ddmu), c(ddgam), c(x_s), c(x_nonend))
  names(grid) <- c("lac", "mu", "K", "gam", "laa", "ssr", "ddlac", "ddmu", "ddgam", "x_s", "x_nonend")
  mainland_n <- (grid["K"] * grid["ssr"] / grid["x_s"])
  names(mainland_n) <- "mainland_n"
  grid <- cbind(grid, mainland_n)
  island_replicates <- list()
  half_lives <- list()
  for (i in 1:nrow(grid)) {
    island_replicates[[i]] <- list()
    half_lives[[i]] <- list()
    print(paste("Parameter set", i))
    for (rep in 1:replicates) {
      print(paste("Replicate", rep))
      island_replicates[[i]][[rep]] <- DAISIE_sim_core(time = time,
                                                       mainland_n = unname(unlist(grid[i, 12])),
                                                       pars = unname(unlist(grid[i, 1:5])), 
                                                       ddmodel = unname(unlist((grid[i, 7:9]))), 
                                                       island_type = island_type, 
                                                       nonoceanic = unname(unlist(grid[i, 10:11])))
    }
    half_lives[[i]] <- DAISIE_calc_half_life(island_replicates = island_replicates[[i]], pars = unname(unlist(grid[i, 1:5])), divdepmodel = "IW") 
  }
  return(island_replicates)
}