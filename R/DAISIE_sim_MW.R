#' @title Simulate multiple islands based on hyperparameters that describe relationships
#' between area, isolation and local parameters
#' @description This function simulates islands of given age, area and isolation.
#' A list of islands/archipelagos with this information is provided as a data frame
#' object (in the example,
#' the archipelago_data from Valente et al 2020 Nature).
#' For each island, local parameters are first calculated based on hyperparameters describing
#' the dependency of rates on island area and isolation. Simulations are then run for each island,
#' from island birth until the given island age based on the
#' local parameters for each island (using the DAISIE_sim function).
#'
#' Returns R list object that contains the simulated islands.
#'
#' @param archipelago_data Object (data.frame) containing a table with columns with
#' the following names:
#' Archipelago: names of the archipelagos or islands to be simulated;
#' Area - the area of the archipelago/island;
#' Age: Age of the oldest island in the archipelago, or age of the island if single island.
#' Distance - distance to the mainland (or other applicable isolation metric).
#' As an example, the dataset from Valente et al 2020 for 41 archipelagos is provided in data(archipelago_data).
#' @param M The size of the mainland pool, i.e the number of species that can
#' potentially colonize the island(s).
#' @param pars Contains the model hyperparameters:
#' \cr \cr \code{pars[1]} corresponds to lambda^c0 (initial cladogenesis rate)
#' \cr \code{pars[2]} corresponds to y (dependency of cladogenesis on area)
#' 0 - no dependency, positive value - cladogenesis increases with area,
#' negative value - cladogenesis decreases with area
#' \cr \code{pars[3]} corresponds to mu_0 (initial extinction rate)
#' \cr \code{pars[4]} corresponds to x (dependency of extinction on area)
#' 0 - no dependency, positive value - extinction increases with area,
#' negative value - extinction decreases with area.
#' \cr \code{pars[5]} corresponds to K_0  (initial
#' carrying capacity). Set K=Inf for non-diversity dependence.
#' \cr \code{pars[6]} corresponds to z
#' (dependency of K on area). 0 - no dependency, positive value - K increases with area,
#' negative value - K decreases with area.
#' \cr \code{pars[7]} corresponds to gamma_0 (initial immigration rate)
#' \cr \code{pars[8]} corresponds to alpha (dependency of immigration on island isolation).
#' 0 - no dependency, positive value - immigration increases with distance,
#' negative value - immigration decreases with distance
#' \cr \code{pars[9]} corresponds to lambda^a_0 (initial anagenesis rate)
#' \cr \code{pars[10]} corresponds to beta (dependency of anagenesis on island isolation).
#' 0 - no dependency,positive value - anagenesis increases with distance, negative value
#' - anagenesis decreases with distance.
#' \cr \code{pars[11]} corresponds to d0, parameter that describes
#' interactive/additive effect of isolation and area on cladogenesis (power models);
#' or that describes the shape of the sigmoidal relationship between
#' isolation and a parameter (sigmoidal models).
#' @param replicates Number of island replicates to be simulated per island in
#' the table archipelago_data. If there are 5 islands in archipelago_data, and number of
#' replicates is set to 10, 50 islands will be simulated.
#' @param divdepmodel Option divdepmodel = 'CS' runs model with clade-specific
#' carrying capacity, where diversity-dependence operates only within single
#' clades, i.e. only among species originating from the same mainland
#' colonist.\cr Option divdepmodel = 'IW' runs model with island-wide carrying
#' capacity, where diversity-dependence operates within and among clades.
#' @param distance_dep Shape of dependency of parameter with distance. Can be of type \cr 'power'
#' or type \cr 'sigmoidal'.
#' @param cladogenesis_dep When distance_dep = 'power' is specified, cladogenesis_dep specifies
#' the type of effect of area and isolation on cladogenesis. The possible options are:
#' \cr 'NULL' (no additive or interactive effect of area on isolation,
#' same as models M1-M14 in Valente et al 2020)
#' \cr 'additive' (additive effect of area and isolation on cladogenesis,
#' same as model M15 in Valente et al 2020)
#' \cr 'interactive' (interactive effect of area and isolation on cladogenesis,
#' same as models M16 and M19 in Valente et al 2020)
#' \cr 'interactive1' (interactive effect of area and isolation on cladogenesis,
#' same as model M17 in Valente et al 2020)
#' \cr 'interactive2' (interactive effect of area and isolation on cladogenesis,
#' same as model M18 in Valente et al 2020)
#' @param sigmoidal_par When distance_dep = 'sigmoidal', sigmoidal_par specifies
#' to which parameter the sigmoidal relationship with distance is applied.
#' Options: \cr 'cladogenesis'\cr 'anagenesis' \cr 'colonisation'.

#' @return For each island listed in the rows of archipelago_data, a given number of islands
#' is simulated depending on the number of replicates specified. The simulations for
#' each island in archipelago_data are an element of the list, which can be called using [[x]].
#' Individual replicates for an island can be called with [[x]][[x]]
#' For example if the object is called global_sims, the 1st replicate for the first island
#' can be called using global_sims[[1]][[1]].
#' The 3rd replicate of the 4th island is called with global_sims[[4]][[3]]
#' Each of the island replicates is a list in itself. The first (e.g.
#' global_sims[[x]][[x]][[1]]) element of that list has the following
#' components: \cr \code{$island_age} - the island or archipelago age \code{$not_present}
#' - the number of mainland lineages that are not present on the island \cr
#' \code{$stt_all} - STT table for all species on the archipelago/island (nI - number of
#' non-endemic species; nA - number of anagenetic species, nC - number of
#' cladogenetic species, present - number of independent colonisations present)
#' \cr \code{$area} Area of the island or archipelago \cr \code{$distance} Distance to the mainland
#' \cr \code{$name} Name of the island or archipelago \cr
#'
#' The subsequent elements of the list each contain information on a single
#' colonist lineage on the island and has 4 components:\cr
#' \code{$branching_times} - island age and stem age of the population/species
#' in the case of Non-endemic, Non-endemic_MaxAge and Endemic with no close extant
#' relatives on the islands
#' species. For cladogenetic species these should be island age and branching
#' times of the radiation including the stem age of the radiation.\cr
#' \code{$stac} - the status of the colonist \cr *
#' Endemic: 2 \cr * Endemic&Non_Endemic: 3 \cr * Non_endemic: 4 \cr
#' \code{$missing_species} - number of island species that were not sampled for
#' particular clade (only applicable for endemic clades) \cr
#' @author Luis Valente, Albert Phillimore, Rampal Etienne
#' @seealso \code{\link{DAISIE_sim_cr}()}.
#' @references Valente, LM, Phillimore AB, Melo M, Warren B, Clegg S, Havenstein K,
#' Tiedemann R, Illera JC, Thebaud C, Aschenbach T and  Etienne RS (2020).
#' A simple dynamic model explain island bird diversity worldwide. Nature, 579, 92-96.
#' @keywords models
#' @examples
#' ## Simulate 10 replicates for the each of the 41 archipelagos in the archipelago_data table,
#' ##based on the hyperparameters of the M19 model (preferred model in Valente et al 2020)
#' ## and the age, area and isolation values given in the archipelago_data table.
#' ## Mainland pool size of 1000, clade-specific carrying capacity. The M19 model is a
#' ## power model with an interactive effect of area and isolation on cladogenesis.
#' ##
#' ## data(archipelago_data)
#' ## result <- DAISIE_sim_MW(
#' ## archipelago_data = archipelago_data,
#' ## M = 1000,
#' ## pars = c(0.040073803,	0,	1.945656546,	0.150429656,	Inf,	0,	67.25643672,
#' ## 0.293635061,	0.059096872,	0.382688527,	0.026510781),
#' ## replicates = 10,
#' ## distance_dep = 'power',
#' ## cladogenesis_dep = 'interactive',
#' ## sigmoidal_par = 'NULL',
#' ## divdepmodel = 'CS')
#' @export DAISIE_sim_MW

DAISIE_sim_MW <- DAISIE_sim_global <- function(
  archipelago_data,
  M,
  pars,
  replicates,
  divdepmodel = 'CS',
  distance_dep = 'power',
  cladogenesis_dep = 'NULL',
  sigmoidal_par = 'NULL')
{

  check_length_pars <- function(length_lim_pars)
  {
    if(length(pars) < length_lim_pars) {
      stop('Insufficient number of parameters specified under pars')
    }
  }

  Archipelago <- archipelago_data$Archipelago
  Area <- archipelago_data$Area
  Distance <- archipelago_data$Distance
  Age <- archipelago_data$Age

  lamc_0 <- pars[1]
  y <- pars[2]
  mu_0 <- pars[3]
  x <- pars[4]
  K_0 <- pars[5]
  z <- pars[6]
  gam_0timesM <- pars[7]
  alpha <- pars[8]
  lama_0 <- pars[9]
  beta <- pars[10]

  if(distance_dep == 'power') {
    lambda_c <- lamc_0 * Area^ y
    mu <- mu_0 * Area^ -x
    K <- K_0 * Area^ z
    gamma <- (gam_0timesM * Distance^ -alpha)/M
    lamda_a <- lama_0 * Distance^ beta

    ## M15 model
    if(cladogenesis_dep == 'additive')
    {
      check_length_pars(11)
      d0 <- pars[11]
      lambda_c <- lamc_0 * (Area^y) * (Distance^d0)
    }

    ## M16 and M19 models
    if(cladogenesis_dep == 'interactive')
    {
      check_length_pars(11)
      d0 <- pars[11]
      lambda_c <- lamc_0 * Area^(y + d0 * log(Distance))
    }

    ## M17 model
    if(cladogenesis_dep == 'interactive1')
    {
      check_length_pars(11)
      d0 <- pars[11]
      lambda_c <- lamc_0 * Area^(y + Distance/d0)
    }

    ## M18 model
    if(cladogenesis_dep == 'interactive2')
    {
      check_length_pars(11)
      d0 <- pars[11]
      lambda_c <- lamc_0 * Area^(y + 1/(1 + d0/Distance))
    }
  }

  if(distance_dep == 'sigmoidal') {
    check_length_pars(11)
    mu <- mu_0 * Area^ -x
    K <- K_0 * Area^ z

    if(sigmoidal_par == 'colonisation'){
      lambda_c <- lamc_0 * Area ^ y
      lamda_a <- lama_0 * Distance ^ beta

      kg <- pars[7]
      xg <- pars[8]
      d_0 <- pars[11]

      g <- function(d,x,d0,k) {k - k * (d/d0)^x/(1 + (d/d0)^x)}

      gamma <- g(Distance,xg,d_0,kg)/M
    }

    if(sigmoidal_par == 'anagenesis'){
      lambda_c <- lamc_0 * Area^ y
      gamma <- (gam_0timesM * Distance^ -alpha)

      kf <- pars[9]
      xf <- pars[10]
      d0_f <- pars[11]

      f <- function(d,x,d0,k) {k * (d/d0)^x/(1 + (d/d0)^x)}
      lamda_a <- f(Distance,xf,d0_f,kf)
    }

    if(sigmoidal_par == 'cladogenesis'){
      lamda_a <- lama_0 * Distance^ beta
      gamma <- (gam_0timesM * Distance^ -alpha)/M

      kf <- pars[1]
      xf <- pars[2]
      d0_f <- pars[11]

      f <- function(d,x,d0,k) {k * (d/d0)^x/(1 + (d/d0)^x)}
      lambda_c <- f(Distance,xf,d0_f,kf)
    }
  }

  Archipelago <- as.character(Archipelago)
  archipelago_specific <- cbind(Archipelago,Area,Distance,Age,lambda_c,mu,K,gamma,lamda_a)
  print('Input data and local parameters for each island:')
  print(archipelago_specific)
  print('Start of the simulations:')

  global_sim <- list()
  for (i in 1:nrow(archipelago_specific)) {
    pars1 <- as.numeric(archipelago_specific[i,5:9])
    the_archipelago <- as.character(archipelago_specific[i,'Archipelago'])
    island_age <- as.numeric(archipelago_specific[i,'Age'])
    area <- archipelago_data[i,"Area"]
    distance <- archipelago_data[i,"Distance"]

    global_sim[[i]] <- DAISIE_sim(time = island_age,
                                  M = M,
                                  pars = pars1,
                                  replicates = replicates,
                                  divdepmodel = 'CS',
                                  plot_sims = FALSE)

    for(a in 1:length(global_sim[[i]])) {
      global_sim[[i]][[a]][[1]]$area <- area
      global_sim[[i]][[a]][[1]]$distance <- distance
      global_sim[[i]][[a]][[1]]$name <- the_archipelago
    }

    print(paste('Simulation of',as.character(archipelago_specific[i,'Archipelago']),'for',island_age,'Myrs each completed',sep=' '))
  }
  return(global_sim)
}
