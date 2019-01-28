context("DAISIE_ext_test")

test_that("test expected species vs simulated with extinction", {
  skip("WIP")
  expect_silent(
    #' Internal function of the DAISIE simulation
    #' @param time simulated amount of time
    #' @param mainland_n number of mainland species, that
    #'   is, the number of species that can potentially colonize the island.
    #'   If \code{\link{DAISIE_sim}} uses a clade-specific diversity dependence,
    #'   this value is set to 1. 
    #'   If \code{\link{DAISIE_sim}} uses an island-specific diversity dependence,
    #'   this value is set to the number of mainland species.
    #' @param pars a numeric vector:
    #' \itemize{
    #'   \item{[1]: cladogenesis rate}
    #'   \item{[2]: extinction rate}
    #'   \item{[3]: carrying capacity}
    #'   \item{[4]: immigration rate}
    #'   \item{[5]: anagenesis rate}
    #' }
    #' @param Apars a named list containing area parameters as created by create_area_params:
    #' \itemize{
    #'   \item{[1]: maximum area}
    #'   \item{[2]: value from 0 to 1 indicating where in the island's history the 
    #'   peak area is achieved}
    #'   \item{[3]: sharpness of peak}
    #'   \item{[4]: total island age}
    #' }
    #' @param Epars a numeric vector:
    #' \itemize{
    #'   \item{[1]: minimum extinction when area is at peak}
    #'   \item{[2]: extinction rate when current area is 0.10 of maximum area}
    #' }
    #' @param island_ontogeny a string describing the type of island ontogeny. Can be \code{NULL},
    #' \code{"quadratic"} for a beta function describing area through time,
    #'  or \code{"linear"} for a linear function
    DAISIE_exinction_test <- function(
      time,
      mainland_n,
      pars,
      Apars = NULL,
      Epars = NULL,
      island_ontogeny = NULL
    ) {
      timeval <- 0
      totaltime <- time
      lac <- pars[1]
      mu <- pars[2]
      K <- pars[3]
      gam <- pars[4]
      laa <- pars[5]
      
      extcutoff <- max(1000, 1000 * (laa + lac + gam))
      ext_multiplier <- 0.5
      stt <- matrix(ncol = 2)
      # if(pars[4] == 0) 
      # {
      #   stop('Rate of colonisation is zero. Island cannot be colonised.')
      # }  
      
      if (are_area_params(Apars) && is.null(island_ontogeny)){
        stop("Apars specified for contant island_ontogeny. Set Apars to NULL")
      }
      
      if (!is.null(island_ontogeny) && island_ontogeny != "linear" && island_ontogeny != "quadratic") {
        stop("Please select valid island ontogeny model. Options are no ontogeny: NULL, 'linear' or 'quadratic'.")
      }
      
      mainland_spec <- seq(1, mainland_n, 1)
      maxspecID <- mainland_n
      
      island_spec = matrix(ncol = 7, nrow = 1000)
      island_spec[,4] = "I"
      stt_table <- matrix(ncol = 4)
      colnames(stt_table) <- c("Time","nI","nA","nC")
      stt_table[1,] <- c(totaltime,0,0,0)
      
      # Pick t_hor (before timeval, to set Amax t_hor)
      t_hor <- get_t_hor(timeval = 0,
                       totaltime = totaltime,
                       Apars = Apars,
                       ext_multiplier = ext_multiplier,
                       island_ontogeny = island_ontogeny,
                       t_hor = NULL)
      
      #### Start Gillespie ####
      while (timeval < totaltime) {
        if (timeval < t_hor) {
          rates <- update_rates(timeval = timeval, totaltime = totaltime, gam = gam,
                                mu = mu, laa = laa, lac = lac, Apars = Apars,
                                Epars = Epars, island_ontogeny = island_ontogeny,
                                extcutoff = extcutoff, K = K,
                                island_spec = island_spec, mainland_n, t_hor)
          if (is.na(timeval) == T) {
            timeval <- totaltime
          } else {
            timeval <- pick_timeval(rates, timeval)
          }
          # Determine event
          # If statement prevents odd behaviour of sample when rates are 0
          if (is.null(island_ontogeny)) {
            possible_event <- sample(
              1:4, 1, 
              prob = c(rates[[1]], rates[[2]], 
                       rates[[3]], rates[[4]]), 
              replace = FALSE)
          } else if (sum(rates[[1]], rates[[2]], 
                         rates[[3]], rates[[4]], 
                         rates[[5]]) > 0) {
            possible_event <- sample(1:5, 1, prob = c(rates[[1]], rates[[2]], 
                                                      rates[[3]], rates[[4]], 
                                                      (rates[[5]] - rates[[2]])),
                                     replace = FALSE)
          }
          if (is.nan(timeval) == T) {
            timeval <- totaltime
          }
          
          if (timeval < totaltime) {
            # Run event
            
            
            new_state <- DAISIE_sim_update_state(timeval = timeval,
                                                 possible_event = possible_event,
                                                 maxspecID = maxspecID,
                                                 mainland_spec = mainland_spec,
                                                 island_spec = island_spec)
            
            island_spec <- new_state$island_spec
            maxspecID <- new_state$maxspecID
            nspec <- nrow(island_spec)
            stt <- rbind(stt, c(nspec, timeval))
          }
          stt_table <- rbind(stt_table,
                             c(totaltime - timeval,
                               length(which(island_spec[,4] == "I")),
                               length(which(island_spec[,4] == "A")),
                               length(which(island_spec[,4] == "C"))))
          
        } else {
          ##### After t_hor is reached ####
          # Recalculate t_hor
          t_hor <- get_t_hor(timeval = timeval, totaltime = totaltime, Apars = Apars,
                           ext_multiplier = ext_multiplier,
                           island_ontogeny = island_ontogeny, t_hor = t_hor)
        }
      }
      
      return(stt)
    }
  )
})



# Integrate test by Giovanni Laudanno
# island_area_for_test <- function(timeval, totaltime, Apars, island_function_shape){
#   testit::assert(are_area_params(Apars))
#   Tmax <- Apars$total_island_age # total time A PARS 1
#   Amax <- Apars$max_area # maximum area
#   Topt <- Apars$proportional_peak_t # peak position in %
#   peak <- Apars$peak_sharpness # peakiness - we specify a value of 1 but this is flexible.
#   proptime <- timeval/Tmax
#   # Constant
#   if (is.null(island_function_shape)){
#     return(Apars$max_area)
#   }
#   # Beta function
#   if(island_function_shape == "quadratic") {
#     
#     f <- Topt / (1 - Topt)
#     a <- f * peak/ ( 1 + f)
#     b <- peak / (1 + f)
#     At <- Amax * proptime^a * (1 - proptime)^ b/ ((a / (a + b))^a * (b
#                                                                      / (a + b))^b)
#     return(At)}
#   
#   #Linear decline
#   if(island_function_shape == "linear") {
#     b <- Amax # intercept (peak area)
#     m <- -(b / Topt) # slope
#     At <- m * timeval + b
#     return(At)
#   }
# }
# 
# test_that("Integrate tests by Giovanni", {
#   skip("WIP")
#   expect_silent(
#     
# # Function to describe changes in extinction rate through time. From
# # Valente et al 2014 ProcB
# get_ext_rate_for_test <- function(timeval, totaltime, mu,
#                                   Apars, Epars,
#                                   island_function_shape,
#                                   extcutoff, N,
#                                   K){
#   # Epars[1] and Epars[2] (mu_min, mu_p) must be user specified
#   testit::assert(are_area_params(Apars))
#   if (is.null(island_function_shape)){
#     extrate <- mu * N
#     
#   } else {
#     
#     X <- log(Epars[1] / Epars[2]) / log(0.1)
#     extrate <- Epars[1]/((island_area_for_test(timeval, totaltime, Apars,
#                                                island_function_shape) / Apars$max_area)^X)
#     extrate[which(extrate > extcutoff)] <- extcutoff
#     extrate <- extrate * N
#     return(extrate)
#   }
# }
# MU <- function(s, N = 1000){
#   totaltime <- 10
#   Apars <- create_area_params(1000, 0.2, 1, totaltime * 1.5)
#   Epars <- c(1.7, 20)
#   island_function_shape <- 'quadratic'
#   extcutoff <- 1000
#   # N <- 1000
#   K <- Inf
#   
#   MU.out <- get_ext_rate_for_test(timeval = s,
#                                   totaltime = totaltime,
#                                   mu = 0,
#                                   Apars = Apars,
#                                   Epars = Epars,
#                                   island_function_shape = island_function_shape,
#                                   extcutoff = extcutoff,
#                                   N = N,
#                                   K = K)
#   return(MU.out)
# }
# RHO0 <- function(t, t0 = 0, P0 = 1000){
#   RHO.out <- integrate(f = MU, lower = t0, upper = t, N = P0)$value
#   return(RHO.out)
# }
# RHO  <- function(t, t0 = 0, P0 = 1000){
#   return(Vectorize(RHO0(t = t, t0 = t0, P0 = P0)))
# }
# Pt   <- function(t, t0 = 0, P0 = 1000){
#   out <- P0 * exp(-RHO(t = t, P0 = P0, t0 = t0))
#   return(out)
# }
# Pt(0, 10, 1000)
# 
# 
# PPt <- tt <- seq(0.04,10,0.02);
# for (i in 1:length(tt))
# {
#   PPt[i] <- Pt(t = tt[i])
# }
#   )
#   
# }
