#' Calculates information criterion from DAISIE ML estimates?
#'
#' @param datalist Data object containing information on colonisation and
#' branching times. This object can be generated using the DAISIE_dataprep
#' function, which converts a user-specified data table into a data object, but
#' the object can of course also be entered directly. It is an R list object
#' with the following elements.\cr The first element of the list has two or
#' three components: \cr \cr \code{$island_age} - the island age \cr Then,
#' depending on whether a distinction between types is made, we have:\cr
#' \code{$not_present} - the number of mainland lineages that are not present
#' on the island \cr or:\cr \code{$not_present_type1} - the number of mainland
#' lineages of type 1 that are not present on the island \cr
#' \code{$not_present_type2} - the number of mainland lineages of type 2 that
#' are not present on the island \cr \cr The remaining elements of the list
#' each contains information on a single colonist lineage on the island and has
#' 5 components:\cr \cr \code{$colonist_name} - the name of the species or
#' clade that colonized the island \cr \code{$branching_times} - island age and
#' stem age of the population/species in the case of Non-endemic,
#' Non-endemic_MaxAge and Endemic anagenetic species. For cladogenetic species
#' these should be island age and branching times of the radiation including
#' the stem age of the radiation.\cr \code{$stac} - the status of the colonist
#' \cr \cr * Non_endemic_MaxAge: 1 \cr * Endemic: 2 \cr * Endemic&Non_Endemic:
#' 3 \cr * Non_Endemic: 4 \cr * Endemic_Singleton_MaxAge: 5 \cr *
#' Endemic_Clade_MaxAge: 6 \cr * Endemic&Non_Endemic_Clade_MaxAge: 7 \cr \cr
#' \code{$missing_species} - number of island species that were not sampled for
#' particular clade (only applicable for endemic clades) \cr \code{$type1or2} -
#' whether the colonist belongs to type 1 or type 2 \cr
#' @param initparsopt The initial values of the parameters that must be
#' optimized
#' @param idparsopt The ids of the parameters that must be optimized. The ids
#' are defined as follows: \cr \cr id = 1 corresponds to lambda^c (cladogenesis
#' rate) \cr id = 2 corresponds to mu (extinction rate) \cr id = 3 corresponds
#' to K (clade-level carrying capacity) \cr id = 4 corresponds to gamma
#' (immigration rate) \cr id = 5 corresponds to lambda^a (anagenesis rate) \cr
#' id = 6 corresponds to lambda^c (cladogenesis rate) for an optional subset of
#' the species \cr id = 7 corresponds to mu (extinction rate) for an optional
#' subset of the species\cr id = 8 corresponds to K (clade-level carrying
#' capacity) for an optional subset of the species\cr id = 9 corresponds to
#' gamma (immigration rate) for an optional subset of the species\cr id = 10
#' corresponds to lambda^a (anagenesis rate) for an optional subset of the
#' species\cr id = 11 corresponds to p_f (fraction of mainland species that
#' belongs to the second subset of species
#' @param idparsfix The ids of the parameters that should not be optimized,
#' e.g. c(1,3) if lambda^c and K should not be optimized.
#' @param parsfix The values of the parameters that should not be optimized
#' @param endmc Numeric for how many simulations should run
#' @param res Sets the maximum number of species for which a probability must
#' be computed, must be larger than the size of the largest clade
#' @param cond cond = 0 : conditioning on island age \cr cond = 1 :
#' conditioning on island age and non-extinction of the island biota \cr
#' @param ddmodel Sets the model of diversity-dependence: \cr \cr ddmodel = 0 :
#' no diversity dependence \cr ddmodel = 1 : linear dependence in speciation
#' rate \cr ddmodel = 11: linear dependence in speciation rate and in
#' immigration rate \cr ddmodel = 2 : exponential dependence in speciation
#' rate\cr ddmodel = 21: exponential dependence in speciation rate and in
#' immigration rate\cr
#'
#' @return List of two numerics with WIC and AICb
#' @export
#'
DAISIE_IC <- function(datalist,
                      initparsopt,
                      idparsopt,
                      parsfix,
                      idparsfix,
                      endmc = 1000,
                      res = 100,
                      cond = 0,
                      ddmodel = 0) {
  MLE_obs <- DAISIE_ML(
     datalist = datalist,
     initparsopt = initparsopt,
     idparsopt = idparsopt,
     parsfix = parsfix,
     idparsfix = idparsfix,
     idparsnoshift = 6:10,
     res = res,
     ddmodel = ddmodel,
     cond = cond,
     eqmodel = 0,
     x_E = 0.95,
     x_I = 0.98,
     tol = c(1e-04, 1e-05, 1e-07),
     maxiter = 1000 * round((1.25) ^ length(idparsopt)),
     methode = "lsodes",
     optimmethod = "subplex"
     )
  sims <- DAISIE_sim_constant_rate(
     time = datalist$island_age,
     M = datalist$not_present, #add the number of species that are present
     pars = MLE_obs[1:5],
     replicates = endmc,
     sample_freq = 1,
     plot_sims = FALSE
     )
  MLE <- rep(0, endmc)
  LL <- rep(0, endmc)
  for (mc in 1:endmc) {
    MLE[[mc]] <- DAISIE_ML(
       datalist = sims[[mc]],
       initparsopt = MLE_obs[idparsopt],
       idparsopt = idparsopt,
       parsfix = MLE_obs[parsfix],
       idparsfix = idparsfix,
       idparsnoshift = 6:10,
       res = res,
       ddmodel = ddmodel,
       cond = cond,
       eqmodel = 0,
       x_E = 0.95,
       x_I = 0.98,
       tol = c(1e-04, 1e-05, 1e-07),
       maxiter = 1000 * round((1.25) ^ length(idparsopt)),
       methode = "lsodes",
       optimmethod = "subplex"
       )$loglik
    LL[[mc]] <- DAISIE_loglik_all(
       pars1 = MLE[[mc]][1:5],
       pars2 = c(res, ddmodel, cond, 0),
       datalist = datalist,
       methode = "lsodes"
       )
  }
  WIC <- -2 * MLE_obs$loglik - 2 * mean(LL - MLE)
  AICb <- -2 * MLE_obs$loglik - 4 * mean(LL - MLE_obs)
  out <- list(WIC = WIC, AICb = AICb)
  return(out)
}
